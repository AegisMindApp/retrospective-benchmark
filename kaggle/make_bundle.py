"""Build a self-contained Kaggle shard bundle.

Design rule: NOTHING that affects the docking protocol is reproduced remotely.

Ligand PDBQTs are generated HERE, with this machine's RDKit 2026.03.3 and Meeko
0.7.1, and shipped. The receptor and the Vina 1.2.7 binary are shipped byte-for-byte.
Kaggle therefore runs exactly one command per compound:

    vina --receptor R.pdbqt --ligand L.pdbqt --center/--size ... --exhaustiveness 32

with no RDKit, no Meeko, and no version negotiation on the remote side. This is
deliberate: Meeko and RDKit build the torsion tree, an RDKit bump alone has already
moved gate AUROC in this project (0.585 -> 0.598), and Amendment 19 exists because
ligand prep silently differed from what was assumed. Pinning versions would be a
claim; shipping the artifacts is a fact.

The bundle manifest records md5s, not version strings, for the same reason.

Shards are taken from the TAIL of the remaining compounds because the local process
works from the head — they meet in the middle. See check_coverage.py for the
stopping rule.

Usage:
    python make_bundle.py --tid CHEMBL612545 --shards 2 --count 1200
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
sys.path.insert(0, PARENT)
sys.path.insert(0, os.path.join(os.path.dirname(PARENT), "amr_glass", "docking"))

VINA_BIN = os.environ.get("VINA_BIN", os.path.join(PARENT, "bin", "vina"))


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tid", required=True)
    ap.add_argument("--shards", type=int, default=2, help="Kaggle allows 2 concurrent")
    ap.add_argument("--count", type=int, default=None,
                    help="how many of the remaining to hand to Kaggle (from the tail); "
                         "default = all remaining")
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--num-modes", type=int, default=5)
    ap.add_argument("--calibration", type=int, default=5,
                    help="already-docked compounds replanted in each shard so the "
                         "remote environment can be checked against the local one")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--also-assigned", nargs="*", default=[],
                    help="manifest.json files whose compounds are already handed to some "
                         "other shard. Without this, 'remaining' is computed from the LOCAL "
                         "cache alone, so a second bundle re-hands work the first already "
                         "holds — which is how 499 compounds of this target ended up in no "
                         "bundle at all while others sat in two.")
    a = ap.parse_args()

    frozen_p = os.path.join(HERE, f"frozen_compounds_{a.tid}.json")
    if not os.path.exists(frozen_p):
        raise SystemExit(f"missing {frozen_p} — run freeze_panel.py first")
    frozen = json.load(open(frozen_p))
    comps = frozen["compounds"]

    cache_p = os.path.join(PARENT, f"vina_cache_{a.tid}.json")
    cache = json.load(open(cache_p)) if os.path.exists(cache_p) else {}
    done = {k[len("bench_"):] if k.startswith("bench_") else k for k in cache}

    assigned = set()
    for mp in a.also_assigned:
        for c in json.load(open(mp))["compounds"]:
            assigned.add(c["name"])
    if assigned:
        print(f"{len(assigned)} compounds already assigned to other shards — excluded")

    remaining = [(i, c) for i, c in enumerate(comps)
                 if c["name"] not in done and c["name"] not in assigned]
    print(f"{len(comps)} total, {len(done)} docked, {len(remaining)} remaining")

    take = a.count or len(remaining)
    tail = remaining[-take:]          # from the tail; local works the head
    print(f"handing {len(tail)} to Kaggle (frozen idx {tail[0][0]}..{tail[-1][0]})")

    # Take the receptor path from receptors.json, NOT from the {tid}_receptor.pdbqt
    # naming convention. On PD-L1 the convention path held a receptor that had FAILED
    # its redocking gate (5J89) while receptors.json pointed at the one that passed
    # (5J8O); the bundle would have docked 455 ligands against the failed structure
    # using the passing structure's box -- a different crystal frame, so the box sat
    # in empty space -- and returned a perfectly well-formed AUROC. Only the receptor
    # md5 check caught it.
    receptor = os.path.join(PARENT, "receptors", f"{a.tid}_receptor.pdbqt")
    recjson = json.load(open(os.path.join(PARENT, "receptors.json")))[a.tid] \
        if os.path.exists(os.path.join(PARENT, "receptors.json")) else None
    if recjson and recjson.get("receptor_pdbqt"):
        receptor = recjson["receptor_pdbqt"]
    if recjson is None:
        raise SystemExit("receptors.json not found — need box_center/box_size")

    outroot = a.outdir or os.path.join(HERE, f"bundle_{a.tid}")
    if os.path.exists(outroot):
        shutil.rmtree(outroot)

    # Ligand prep, locally, with this machine's toolchain.
    import run_vina_docking as rvd            # noqa: E402
    print(f"\npreparing ligand PDBQTs locally (RDKit+Meeko stay on THIS machine) ...")

    # Calibration controls: compounds ALREADY docked locally, replanted in every
    # shard. Vina has no --seed, so these will not match to 2dp — but a local and a
    # remote score for the same ligand should agree within ~0.3 kcal/mol. A wider
    # gap means the environments differ by more than the RNG, and it is much better
    # to learn that from 5 compounds than from 1,100.
    calib = [(i, c) for i, c in enumerate(comps) if c["name"] in done][:a.calibration]
    if len(calib) < a.calibration:
        print(f"WARNING: only {len(calib)} calibration compounds available")
    print(f"calibration controls per shard: {[c['name'] for _, c in calib]}")

    # Calibration FIRST, not last. Appending them meant they would be docked after
    # ~60h of compute, which defeats the entire purpose — the point is to learn the
    # environments differ from 5 compounds rather than from 1,100.
    per = [calib + tail[i::a.shards] for i in range(a.shards)]  # stride: even difficulty
    for s, items in enumerate(per):
        sdir = os.path.join(outroot, f"shard{s}")
        os.makedirs(os.path.join(sdir, "ligands"), exist_ok=True)
        ok, failed = [], []
        for n, (idx, c) in enumerate(items):
            lp = os.path.join(sdir, "ligands", f"bench_{c['name']}.pdbqt")
            try:
                if rvd.smiles_to_pdbqt(c["smiles"], f"bench_{c['name']}", lp):
                    ok.append({"idx": idx, "name": c["name"], "label": c["label"]})
                else:
                    failed.append(c["name"])
            except Exception as e:
                failed.append(f"{c['name']}({e})")
            if (n + 1) % 100 == 0:
                print(f"  shard{s}: {n+1}/{len(items)} prepared", flush=True)

        shutil.copy2(VINA_BIN, os.path.join(sdir, "vina"))
        os.chmod(os.path.join(sdir, "vina"), 0o755)
        shutil.copy2(receptor, os.path.join(sdir, "receptor.pdbqt"))

        manifest = {
            "target": a.tid, "shard": s, "n_shards": a.shards,
            "compounds": ok, "prep_failed": failed,
            "calibration": [c["name"] for _, c in calib],
            "list_md5": frozen.get("list_md5"),
            # md5s, not version strings: a version is a claim, an md5 is a fact
            "protocol": {
                "exhaustiveness": a.exhaustiveness,
                "num_modes": a.num_modes,
                "box_center": [round(float(v), 2) for v in recjson["box_center"]],
                "box_size": [round(float(v), 2) for v in recjson["box_size"]],
                "flexible_ligand": True,
                "vina_md5": md5(VINA_BIN),
                "receptor_md5": md5(receptor),
            },
        }
        json.dump(manifest, open(os.path.join(sdir, "manifest.json"), "w"), indent=1)
        shutil.copy2(os.path.join(HERE, "kaggle_shard_worker.py"),
                     os.path.join(sdir, "kaggle_shard_worker.py"))
        print(f"shard{s}: {len(ok)} ligands ready, {len(failed)} prep-failed -> {sdir}")

    print(f"\nbundle at {outroot}")
    print("upload each shard<N>/ as a Kaggle Dataset, then run kaggle_shard_worker.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
