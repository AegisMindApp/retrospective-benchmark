"""Kaggle-side docking worker. Runs one shard of a frozen compound list.

Deliberately has NO scientific dependencies — no RDKit, no Meeko, no OpenBabel.
Ligand PDBQTs, the receptor and the Vina binary all arrive pre-built in the shard
bundle, so this script only shells out to Vina and parses its output. That is what
makes a Kaggle result mergeable with a local one: there is no remote toolchain that
could differ.

Session handling
----------------
Kaggle CPU sessions are capped at 12h and the filesystem is ephemeral. The cache is
written after EVERY compound to /kaggle/working, and the worker stops cleanly at
--budget-hours (default 11.3) so the notebook can commit its output before the
session is killed. Re-running with the same output dataset attached resumes.

Usage (in a Kaggle notebook cell):
    !python kaggle_shard_worker.py --bundle /kaggle/input/<dataset> --workers 4
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

SCORE = re.compile(r"^\s*1\s+(-?\d+\.\d+)", re.M)   # first ranked mode


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def dock_one(a):
    vina, receptor, lig, name, ctr, size, exh, modes, outdir = a
    out = os.path.join(outdir, f"{name}_out.pdbqt")
    cmd = [vina, "--receptor", receptor, "--ligand", lig,
           "--center_x", str(ctr[0]), "--center_y", str(ctr[1]), "--center_z", str(ctr[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(exh), "--num_modes", str(modes), "--out", out]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        m = SCORE.search(r.stdout)
        return (name, float(m.group(1))) if m else (name, None)
    except Exception:
        return (name, None)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True, help="/kaggle/input/<dataset-name>")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--budget-hours", type=float, default=11.3,
                    help="stop cleanly before Kaggle's 12h kill")
    ap.add_argument("--workdir", default="/kaggle/working")
    a = ap.parse_args()

    # kagglehub zips any upload with >50 files, and Kaggle may present it either
    # extracted at the dataset root, nested in a subdirectory, or still zipped.
    # Locate the manifest rather than assuming a layout.
    bundle = a.bundle
    if not os.path.exists(os.path.join(bundle, "manifest.json")):
        hits = glob.glob(os.path.join(bundle, "**", "manifest.json"), recursive=True)
        if hits:
            bundle = os.path.dirname(hits[0])
        else:
            zips = glob.glob(os.path.join(bundle, "**", "*.zip"), recursive=True)
            if not zips:
                raise SystemExit(f"no manifest.json and no zip under {bundle}")
            import zipfile
            bundle = os.path.join(a.workdir, "bundle")
            os.makedirs(bundle, exist_ok=True)
            with zipfile.ZipFile(zips[0]) as z:
                z.extractall(bundle)
            hits = glob.glob(os.path.join(bundle, "**", "manifest.json"), recursive=True)
            if not hits:
                raise SystemExit(f"manifest.json missing from {zips[0]}")
            bundle = os.path.dirname(hits[0])
        print(f"bundle resolved to {bundle}")
    a.bundle = bundle

    man = json.load(open(os.path.join(a.bundle, "manifest.json")))
    p = man["protocol"]
    tid, shard = man["target"], man["shard"]

    # Kaggle input is read-only and cannot hold the +x bit: stage the binary.
    vina = os.path.join(a.workdir, "vina")
    shutil.copy2(os.path.join(a.bundle, "vina"), vina)
    os.chmod(vina, 0o755)
    receptor = os.path.join(a.bundle, "receptor.pdbqt")

    # The protocol is only meaningful if the artifacts are the ones it describes.
    if md5(vina) != p["vina_md5"]:
        raise SystemExit(f"vina md5 mismatch: {md5(vina)} != {p['vina_md5']}")
    if md5(receptor) != p["receptor_md5"]:
        raise SystemExit(f"receptor md5 mismatch")
    print(f"artifacts verified: vina {p['vina_md5'][:8]} receptor {p['receptor_md5'][:8]}")

    cache_p = os.path.join(a.workdir, f"vina_cache_{tid}_shard{shard}.json")

    # Resume across sessions. Kaggle caps a session at 12h, so a shard spans
    # several runs; each one carries the previous cache back in as an attached
    # dataset. Search ALL of /kaggle/input rather than just the bundle, because
    # the progress cache arrives as its own dataset, and MERGE every match — with
    # more than one prior dataset attached, taking the first found would silently
    # discard a session's work.
    cache = json.load(open(cache_p)) if os.path.exists(cache_p) else {}
    found = sorted(set(glob.glob(f"/kaggle/input/**/vina_cache_{tid}_shard{shard}.json",
                                 recursive=True))
                   | set(glob.glob(os.path.join(a.bundle, f"vina_cache_{tid}_shard{shard}.json"))))
    for prev in found:
        try:
            got = json.load(open(prev))
            new = {k: v for k, v in got.items() if k not in cache}
            cache.update(got)
            print(f"  [resume] {prev}: {len(got)} results (+{len(new)} new)")
        except Exception as e:
            print(f"  [resume] WARNING: could not read {prev}: {e}")
    if cache:
        json.dump(cache, open(cache_p, "w"))
    print(f"cache holds {len(cache)} results from {len(found)} prior session(s)")

    outdir = os.path.join(a.workdir, "poses")
    os.makedirs(outdir, exist_ok=True)

    todo = [c for c in man["compounds"] if f"bench_{c['name']}" not in cache]
    print(f"shard {shard}/{man['n_shards']}: {len(man['compounds'])} assigned, "
          f"{len(todo)} to do, exh={p['exhaustiveness']}")

    args = [(vina, receptor,
             os.path.join(a.bundle, "ligands", f"bench_{c['name']}.pdbqt"),
             f"bench_{c['name']}", p["box_center"], p["box_size"],
             p["exhaustiveness"], p["num_modes"], outdir) for c in todo]

    t0, done, stopped = time.time(), 0, False
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(dock_one, x): x[3] for x in args}
        for fut in as_completed(futs):
            name, score = fut.result()
            if score is not None:
                cache[name] = score
                json.dump(cache, open(cache_p + ".tmp", "w"))
                os.replace(cache_p + ".tmp", cache_p)     # atomic, survives a kill
            done += 1
            el = (time.time() - t0) / 3600
            if done % 10 == 0:
                rate = done / el if el else 0
                print(f"  {done}/{len(args)}  {rate:.1f}/h  {el:.2f}h elapsed", flush=True)
            if el > a.budget_hours:
                print(f"\nbudget {a.budget_hours}h reached — stopping cleanly. "
                      f"{len(cache)} results saved.")
                stopped = True
                for f in futs:
                    f.cancel()
                break

    json.dump({"target": tid, "shard": shard, "n_results": len(cache),
               "protocol": p, "list_md5": man.get("list_md5"),
               "complete": not stopped and len(cache) >= len(man["compounds"])},
              open(os.path.join(a.workdir, f"shard{shard}_report.json"), "w"), indent=1)
    print(f"\n{len(cache)} results in {cache_p}")
    print("Save/commit this notebook so /kaggle/working persists, then attach it as "
          "an input dataset next session to resume.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
