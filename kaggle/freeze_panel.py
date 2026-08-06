"""Freeze the exact compound list a docking run works through.

Why this is the first step
--------------------------
`run_pilot_fast.py` builds its compound list at runtime: actives are subsampled
from the panel (seeded), but DECOYS are pulled live from the ChEMBL REST API.
Amendment 11 made the API responses deterministically *ordered*, which is not the
same as making them *stable* — ChEMBL is a living database, so a query issued from
Kaggle next week can legitimately return different molecules.

That makes regenerating the list on the remote side unsafe: two environments would
dock two different 2,526-compound sets and the shards would not correspond. So the
list is frozen to disk here, once, and shipped.

The list is written in the EXACT order the runner produces it (actives + decoys +
real inactives, then `random.Random(0).shuffle`), so a frozen index means the same
compound everywhere.

Verification against a live run
-------------------------------
For a target already part-way through docking, every compound name in the existing
cache must appear in the regenerated list. If any is missing, ChEMBL has moved
underneath us and the frozen list does NOT describe the running job — in which case
that target cannot be sharded and the script exits non-zero.

Usage:
    python freeze_panel.py panel_pdl1_only.json CHEMBL612545 \
        [--verify-cache ../vina_cache_CHEMBL612545.json]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
sys.path.insert(0, PARENT)

import decoys                      # noqa: E402
from run_pilot_fast import SHUFFLE_SEED   # noqa: E402  (single source of truth)


def freeze(panel_path: str, tid: str) -> list[dict]:
    """Reproduce run_pilot_fast's compound list for one target, in order."""
    panel = json.load(open(panel_path))
    comps = panel["targets"][tid]
    actives = [c for c in comps if c["label"] == 1]
    real_inact = [c for c in comps if c["label"] == 0]

    print(f"building decoys ({len(actives)} actives -> max {decoys.MAX_ACTIVES}) ...",
          flush=True)
    dres = decoys.build_decoys(
        [{"molecule_chembl_id": a["molecule_chembl_id"], "smiles": a["smiles"]}
         for a in actives], tid)
    if dres["status"] != "ok":
        raise SystemExit(f"decoy build failed: {dres['status']}")
    print(f"decoys: n={dres.get('n_decoys')} auroc={dres.get('debias_auroc')}", flush=True)

    scored = (
        [{"smiles": a["smiles"], "molecule_chembl_id": a["molecule_chembl_id"], "label": 1}
         for a in dres["actives"]]
        + [{"smiles": d["smiles"], "molecule_chembl_id": d["molecule_chembl_id"], "label": 0}
           for d in dres["decoys"]]
        + [{"smiles": c["smiles"], "molecule_chembl_id": c["molecule_chembl_id"], "label": 0}
           for c in real_inact]
    )
    random.Random(SHUFFLE_SEED).shuffle(scored)      # identical to the runner
    return scored


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("panel")
    ap.add_argument("tid")
    ap.add_argument("--verify-cache", default=None,
                    help="cache from a live run; every name in it must be present")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    scored = freeze(a.panel, a.tid)
    names = [c["molecule_chembl_id"] for c in scored]
    print(f"frozen list: {len(names)} compounds "
          f"({sum(c['label'] for c in scored)} actives)")

    if a.verify_cache:
        cache = json.load(open(a.verify_cache)) if os.path.exists(a.verify_cache) else {}
        docked = {k[len("bench_"):] if k.startswith("bench_") else k for k in cache}
        missing = sorted(docked - set(names))
        print(f"\nverification against {len(docked)} already-docked compounds:")
        if missing:
            print(f"  MISSING {len(missing)}: {missing[:10]}"
                  f"{' ...' if len(missing) > 10 else ''}")
            print("  ChEMBL has changed since the live run started. The frozen list does\n"
                  "  NOT describe the running job -- this target cannot be sharded.")
            return 1
        print(f"  all {len(docked)} present -> ChEMBL is stable, list matches the live run")

    out = a.out or os.path.join(HERE, f"frozen_compounds_{a.tid}.json")
    payload = {"target": a.tid, "shuffle_seed": SHUFFLE_SEED,
               "max_actives": decoys.MAX_ACTIVES, "n": len(scored),
               "compounds": [{"name": c["molecule_chembl_id"], "smiles": c["smiles"],
                              "label": c["label"]} for c in scored]}
    body = json.dumps(payload, indent=1)
    json.dump({**payload,
               "list_md5": hashlib.md5(
                   json.dumps(payload["compounds"], sort_keys=True).encode()).hexdigest()},
              open(out, "w"), indent=1)
    print(f"\nwrote {out}  ({len(body)//1024} KB)")
    print(f"list_md5 = {hashlib.md5(json.dumps(payload['compounds'], sort_keys=True).encode()).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
