"""Is the union of local + Kaggle results complete? Should the local run be stopped?

The local PD-L1 process runs pre-fix code: it will grind through all 2,526 compounds
regardless of what Kaggle returns, including ones Kaggle has already done. That is
wasted work, not wrong work — but it needs an explicit stopping rule, because the
run is only finished when the UNION covers the frozen list, not when either side
individually does.

Local works the head of the list, Kaggle the tail. They meet in the middle.

Exit codes:  0 = complete, stop local   |   2 = still gaps

Usage:
    python check_coverage.py --tid CHEMBL612545 [shard caches...]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
sys.path.insert(0, PARENT)
import vina_cache as vc            # noqa: E402


def strip(keys) -> set:
    return {k[len("bench_"):] if k.startswith("bench_") else k for k in keys}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tid", required=True)
    ap.add_argument("shards", nargs="*", help="shard caches returned from Kaggle")
    ap.add_argument("--log", default=None, help="live log to scrape as well")
    a = ap.parse_args()

    frozen = json.load(open(os.path.join(HERE, f"frozen_compounds_{a.tid}.json")))
    want = [c["name"] for c in frozen["compounds"]]

    local = vc.load(os.path.join(PARENT, f"vina_cache_{a.tid}.json"))
    if a.log and os.path.exists(a.log):
        local = {**vc.recover_from_log(a.log), **local}
    have_local = strip(local)

    have_remote: set = set()
    for sp in a.shards:
        have_remote |= strip(vc.load(sp))

    union = have_local | have_remote
    missing = [n for n in want if n not in union]

    print(f"frozen list      {len(want)}")
    print(f"local            {len(have_local)}")
    print(f"kaggle shards    {len(have_remote)}")
    print(f"overlap (waste)  {len(have_local & have_remote)}")
    print(f"union            {len(union)}  ({100*len(union)/len(want):.1f}%)")
    print(f"missing          {len(missing)}")

    if missing:
        # where the gap sits tells you whether the two ends have met
        idx = [i for i, n in enumerate(want) if n in set(missing)]
        print(f"gap spans frozen idx {idx[0]}..{idx[-1]}")
        print("\nNOT complete — leave the local run going.")
        return 2

    print("\nCOMPLETE — the union covers the frozen list. Stop the local run:")
    print("  pkill -f panel_pdl1_only     # then merge_shards.py --apply")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
