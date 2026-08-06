"""Merge Kaggle shard results into the local cache — but only if they check out.

Three gates, all of which must pass before a single remote score is accepted:

1. **Protocol identity.** The shard's manifest records md5s of the Vina binary and
   the receptor, plus exhaustiveness / box / num_modes. These must equal the local
   protocol sidecar. md5s rather than version strings: a version is a claim.

2. **Compound-list identity.** `list_md5` must match the frozen list, so both sides
   were working through the same 2,526 compounds in the same order.

3. **Cross-environment calibration.** Each shard re-docks a handful of compounds
   already scored locally. Vina has no fixed seed, so these will not agree exactly —
   but they should agree within --tolerance kcal/mol. A wider spread means the
   environments differ by more than the RNG and the merge is refused.

Gate 3 is the one that would actually catch a silently different environment, which
is the failure this whole design exists to prevent.

Usage:
    python merge_shards.py --tid CHEMBL612545 shard0_cache.json shard1_cache.json
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
sys.path.insert(0, PARENT)
import vina_cache as vc            # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tid", required=True)
    ap.add_argument("shards", nargs="+", help="shard cache JSONs returned from Kaggle")
    ap.add_argument("--manifest-dir", default=None,
                    help="dir holding shard<N>/manifest.json (default: bundle_<tid>)")
    ap.add_argument("--tolerance", type=float, default=0.3,
                    help="max |local - remote| kcal/mol on calibration compounds")
    ap.add_argument("--apply", action="store_true",
                    help="actually write; default is a dry run")
    a = ap.parse_args()

    local_p = os.path.join(PARENT, f"vina_cache_{a.tid}.json")
    local = vc.load(local_p)
    local_proto = vc.load(vc.protocol_path(local_p))
    if not local_proto:
        raise SystemExit(f"local cache has no protocol sidecar — refusing to merge")

    mdir = a.manifest_dir or os.path.join(HERE, f"bundle_{a.tid}")
    accepted, rejected = {}, []

    for sp in a.shards:
        sc = vc.load(sp)
        if not sc:
            rejected.append((sp, "empty or unreadable")); continue

        # locate the manifest that describes this shard
        man = None
        for d in sorted(os.listdir(mdir)) if os.path.isdir(mdir) else []:
            mp = os.path.join(mdir, d, "manifest.json")
            if os.path.exists(mp):
                m = json.load(open(mp))
                if f"shard{m['shard']}" in os.path.basename(sp):
                    man = m; break
        if man is None:
            rejected.append((sp, "no matching manifest")); continue

        # gate 1 — protocol identity (compare only keys both sides record)
        rp = man["protocol"]
        shared = set(rp) & set(local_proto)
        diff = {k: (local_proto[k], rp[k]) for k in shared if local_proto[k] != rp[k]}
        if diff:
            rejected.append((sp, f"protocol mismatch: {diff}")); continue

        # gate 2 — same compound list
        frozen = json.load(open(os.path.join(HERE, f"frozen_compounds_{a.tid}.json")))
        if man.get("list_md5") and man["list_md5"] != frozen.get("list_md5"):
            rejected.append((sp, "compound list md5 differs")); continue

        # gate 3 — cross-environment calibration
        deltas = []
        for name in man.get("calibration", []):
            k = f"bench_{name}"
            if k in local and k in sc:
                deltas.append((name, local[k], sc[k], abs(local[k] - sc[k])))
        if not deltas:
            rejected.append((sp, "no calibration overlap — cannot verify environment"))
            continue
        worst = max(d[3] for d in deltas)
        print(f"\n{os.path.basename(sp)} calibration ({len(deltas)} compounds):")
        for n, lo, re_, d in deltas:
            print(f"  {n:<18} local {lo:>7.2f}  remote {re_:>7.2f}  |d| {d:.2f}"
                  f"{'  <-- OVER' if d > a.tolerance else ''}")
        print(f"  worst {worst:.2f}, mean {statistics.mean(d[3] for d in deltas):.2f} "
              f"(tolerance {a.tolerance})")
        if worst > a.tolerance:
            rejected.append((sp, f"calibration worst |d|={worst:.2f} > {a.tolerance}"))
            continue

        # accept everything except the calibration replants (already scored locally)
        calib_keys = {f"bench_{n}" for n in man.get("calibration", [])}
        new = {k: v for k, v in sc.items() if k not in calib_keys and k not in local}
        accepted.update(new)
        print(f"  ACCEPTED {len(new)} new results")

    print(f"\n{'='*60}")
    for sp, why in rejected:
        print(f"REJECTED {os.path.basename(sp)}: {why}")
    print(f"accepted {len(accepted)} new results; local had {len(local)}")

    if not a.apply:
        print("\nDRY RUN — re-run with --apply to write")
        return 1 if rejected else 0

    merged = {**local, **accepted}
    vc.save(merged, local_p)
    print(f"local cache now {len(merged)} -> {local_p}")
    return 1 if rejected else 0


if __name__ == "__main__":
    raise SystemExit(main())
