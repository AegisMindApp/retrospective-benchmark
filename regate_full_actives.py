"""Amendment 16 re-gate: apply the registered debias gate at full active count.

The 13 non-Mpro exclusions were all measured at MAX_ACTIVES=30. Mpro passed at
n=30 (AUROC 0.543) and failed at n=85 (0.665), so subsampling made qualification
easier for that target. This re-runs the gate at MAX_ACTIVES=100 for every target
with local panel data, to test whether that direction is general.

No docking. Decoy build + the registered 7-descriptor logistic regression only.
"""
import glob, json, sys, traceback
import decoys

TARGETS = {}
for f in sorted(glob.glob("panel_*.json")):
    if "wiretest" in f or "dryrun" in f:
        continue
    try:
        for tid, comps in (json.load(open(f)).get("targets") or {}).items():
            TARGETS.setdefault(tid, comps)
    except Exception:
        pass

print(f"MAX_ACTIVES = {decoys.MAX_ACTIVES}  (Amendment 16)")
print(f"targets with local panel data: {len(TARGETS)}\n", flush=True)

out = {}
for tid, comps in TARGETS.items():
    act = [c for c in comps if c["label"] == 1]
    print(f"--- {tid}: {len(act)} actives ---", flush=True)
    try:
        d = decoys.build_decoys(
            [{"molecule_chembl_id": a["molecule_chembl_id"], "smiles": a["smiles"]} for a in act], tid)
        rec = {k: d.get(k) for k in ("status", "n_actives", "n_decoys", "debias_auroc", "scaffold_count")}
    except Exception as e:
        rec = {"status": f"error: {e}"}
        traceback.print_exc()
    out[tid] = rec
    au = rec.get("debias_auroc")
    print(f"    status={rec['status']} n_act={rec.get('n_actives')} "
          f"auroc={au if au is None else round(au,4)}\n", flush=True)
    json.dump(out, open("regate_full_actives_report.json", "w"), indent=2)

print("=" * 62)
print(f"{'target':<16}{'n_act':>6}{'gate AUROC':>12}  status")
for tid, r in out.items():
    au = r.get("debias_auroc")
    print(f"{tid:<16}{str(r.get('n_actives','-')):>6}{('-' if au is None else f'{au:.4f}'):>12}  {r['status']}")
