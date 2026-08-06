#!/usr/bin/env python3
"""
run_pilot.py — orchestrator. Enforces the blind→unblind order.

Sequence (each step is a hard checkpoint; see PREREGISTRATION.md):
  0. PREREGISTRATION.md is committed (freeze) — done out-of-band, in git.
  1. build_panel() froze panel_frozen.json (data-availability only; platform never consulted).
  2. This script produces ALL method rankings BLIND (labels never touched during scoring),
     writes predictions_<method>.json.
  3. Only then does it unblind: join to labels, compute metrics + the paired panel test,
     write report.json.

Any method whose inputs are missing (e.g. no receptor for the Vina/platform paths, or no
API keys) is SKIPPED and reported as skipped — never silently substituted, never faked.
"""
from __future__ import annotations

import asyncio
import json
import os
import random
from typing import Dict, List

import numpy as np

SHUFFLE_SEED = 0   # deterministic: ties must not favour the actives-first assembly order

from metrics import score_ranking, bootstrap_ci, paired_panel_test, EF_FRACTIONS
from baselines import chance_scores, fp_similarity_scores
import platform_rank
import decoys

PRIMARY_METRIC = "EF@1%"          # frozen primary endpoint (see PREREGISTRATION.md)


def _load_panel(path="panel_frozen.json") -> dict:
    with open(path) as fh:
        return json.load(fh)


def _load_receptors(path="receptors.json") -> dict:
    """Optional: {target_chembl_id: {receptor_pdbqt, box_center:[x,y,z], box_size:[x,y,z],
    target_name, context}}. Absent → Vina/platform paths are skipped, chance/fp still run."""
    return json.load(open(path)) if os.path.exists(path) else {}


async def run(panel_path="panel_frozen.json"):
    panel = _load_panel(panel_path)["targets"]
    receptors = _load_receptors()
    per_target_metrics: Dict[str, Dict[str, dict]] = {}   # target -> method -> metrics
    decoy_report: Dict[str, dict] = {}                    # target -> decoy/debias summary
    skipped: List[str] = []

    for tid, comps in panel.items():
        actives = [c for c in comps if c["label"] == 1]
        real_inactives = [c for c in comps if c["label"] == 0]   # keep genuine post-cutoff inactives

        # Amendment 1: supply negatives as property-matched debiased decoys.
        dres = decoys.build_decoys(
            [{"molecule_chembl_id": a["molecule_chembl_id"], "smiles": a["smiles"]} for a in actives],
            tid)
        decoy_report[tid] = {k: dres.get(k) for k in
                             ("status", "n_actives", "n_decoys", "debias_auroc", "scaffold_count")}
        if dres["status"] != "ok":
            # property-separable (or too few decoys) → excluded, never shipped
            skipped.append(f"{tid}:decoy-excluded ({dres['status']}, "
                           f"debias_auroc={dres.get('debias_auroc')})")
            continue

        scored = ([{"smiles": a["smiles"], "molecule_chembl_id": a["molecule_chembl_id"], "label": 1}
                   for a in dres["actives"]]                    # possibly subsampled to MAX_ACTIVES
                  + [{"smiles": d["smiles"], "molecule_chembl_id": d["molecule_chembl_id"], "label": 0}
                     for d in dres["decoys"]]
                  + [{"smiles": c["smiles"], "molecule_chembl_id": c["molecule_chembl_id"], "label": 0}
                     for c in real_inactives])
        random.Random(SHUFFLE_SEED).shuffle(scored)   # break actives-first tie bias
        smiles = [c["smiles"] for c in scored]
        labels = [c["label"] for c in scored]
        names = [c["molecule_chembl_id"] for c in scored]
        rec = receptors.get(tid)
        methods_scores: Dict[str, List[float]] = {}

        # --- BLIND scoring (labels not consulted) ---
        methods_scores["chance"] = chance_scores(len(smiles), seed=0)
        try:
            methods_scores["fp_similarity"] = fp_similarity_scores(smiles, tid)
        except Exception as e:
            skipped.append(f"{tid}:fp_similarity ({type(e).__name__})")

        if rec:
            try:
                methods_scores["vina"] = await platform_rank.platform_scores(
                    smiles, rec["target_name"], rec.get("context", ""),
                    rec["receptor_pdbqt"], tuple(rec["box_center"]), tuple(rec["box_size"]),
                    names=names, use_ensemble=False)          # docking-only
                methods_scores["platform"] = await platform_rank.platform_scores(
                    smiles, rec["target_name"], rec.get("context", ""),
                    rec["receptor_pdbqt"], tuple(rec["box_center"]), tuple(rec["box_size"]),
                    names=names, use_ensemble=True)            # docking + ensemble
            except Exception as e:
                skipped.append(f"{tid}:vina/platform ({type(e).__name__}: {str(e)[:60]})")
        else:
            skipped.append(f"{tid}:vina/platform (no receptor)")

        # persist blind predictions before unblinding
        json.dump({m: s for m, s in methods_scores.items()},
                  open(f"predictions_{tid}.json", "w"))

        # --- UNBLIND: metrics ---
        per_target_metrics[tid] = {}
        for m, s in methods_scores.items():
            met = score_ranking(s, labels)
            p, lo, hi = bootstrap_ci(s, labels, PRIMARY_METRIC)
            met[f"{PRIMARY_METRIC}_ci"] = [p, lo, hi]
            per_target_metrics[tid][m] = met

    # --- Paired panel test: the PRIMARY claim ---
    report = {"primary_metric": PRIMARY_METRIC, "per_target": per_target_metrics,
              "decoy_report": decoy_report, "skipped": skipped, "panel_tests": {}}
    def col(method):
        return {t: per_target_metrics[t][method][PRIMARY_METRIC]
                for t in per_target_metrics if method in per_target_metrics[t]}
    if any("platform" in per_target_metrics[t] for t in per_target_metrics):
        report["panel_tests"]["platform_vs_vina"] = paired_panel_test(col("platform"), col("vina"))
        report["panel_tests"]["platform_vs_chance"] = paired_panel_test(col("platform"), col("chance"))
        report["panel_tests"]["vina_vs_chance"] = paired_panel_test(col("vina"), col("chance"))

    json.dump(report, open("report.json", "w"), indent=2)
    print(json.dumps({"primary_metric": PRIMARY_METRIC,
                      "targets_scored": len(per_target_metrics),
                      "skipped": skipped,
                      "panel_tests": report["panel_tests"]}, indent=2))
    return report


if __name__ == "__main__":
    asyncio.run(run())
