#!/usr/bin/env python3
"""
metrics.py — enrichment metrics for the blind retrospective benchmark.

Uses RDKit's validated scoring implementation (rdkit.ML.Scoring) for EF / BEDROC /
AUC rather than hand-rolling them, then adds bootstrap CIs and a paired
cross-target significance test (the primary claim is about the *distribution*
across the target panel, not any single hero target).

A "ranking" here is a list of (score, label) with label ∈ {1 active, 0 inactive},
where a HIGHER score = the method's stronger belief the compound is active.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple, Dict
import numpy as np
from rdkit.ML.Scoring.Scoring import CalcBEDROC, CalcEnrichment, CalcAUC

BEDROC_ALPHA = 20.0            # ~80% of the weight in the first ~8% — early recognition
EF_FRACTIONS = (0.01, 0.05)   # EF@1%, EF@5%


def _sorted_labels(scores: Sequence[float], labels: Sequence[int]) -> List[List[int]]:
    """RDKit scoring wants rows sorted by score DESC; each row exposes the active
    flag at a column index. We hand it single-column rows of [label]."""
    order = np.argsort(-np.asarray(scores, dtype=float))
    return [[int(labels[i])] for i in order]


def score_ranking(scores: Sequence[float], labels: Sequence[int]) -> Dict[str, float]:
    """All primary/secondary metrics for one target's ranking."""
    rows = _sorted_labels(scores, labels)
    out: Dict[str, float] = {}
    ef = CalcEnrichment(rows, 0, list(EF_FRACTIONS))
    for frac, val in zip(EF_FRACTIONS, ef):
        out[f"EF@{int(frac*100)}%"] = float(val)
    out["BEDROC"] = float(CalcBEDROC(rows, 0, BEDROC_ALPHA))
    out["AUROC"] = float(CalcAUC(rows, 0))
    out["n"] = len(labels)
    out["n_active"] = int(np.sum(labels))
    return out


def bootstrap_ci(scores: Sequence[float], labels: Sequence[int], metric: str,
                 n_boot: int = 2000, seed: int = 0, ci: float = 0.95
                 ) -> Tuple[float, float, float]:
    """Point estimate + (lo, hi) percentile CI for one metric on one target,
    resampling compounds with replacement. `seed` is explicit (no wall-clock RNG)
    so the interval is reproducible and part of the frozen protocol."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=int)
    point = score_ranking(scores, labels)[metric]
    rng = np.random.default_rng(seed)
    N = len(labels)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, N, N)
        if labels[idx].sum() == 0 or labels[idx].sum() == len(idx):
            continue  # a degenerate resample (all one class) — skip
        vals.append(score_ranking(scores[idx], labels[idx])[metric])
    lo, hi = np.percentile(vals, [(1-ci)/2*100, (1+ci)/2*100]) if vals else (np.nan, np.nan)
    return float(point), float(lo), float(hi)


def paired_panel_test(per_target_method: Dict[str, float],
                      per_target_baseline: Dict[str, float]) -> Dict[str, float]:
    """The PRIMARY test: does the method beat the baseline across the panel?
    Wilcoxon signed-rank on paired per-target metric values (platform vs Vina).
    Non-parametric, no normality assumption, robust to a small panel."""
    from scipy.stats import wilcoxon
    keys = sorted(set(per_target_method) & set(per_target_baseline))
    a = np.array([per_target_method[k] for k in keys], dtype=float)
    b = np.array([per_target_baseline[k] for k in keys], dtype=float)
    diff = a - b
    res = {"n_targets": len(keys), "mean_delta": float(np.mean(diff)),
           "median_delta": float(np.median(diff)),
           "wins": int(np.sum(diff > 0)), "losses": int(np.sum(diff < 0))}
    nonzero = diff[diff != 0]
    if len(nonzero) >= 1:
        try:
            stat, p = wilcoxon(a[diff != 0], b[diff != 0])
            res["wilcoxon_stat"], res["p_value"] = float(stat), float(p)
        except ValueError:
            res["wilcoxon_stat"], res["p_value"] = float("nan"), float("nan")
    else:
        res["wilcoxon_stat"], res["p_value"] = float("nan"), 1.0
    return res


if __name__ == "__main__":
    # Self-test: a perfect ranking should max out; a random ranking should sit near baseline.
    rng = np.random.default_rng(42)
    N, n_act = 500, 25
    labels = np.array([1]*n_act + [0]*(N-n_act))
    perfect = np.array([1.0]*n_act + [0.0]*(N-n_act))       # actives ranked first
    random_ = rng.random(N)
    worst = np.array([0.0]*n_act + [1.0]*(N-n_act))         # actives ranked last

    print("PERFECT:", {k: round(v, 3) for k, v in score_ranking(perfect, labels).items()})
    print("RANDOM :", {k: round(v, 3) for k, v in score_ranking(random_, labels).items()})
    print("WORST  :", {k: round(v, 3) for k, v in score_ranking(worst, labels).items()})
    p, lo, hi = bootstrap_ci(random_, labels, "AUROC", n_boot=500)
    print(f"RANDOM AUROC bootstrap: {p:.3f} [{lo:.3f}, {hi:.3f}]  (expect ~0.5)")
    # paired panel: platform slightly better than baseline on 3 of 4 targets
    pm = {"T1": 5.0, "T2": 3.0, "T3": 2.0, "T4": 1.0}
    bl = {"T1": 2.0, "T2": 2.5, "T3": 1.0, "T4": 1.5}
    print("PANEL TEST:", {k: (round(v, 3) if isinstance(v, float) else v)
                          for k, v in paired_panel_test(pm, bl).items()})
