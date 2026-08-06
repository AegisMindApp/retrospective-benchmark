#!/usr/bin/env python3
"""
platform_rank.py — the definition of "the platform's prediction" under test.

The whole benchmark turns on one question: does the FULL platform beat the free
Vina it wraps? So the platform score must be docking PLUS the multi-model
reasoning layer — otherwise we'd just be re-testing Vina against itself.

Platform score for a (compound, target) pair =
    z(−vina_affinity)  +  z(ensemble_activity_likelihood)
i.e. the docking signal combined with an independent multi-model assessment of
how likely the compound is active given the target context. If the reasoning
layer adds nothing, this collapses toward Vina and the paired panel test will
(correctly) show no edge — that is the honest outcome the design is built to detect.

LEAKAGE NOTE (see PREREGISTRATION.md §2): the ensemble prompt gives structure +
target context but is NEVER told the measured label, and the compounds are
post-training-cutoff, so the ensemble cannot simply recall the answer. If the
ensemble's only edge were memorisation, the temporal split blunts it — which is
exactly what we want to measure honestly.
"""
from __future__ import annotations

import os
import sys
from typing import List, Sequence, Optional

import numpy as np

from baselines import vina_scores


_ACTIVITY_PROMPT = (
    "You are assessing whether a small molecule is likely to be a potent binder/inhibitor "
    "of a specific protein target, for a virtual-screening prioritisation. Judge ONLY from "
    "the chemistry and the target context below — do not rely on recalling any specific "
    "published assay result.\n\n"
    "Target: {target}\nContext: {context}\nCompound SMILES: {smiles}\n\n"
    "Reply with ONLY a number from 0.0 to 1.0 = your estimated probability the compound is a "
    "potent (sub-µM) binder of this target, then a 5-word reason."
)


def _zscore(x: Sequence[float]) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    m = np.nanmean(a)
    s = np.nanstd(a)
    a = np.where(np.isnan(a), m, a)          # impute failed dockings to the mean (rank-neutral)
    return (a - m) / s if s > 0 else np.zeros_like(a)


async def ensemble_activity_scores(smiles_list: Sequence[str], target: str, context: str
                                   ) -> List[float]:
    """Independent multi-model activity-likelihood per compound (0..1), averaged
    across the ensemble. Reuses the platform's own ensemble so this is genuinely
    'the platform', not a bespoke model. Needs provider API keys."""
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
    from aegismind import AegisMind                       # noqa
    from autonomous.model_utils import get_available_models
    aegis = AegisMind(models=get_available_models(include_mistral=True, quiet=True),
                      max_iterations=1, enable_memory=False, quiet_init=True)
    import re
    scores = []
    for smi in smiles_list:
        prompt = _ACTIVITY_PROMPT.format(target=target, context=context, smiles=smi)
        resp = await aegis.ensemble.generate_parallel(prompt, max_tokens=1000, temperature=0.1)
        vals = []
        for text in (resp or {}).values():
            t = text.get("content", "") if isinstance(text, dict) else str(text)
            m = re.search(r"(0?\.\d+|[01]\.?\d*)", t)
            if m:
                try:
                    v = float(m.group(1))
                    if 0.0 <= v <= 1.0:
                        vals.append(v)
                except ValueError:
                    pass
        scores.append(float(np.mean(vals)) if vals else 0.5)
    return scores


async def platform_scores(smiles_list: Sequence[str], target: str, context: str,
                          receptor_pdbqt: str, box_center, box_size,
                          names: Optional[Sequence[str]] = None,
                          use_ensemble: bool = True) -> List[float]:
    """The frozen platform-score definition: z(docking) + z(ensemble)."""
    vina = vina_scores(smiles_list, receptor_pdbqt, box_center, box_size, names)
    z_dock = _zscore(vina)
    if not use_ensemble:
        return list(z_dock)                              # docking-only ablation
    ens = await ensemble_activity_scores(smiles_list, target, context)
    return list(z_dock + _zscore(ens))


if __name__ == "__main__":
    print("platform_rank defines: platform = z(-vina) + z(ensemble_activity).")
    print("Run via run_pilot.py once panel_frozen.json + the vina binary are in place.")
