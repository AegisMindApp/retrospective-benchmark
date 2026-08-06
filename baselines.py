#!/usr/bin/env python3
"""
baselines.py — the reference methods every target is scored against, blind.

Three baselines, in ascending order of what they'd force the platform to beat:
  1. chance          — random scores (fixed seed) → the floor / null distribution.
  2. fp_similarity   — max Tanimoto to PRE-cutoff known actives (RDKit Morgan).
                       A cheap ML baseline; FAIR because it uses only pre-cutoff
                       knowledge, never the post-cutoff test labels.
  3. vina_docking    — −(best Vina affinity). The one that matters: does the full
                       platform beat the free docking tool it wraps? Reuses the
                       existing wrapper in vina_tools.py
                       (needs the vina binary at /tmp/vina + a prepared receptor).

Each baseline returns scores aligned to the input compound order, higher = more
likely active, so metrics.score_ranking() consumes them directly.
"""
from __future__ import annotations

import os
import sys
from typing import List, Sequence, Optional

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from data_split import _get, CHEMBL, CUTOFF_YEAR, ACTIVE_PCHEMBL, _doc_years


# ── 1. chance ─────────────────────────────────────────────────────────────────
def chance_scores(n: int, seed: int = 0) -> List[float]:
    return list(np.random.default_rng(seed).random(n))


# ── 2. fingerprint similarity to pre-cutoff actives ───────────────────────────
def _morgan(smiles: str, radius: int = 2, nbits: int = 2048):
    m = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(m, radius, nbits) if m else None


def _pre_cutoff_active_fps(target_chembl_id: str, max_records: int = 4000):
    """Morgan FPs of the target's KNOWN actives published on/before the cutoff —
    the only knowledge this baseline is allowed to use."""
    fps = []
    offset, limit = 0, 1000
    raw = []
    while offset < max_records:
        data = _get(f"{CHEMBL}/activity.json", {
            "target_chembl_id": target_chembl_id, "pchembl_value__isnull": "false",
            "standard_type__in": "IC50,Ki,Kd,EC50", "order_by": "-activity_id",
            "limit": limit, "offset": offset})
        acts = data.get("activities", [])
        if not acts:
            break
        raw.extend(acts); offset += limit
        if not data.get("page_meta", {}).get("next"):
            break
    years = _doc_years(sorted({a["document_chembl_id"] for a in raw
                               if a.get("document_chembl_id")}))
    for a in raw:
        yr = years.get(a.get("document_chembl_id"))
        try:
            p = float(a.get("pchembl_value"))
        except (TypeError, ValueError):
            continue
        if yr and yr <= CUTOFF_YEAR and p >= ACTIVE_PCHEMBL and a.get("canonical_smiles"):
            fp = _morgan(a["canonical_smiles"])
            if fp is not None:
                fps.append(fp)
    return fps


def fp_similarity_scores(smiles_list: Sequence[str], target_chembl_id: str) -> List[float]:
    ref = _pre_cutoff_active_fps(target_chembl_id)
    if not ref:
        return [0.0] * len(smiles_list)   # no prior actives → baseline is uninformative
    out = []
    for smi in smiles_list:
        fp = _morgan(smi)
        out.append(max(DataStructs.BulkTanimotoSimilarity(fp, ref)) if fp is not None else 0.0)
    return out


# ── 3. Vina docking (reuse existing wrapper) ──────────────────────────────────
def vina_scores(smiles_list: Sequence[str], receptor_pdbqt: str,
                box_center, box_size, names: Optional[Sequence[str]] = None
                ) -> List[float]:
    """−(best Vina affinity) per compound. Requires the vina binary (default
    $VINA_BIN) and a prepared receptor. Reuses vina_tools.smiles_to_pdbqt /
    run_vina so the docking is identical to the platform's own docking path."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import vina_tools as rvd  # type: ignore
    from vina_tools import smiles_to_pdbqt, run_vina, VINA_BIN, RESULTS_DIR  # type: ignore
    if not os.path.exists(VINA_BIN):
        raise FileNotFoundError(
            f"Vina binary not found at {VINA_BIN}. Install it (or set VINA_BIN) before "
            f"running the docking baseline — see README.md.")
    # run_vina reads the search box from module globals — set them per target here.
    rvd.BOX_CENTER = tuple(box_center)
    rvd.BOX_SIZE = tuple(box_size)
    scores = []
    names = names or [f"cmpd_{i}" for i in range(len(smiles_list))]
    for smi, name in zip(smiles_list, names):
        lig = str(RESULTS_DIR / f"{name}.pdbqt")
        if not smiles_to_pdbqt(smi, name, lig):
            scores.append(float("nan")); continue
        res = run_vina(receptor_pdbqt, lig, name)
        aff = res.get("best_affinity_kcal_mol")
        scores.append(-float(aff) if aff is not None else float("nan"))
    return scores


if __name__ == "__main__":
    # verifiable offline: chance + fp-similarity shape check on toy SMILES
    smis = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
    print("chance:", [round(x, 3) for x in chance_scores(len(smis))])
    fp = _morgan("CC(=O)Oc1ccccc1C(=O)O")
    print("morgan fp built:", fp is not None, "bits set:", int(fp.GetNumOnBits()) if fp else None)
