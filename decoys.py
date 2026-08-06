#!/usr/bin/env python3
"""
decoys.py — property-matched, debiased decoys (PRE-REGISTRATION Amendment 1).

Recent bioactivity data is active-dominated, so an enrichment benchmark needs
synthetic negatives. These decoys are: (1) property-matched to the actives,
(2) topologically dissimilar to every active (no latent-active / analogue
leakage), (3) free of any recorded activity against the target. A hostile-reviewer
DEBIAS GATE then verifies a physicochemical-only classifier CANNOT separate actives
from decoys — otherwise "enrichment" would be a property artifact, not recognition.

All frozen parameters mirror PREREGISTRATION.md Amendment 1 §A1.4/§A1.5.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

from data_split import _get, CHEMBL

# ── FROZEN (Amendment 1 §A1.4/§A1.5; TPSA window added Amendment 7; server PSA filter removed Amendment 9) ──
DECOY_RATIO      = 25        # decoys per active
MAX_ACTIVES      = 100       # Amendment 16 (2026-08-04): was 30; the cap discarded 55 of
                             # 85 Mpro actives and weakened the primary endpoint. Seed 0.
MIN_ACTIVES      = 15        # panel eligibility (supersedes the original ≥10/≥10)
MW_WIN, LOGP_WIN = 25.0, 0.5
HBD_WIN, HBA_WIN, RTB_WIN = 1, 1, 1     # ± windows; formal charge must match EXACTLY
TPSA_WIN         = 20.0      # Amendment 7: ±20 Å² TPSA window (matches DUDE-Z standard)
MAX_SIM          = 0.35      # decoy↔active max Tanimoto (Morgan r2, 2048)
DEBIAS_AUROC_MAX = 0.60      # accept decoys only if property-only classifier AUROC ≤ this
SEED             = 0
FP_R, FP_BITS    = 2, 2048
_DESC_KEYS       = ("mw", "clogp", "hbd", "hba", "rtb", "tpsa", "charge")


def descriptors(smiles: str) -> Optional[Dict[str, float]]:
    m = Chem.MolFromSmiles(smiles or "")
    if m is None:
        return None
    return {"mw": Descriptors.MolWt(m), "clogp": Crippen.MolLogP(m),
            "hbd": Lipinski.NumHDonors(m), "hba": Lipinski.NumHAcceptors(m),
            "rtb": Descriptors.NumRotatableBonds(m), "tpsa": Descriptors.TPSA(m),
            "charge": Chem.GetFormalCharge(m)}


def _fp(smiles: str):
    m = Chem.MolFromSmiles(smiles or "")
    return AllChem.GetMorganFingerprintAsBitVect(m, FP_R, FP_BITS) if m else None


def scaffold_count(smiles_list: Sequence[str]) -> int:
    """Distinct Bemis–Murcko scaffolds among the actives (analogue-bias transparency)."""
    scafs = set()
    for smi in smiles_list:
        try:
            s = MurckoScaffold.MurckoScaffoldSmiles(smi)
            if s:
                scafs.add(s)
        except Exception:
            pass
    return len(scafs)


def property_match(a: Dict[str, float], c: Dict[str, float]) -> bool:
    return (abs(a["mw"] - c["mw"]) <= MW_WIN and abs(a["clogp"] - c["clogp"]) <= LOGP_WIN
            and abs(a["hbd"] - c["hbd"]) <= HBD_WIN and abs(a["hba"] - c["hba"]) <= HBA_WIN
            and abs(a["rtb"] - c["rtb"]) <= RTB_WIN and a["charge"] == c["charge"]
            and abs(a["tpsa"] - c["tpsa"]) <= TPSA_WIN)  # Amendment 7


def _fetch_candidates(a: Dict[str, float], limit: int = 500) -> List[Tuple[str, str]]:
    """Drug-like ChEMBL molecules inside the active's property windows (server-side)."""
    # Amendment 9: PSA server filter removed — ChEMBL molecule_properties__psa uses a
    # different computation than RDKit Descriptors.TPSA, causing inconsistent filtering
    # that worsened gate AUROC across all targets in run9 vs run8.  Client-side RDKit
    # TPSA check in property_match() provides consistent matching.
    # Amendment 11: add order_by=molecule_chembl_id for deterministic API responses.
    # Without ordering, ChEMBL returns an arbitrary 500-compound subset on each call,
    # making the decoy set non-reproducible even with a fixed rng seed.
    try:
        data = _get(f"{CHEMBL}/molecule.json", {
            "molecule_properties__full_mwt__gte": round(a["mw"] - MW_WIN, 1),
            "molecule_properties__full_mwt__lte": round(a["mw"] + MW_WIN, 1),
            "molecule_properties__alogp__gte": round(a["clogp"] - LOGP_WIN, 2),
            "molecule_properties__alogp__lte": round(a["clogp"] + LOGP_WIN, 2),
            "molecule_properties__hbd__gte": int(a["hbd"] - HBD_WIN),
            "molecule_properties__hbd__lte": int(a["hbd"] + HBD_WIN),
            "molecule_properties__hba__gte": int(a["hba"] - HBA_WIN),
            "molecule_properties__hba__lte": int(a["hba"] + HBA_WIN),
            "molecule_structures__canonical_smiles__isnull": "false",
            "order_by": "molecule_chembl_id",   # Amendment 11: deterministic ordering
            "limit": limit,
        })
    except Exception:
        return []
    out = []
    for m in data.get("molecules", []):
        smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
        cid = m.get("molecule_chembl_id")
        if smi and cid:
            out.append((cid, smi))
    return out


def debias_auroc(active_props: List[Dict], decoy_props: List[Dict]) -> float:
    """5-fold-CV AUROC of a logistic regression on the 7 physicochemical descriptors,
    separating actives from decoys. LOW (≤0.60) = decoys are NOT property-separable."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.model_selection import StratifiedKFold, cross_val_score
    X = np.array([[p[k] for k in _DESC_KEYS] for p in (active_props + decoy_props)], float)
    y = np.array([1] * len(active_props) + [0] * len(decoy_props))
    n_splits = int(min(5, np.bincount(y).min()))
    if n_splits < 2:
        return float("nan")
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    return float(np.mean(cross_val_score(clf, X, y, cv=cv, scoring="roc_auc")))


def _select(actives, active_fps, exclude_ids, attempt: int) -> List[Dict]:
    """Draw a decoy set (deterministic per attempt via SEED+attempt)."""
    rng = np.random.default_rng(SEED + attempt)
    need = DECOY_RATIO * len(actives)
    chosen: Dict[str, Dict] = {}
    order = list(range(len(actives)))
    rng.shuffle(order)
    for i in order:
        a = actives[i]
        cands = _fetch_candidates(a["props"])
        rng.shuffle(cands)
        for cid, smi in cands:
            if cid in exclude_ids or cid in chosen:
                continue
            cp = descriptors(smi)
            if cp is None or not property_match(a["props"], cp):
                continue
            fp = _fp(smi)
            if fp is None:
                continue
            if max(DataStructs.BulkTanimotoSimilarity(fp, active_fps)) >= MAX_SIM:
                continue                      # too similar to some active → drop
            chosen[cid] = {"molecule_chembl_id": cid, "smiles": smi, "props": cp, "label": 0}
            if len(chosen) >= need:
                return list(chosen.values())
    return list(chosen.values())


def build_decoys(actives: List[Dict], target_chembl_id: str) -> Dict:
    """actives: [{molecule_chembl_id, smiles, ...}]. Returns decoys + debias/scaffold
    report + accept/exclude status per Amendment 1. Subsamples to MAX_ACTIVES (seed 0)."""
    rng = np.random.default_rng(SEED)
    if len(actives) > MAX_ACTIVES:
        idx = rng.choice(len(actives), MAX_ACTIVES, replace=False)
        actives = [actives[i] for i in sorted(idx)]
    for a in actives:
        a.setdefault("props", descriptors(a["smiles"]))
        a["label"] = 1
    # Drop macromolecules (MW>800 Da) — virtual screening targets small molecules only;
    # macromolecule actives would cause ChEMBL API timeouts and no property-matched decoys exist.
    actives = [a for a in actives if a["props"] is not None and a["props"]["mw"] <= 800]
    active_fps = [_fp(a["smiles"]) for a in actives if _fp(a["smiles"]) is not None]
    exclude_ids = {a["molecule_chembl_id"] for a in actives}
    scaf = scaffold_count([a["smiles"] for a in actives])

    for attempt in range(3):                  # ≤3 re-draws (Amendment 1 §A1.5)
        decoys = _select(actives, active_fps, exclude_ids, attempt)
        if len(decoys) < DECOY_RATIO * len(actives) * 0.5:
            continue                          # couldn't source enough — try another draw
        auroc = debias_auroc([a["props"] for a in actives], [d["props"] for d in decoys])
        if not (auroc == auroc) or auroc <= DEBIAS_AUROC_MAX:   # NaN or passes gate
            return {"target": target_chembl_id, "status": "ok", "n_actives": len(actives),
                    "n_decoys": len(decoys), "debias_auroc": auroc, "scaffold_count": scaf,
                    "attempt": attempt, "actives": actives, "decoys": decoys}
    return {"target": target_chembl_id, "status": "excluded_property_separable",
            "n_actives": len(actives), "debias_auroc": auroc, "scaffold_count": scaf,
            "actives": actives, "decoys": []}


if __name__ == "__main__":
    # Offline self-test of the local logic (no network).
    drugs = ["CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
             "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "CCN(CC)CCNC(=O)c1cc(Cl)c(N)cc1OC"]
    print("descriptors(aspirin):", {k: round(v, 2) for k, v in descriptors(drugs[0]).items()})
    print("scaffold_count(4 drugs):", scaffold_count(drugs))
    a0, a1 = descriptors(drugs[0]), descriptors(drugs[1])
    print("property_match(aspirin,caffeine):", property_match(a0, a1))
    # debias gate: separable groups → high AUROC (REJECT); overlapping → ~0.5 (ACCEPT)
    rng = np.random.default_rng(0)
    def synth(mu):
        return [{"mw": mu+rng.normal(0,5), "clogp": mu/100+rng.normal(0,.2), "hbd": 2, "hba": 4,
                 "rtb": 5, "tpsa": 60, "charge": 0} for _ in range(40)]
    print("debias AUROC (separable mw 200 vs 500):", round(debias_auroc(synth(200), synth(500)), 3), "→ expect HIGH (reject)")
    print("debias AUROC (overlapping mw 300 vs 300):", round(debias_auroc(synth(300), synth(300)), 3), "→ expect ~0.5 (accept)")
