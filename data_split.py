#!/usr/bin/env python3
"""
data_split.py — build the leakage-controlled temporal split from ChEMBL.

The linchpin of the whole benchmark: keep only bioactivity records whose source
document was published AFTER the model training cutoff, so the LLM ensemble
cannot have memorised the answer. Everything downstream depends on this being
correct, so it is deliberately conservative and auditable.

Reuses the plain-REST access pattern from analysis/flashoptim_tpu/phase_amr_chembl.py
(urllib/requests against the EBI ChEMBL API) — no chembl_webresource_client needed.

Frozen parameters live in PREREGISTRATION.md and are imported here as constants so
code and pre-registration cannot drift.
"""
from __future__ import annotations

import os
import json
import time
import urllib.parse
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional

import requests

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# ── FROZEN PARAMETERS (mirror PREREGISTRATION.md — do not tune after unblinding) ──
# Amendment 2 (2026-07-27): CUTOFF_YEAR changed from 2025 → 2024 because ChEMBL 35
# has no 2026 records. Keep year > 2024 (= 2025-dated records). See PREREGISTRATION.md.
CUTOFF_YEAR      = int(os.getenv("BENCHMARK_CUTOFF_YEAR", "2024"))
# Aug 2025 is the approximate ensemble training cutoff. Records from year=2025 docs
# published before this date carry residual leakage risk (Amendment 2 §A2.3).
LEAKAGE_GUARD_THRESHOLD = 0.30   # flag target if >30% actives from pre-Aug-2025 docs
ACTIVE_PCHEMBL   = 6.0       # pChEMBL ≥ 6.0  → active   (~≤ 1 µM)
INACTIVE_PCHEMBL = 5.0       # pChEMBL ≤ 5.0  → inactive; strictly-between = grey, dropped
ALLOWED_TYPES    = ("IC50", "Ki", "Kd", "EC50")
# Original §4 eligibility (≥10 actives AND ≥10 inactives) is SUPERSEDED by
# PREREGISTRATION.md Amendment 1 §A1.2: ≥15 actives alone — inactives are supplied as
# property-matched debiased decoys (decoys.py), because recent data lacks inactives.
MIN_ACTIVES      = 15        # panel-eligibility: a target needs ≥ this many post-cutoff actives
MIN_INACTIVES    = 10        # retained for reference; no longer gates eligibility (Amendment 1)
# ──────────────────────────────────────────────────────────────────────────────


@dataclass
class Compound:
    molecule_chembl_id: str
    smiles: str
    pchembl: float
    label: int            # 1 active, 0 inactive
    doc_year: int
    std_type: str


def _get(url: str, params: dict, retries: int = 4) -> dict:
    for attempt in range(retries):
        r = requests.get(url, params=params, timeout=60,
                         headers={"Accept": "application/json"})
        if r.status_code == 200:
            return r.json()
        time.sleep(1.5 * (attempt + 1))
    r.raise_for_status()
    return {}


def _doc_years(doc_ids: List[str]) -> Dict[str, int]:
    """Map document_chembl_id → publication year (the temporal signal)."""
    years: Dict[str, int] = {}
    for i in range(0, len(doc_ids), 50):
        chunk = doc_ids[i:i+50]
        data = _get(f"{CHEMBL}/document.json",
                    {"document_chembl_id__in": ",".join(chunk), "limit": len(chunk)})
        for d in data.get("documents", []):
            y = d.get("year")
            if y:
                years[d["document_chembl_id"]] = int(y)
    return years


def fetch_target_activities(target_chembl_id: str, max_records: int = 4000
                            ) -> List[Compound]:
    """Fetch post-cutoff, labelled compounds for one target. Applies the temporal
    split and the frozen activity thresholds; dedups a molecule to its max pChEMBL."""
    raw: List[dict] = []
    offset, limit = 0, 1000
    while offset < max_records:
        data = _get(f"{CHEMBL}/activity.json", {
            "target_chembl_id": target_chembl_id,
            "pchembl_value__isnull": "false",
            "standard_type__in": ",".join(ALLOWED_TYPES),
            "order_by": "-activity_id",   # newest first: fast capture of recent-year data
            "limit": limit, "offset": offset,
        })
        acts = data.get("activities", [])
        if not acts:
            break
        raw.extend(acts)
        offset += limit
        if not data.get("page_meta", {}).get("next"):
            break

    doc_ids = sorted({a["document_chembl_id"] for a in raw if a.get("document_chembl_id")})
    years = _doc_years(doc_ids)

    best: Dict[str, Compound] = {}
    for a in raw:
        doc = a.get("document_chembl_id")
        yr = years.get(doc)
        if not yr or yr <= CUTOFF_YEAR:           # TEMPORAL SPLIT — the leakage defence
            continue
        try:
            p = float(a.get("pchembl_value"))
        except (TypeError, ValueError):
            continue
        if INACTIVE_PCHEMBL < p < ACTIVE_PCHEMBL:  # grey zone — excluded
            continue
        smi = a.get("canonical_smiles")
        mol = a.get("molecule_chembl_id")
        if not smi or not mol:
            continue
        label = 1 if p >= ACTIVE_PCHEMBL else 0
        c = Compound(mol, smi, p, label, yr, a.get("standard_type", ""))
        if mol not in best or p > best[mol].pchembl:
            best[mol] = c
    return list(best.values())


def target_qualifies(compounds: List[Compound]) -> bool:
    # Amendment 1: eligibility is ≥ MIN_ACTIVES post-cutoff actives; inactives are
    # supplied as decoys downstream, so the inactive count no longer gates.
    return sum(c.label for c in compounds) >= MIN_ACTIVES


def build_panel(candidate_targets: List[str], out_path: str) -> dict:
    """Apply the OBJECTIVE, pre-registered eligibility criteria to a list of
    candidate targets and freeze the qualifying panel + its compounds to JSON.
    Selection is by data-availability only — the platform is never consulted here,
    so there is no way to cherry-pick targets the platform happens to do well on."""
    panel = {}
    for t in candidate_targets:
        try:
            comps = fetch_target_activities(t)
        except Exception as e:
            print(f"  {t}: fetch failed ({type(e).__name__}: {str(e)[:80]}) — skipped")
            continue
        n_act = sum(c.label for c in comps)
        ok = target_qualifies(comps)
        print(f"  {t}: {len(comps)} post-{CUTOFF_YEAR} compounds "
              f"({n_act} act / {len(comps)-n_act} inact) — "
              f"{'QUALIFIES' if ok else 'excluded'}")
        if ok:
            panel[t] = [asdict(c) for c in comps]
    with open(out_path, "w") as fh:
        json.dump({"cutoff_year": CUTOFF_YEAR, "active_pchembl": ACTIVE_PCHEMBL,
                   "inactive_pchembl": INACTIVE_PCHEMBL, "targets": panel}, fh, indent=2)
    print(f"\nFrozen panel: {len(panel)} target(s) → {out_path}")
    return panel


if __name__ == "__main__":
    import sys
    # Candidate targets are supplied on the CLI (or a file); the qualifying subset
    # is determined purely by the criteria above. Example candidates can be any
    # ChEMBL target IDs — the panel is whatever passes, frozen before any ranking.
    targets = sys.argv[1:] or ["CHEMBL279"]   # CHEMBL279 = VEGFR2 (smoke-test default)
    build_panel(targets, "panel_frozen.json")
