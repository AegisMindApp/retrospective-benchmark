#!/usr/bin/env python3
"""
prepare_receptors.py — download PDB structures, strip non-protein, convert to PDBQT.

Writes receptors/<CHEMBL_ID>_receptor.pdbqt and populates receptors.json with
box centre/size derived from co-crystallised ligand centroid.
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import numpy as np

BENCH  = Path(__file__).resolve().parent
RECEP  = BENCH / "receptors"
RECEP.mkdir(exist_ok=True)

# Use the manual PDBQT writer from the AMR docking module — it produces a
# properly rigid receptor (no ROOT/BRANCH/TORSDOF ligand-tree markers) that
# AutoDock Vina 1.2 accepts via --receptor.  obabel's PDB→PDBQT treats the
# protein as a flexible molecule and inserts ROOT/BRANCH records, which Vina
# rejects: "PDBQT parsing error: Unknown or inappropriate tag found in rigid
# receptor."
from vina_tools import pdb_to_pdbqt_receptor

# PDB codes and known binding-box parameters for each target.
# Box centre is derived from co-crystallised ligand centroid; size 22 Å covers the ATP pocket.
RECEPTOR_MAP = {
    # ── kinase panel (original) ──────────────────────────────────────────────────
    "CHEMBL279":    {"pdb": "4ASD", "box_size": [22, 22, 22],
                     "target_name": "VEGFR2 (KDR)",
                     "context": "vascular endothelial growth factor receptor 2, kinase domain, ATP site"},
    "CHEMBL203":    {"pdb": "4ZAU", "box_size": [22, 22, 22],
                     "target_name": "EGFR",
                     "context": "epidermal growth factor receptor, kinase domain, erlotinib binding site"},
    "CHEMBL301":    {"pdb": "3QQK", "box_size": [20, 20, 20],
                     "target_name": "CDK2",
                     "context": "cyclin-dependent kinase 2, ATP binding site, cell cycle regulator"},
    "CHEMBL2971":   {"pdb": "3LPB", "box_size": [22, 22, 22],
                     "target_name": "JAK2",
                     "context": "Janus kinase 2, ATP binding site, involved in cytokine signalling"},
    "CHEMBL205":    {"pdb": "1UWH", "box_size": [22, 22, 22],
                     "target_name": "BRAF",
                     "context": "BRAF serine/threonine kinase, ATP binding site, MAPK pathway"},
    "CHEMBL325":    {"pdb": "2BCJ", "box_size": [20, 20, 20],
                     "target_name": "Aurora A kinase",
                     "context": "Aurora kinase A, ATP binding site, mitotic spindle assembly"},
    "CHEMBL2147":   {"pdb": "2XP2", "box_size": [22, 22, 22],
                     "target_name": "ALK",
                     "context": "anaplastic lymphoma kinase, ATP binding site"},
    "CHEMBL1075104":{"pdb": "3GEN", "box_size": [22, 22, 22],
                     "target_name": "BTK",
                     "context": "Bruton's tyrosine kinase, ATP binding site, B-cell signalling"},
    # ── diverse non-kinase panel (Amendment 4) ───────────────────────────────────
    "CHEMBL4523582":{"pdb": "7K40", "box_size": [22, 22, 22],
                     "target_name": "SARS-CoV-2 main protease (Mpro)",
                     "context": "SARS-CoV-2 Mpro cysteine protease, substrate binding site S1/S2 pockets"},
    "CHEMBL244":    {"pdb": "2J34", "box_size": [22, 22, 22],
                     "target_name": "Factor Xa",
                     "context": "coagulation factor Xa serine protease, S1 binding pocket, anticoagulant target"},
    "CHEMBL1163125":{"pdb": "3MXF", "box_size": [20, 20, 20],
                     "target_name": "BRD4 bromodomain 1",
                     "context": "BRD4 bromodomain 1, acetyl-lysine binding site, epigenetic reader"},
    # ── high-diversity panel (Amendment 5) ───────────────────────────────────────
    "CHEMBL612545": {"pdb": "5J89", "box_size": [22, 22, 22],
                     "target_name": "PD-L1",
                     "context": "PD-L1 immune checkpoint protein, C-D loop beta-sheet PPI inhibitor binding site"},
    "CHEMBL3717":   {"pdb": "4EY6", "box_size": [22, 22, 22],
                     "target_name": "Acetylcholinesterase",
                     "context": "acetylcholinesterase active site gorge, catalytic triad Ser203/His447/Glu334"},
    "CHEMBL4005":   {"pdb": "2IKH", "box_size": [20, 20, 20],
                     "target_name": "Aldose reductase",
                     "context": "aldose reductase NADPH-binding site, anion binding subsite, diabetic complication target"},
    "CHEMBL3778":   {"pdb": "4JRG", "box_size": [20, 20, 20],
                     "target_name": "MDM2",
                     "context": "MDM2 p53-binding domain, Phe19/Trp23/Leu26 hydrophobic hotspot, PPI inhibitor site"},
}


def download_pdb(pdb_id: str, out_path: Path) -> bool:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, str(out_path))
        return out_path.exists()
    except Exception as e:
        print(f"  download failed: {e}")
        return False


def extract_ligand_center(pdb_path: Path) -> list[float] | None:
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() not in ("HOH", "WAT", "EDO", "GOL", "PEG"):
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    pass
    if not coords:
        return None
    arr = np.array(coords)
    return [round(float(v), 2) for v in arr.mean(axis=0)]


def strip_protein(pdb_path: Path, out_path: Path):
    with open(pdb_path) as f, open(out_path, "w") as out:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM":
                out.write(line)
            elif rec in ("TER", "END"):
                out.write(line)


def to_pdbqt(pdb_path: Path, out_path: Path) -> bool:
    """Convert protein PDB to rigid receptor PDBQT (no tree markup)."""
    try:
        pdb_to_pdbqt_receptor(str(pdb_path), str(out_path))
        return out_path.exists() and out_path.stat().st_size > 0
    except Exception as e:
        print(f"  pdb_to_pdbqt_receptor failed: {e}")
        return False


def prepare_receptors(targets: list[str]) -> dict:
    result = {}
    for tid in targets:
        if tid not in RECEPTOR_MAP:
            print(f"{tid}: no PDB mapping — skipped")
            continue
        cfg = RECEPTOR_MAP[tid]
        pdb_id = cfg["pdb"]
        pdb_raw  = RECEP / f"{pdb_id}_raw.pdb"
        pdb_prot = RECEP / f"{pdb_id}_protein.pdb"
        pdbqt    = RECEP / f"{tid}_receptor.pdbqt"

        if pdbqt.exists():
            print(f"{tid}: already prepared — {pdbqt.stat().st_size // 1024}KB")
        else:
            print(f"{tid}: downloading {pdb_id}...")
            if not pdb_raw.exists():
                if not download_pdb(pdb_id, pdb_raw):
                    print(f"  {tid}: FAILED download")
                    continue
            strip_protein(pdb_raw, pdb_prot)
            if not to_pdbqt(pdb_prot, pdbqt):
                print(f"  {tid}: FAILED pdbqt conversion")
                continue
            print(f"  {tid}: prepared {pdbqt.stat().st_size // 1024}KB")

        center = cfg.get("box_center") or extract_ligand_center(pdb_raw) or [0.0, 0.0, 0.0]
        result[tid] = {
            "receptor_pdbqt": str(pdbqt),
            "box_center": center,
            "box_size": cfg["box_size"],
            "pdb_id": pdb_id,
            "target_name": cfg.get("target_name", tid),
            "context": cfg.get("context", ""),
        }

    out = BENCH / "receptors.json"
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nreceptors.json written: {len(result)} targets")
    return result


if __name__ == "__main__":
    import sys
    targets = sys.argv[1:] or list(RECEPTOR_MAP.keys())
    prepare_receptors(targets)
