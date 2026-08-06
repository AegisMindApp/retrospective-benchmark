#!/usr/bin/env python3
"""
prepare_ligands.py — SMILES → 3D → PDBQT for each compound in the panel.

Uses RDKit for conformer generation + obabel for PDBQT conversion.
Writes ligand PDBQT files to ligands/<target>/<molecule_id>.pdbqt
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

OBABEL = os.getenv("OBABEL", "obabel")
BENCH  = Path(__file__).parent
LIGS   = BENCH / "ligands"


def smiles_to_pdbqt(smiles: str, mol_id: str, out_path: Path) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result != 0:
        return False
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tmp:
        tmp_sdf = tmp.name
    writer = Chem.SDWriter(tmp_sdf)
    writer.write(mol)
    writer.close()

    r = subprocess.run(
        [OBABEL, tmp_sdf, "-O", str(out_path),
         "--partialcharge", "gasteiger", "-h"],
        capture_output=True, text=True, timeout=30
    )
    os.unlink(tmp_sdf)
    return out_path.exists() and out_path.stat().st_size > 0


def prepare_panel_ligands(panel_path="panel_frozen.json") -> dict:
    with open(panel_path) as f:
        panel = json.load(f)

    stats = {}
    for tid, compounds in panel["targets"].items():
        out_dir = LIGS / tid
        out_dir.mkdir(parents=True, exist_ok=True)
        ok = fail = skip = 0
        for c in compounds:
            mol_id = c["molecule_chembl_id"]
            out_pdbqt = out_dir / f"{mol_id}.pdbqt"
            if out_pdbqt.exists():
                ok += 1
                continue
            if smiles_to_pdbqt(c["smiles"], mol_id, out_pdbqt):
                ok += 1
            else:
                fail += 1
                print(f"  FAIL: {mol_id} ({c['smiles'][:40]})")
        stats[tid] = {"ok": ok, "fail": fail}
        print(f"{tid}: {ok} ok, {fail} failed")
    return stats


if __name__ == "__main__":
    stats = prepare_panel_ligands()
    print("\nLigand prep complete:", stats)
