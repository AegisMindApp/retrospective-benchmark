"""Generic AutoDock Vina helpers used by the benchmark.

Extracted from a larger internal docking script so this repository is self-contained.
Only protocol-neutral machinery is included: receptor and ligand preparation, the Vina
invocation, and version reporting. Target-specific constants, candidate compound
structures and patent-evidence generation are deliberately NOT part of this repository.

Ligand preparation is a hard gate. If Meeko cannot build a torsion tree the ligand is
EXCLUDED and counted, never silently written rigid — an earlier version fell back to a
writer emitting TORSDOF 0 with no BRANCH records, so every ligand was docked without
conformational search, invalidating an entire screen before it was noticed.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

DOCKING_DIR = Path(__file__).parent
RESULTS_DIR = DOCKING_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# Overridden per target by the caller; no meaningful default.
BOX_CENTER = (0.0, 0.0, 0.0)
BOX_SIZE = (22.0, 22.0, 22.0)

VINA_BIN = os.environ.get("VINA_BIN", str(DOCKING_DIR / "bin" / "vina"))

def pdb_to_pdbqt_receptor(pdb_path: str, pdbqt_path: str) -> None:
    """
    Minimal receptor PDBQT preparation:
    - Keep ATOM records only (strip HETATM, waters, ligands)
    - Assign AD4 atom types from element column
    - Set partial charge to 0.000 (Vina ignores receptor charges)
    """
    lines_out = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            # Skip water oxygens
            resname = line[17:20].strip()
            if resname in ("HOH", "WAT", "DOD"):
                continue
            atom_name = line[12:16].strip()
            element   = line[76:78].strip().upper() if len(line) >= 78 else atom_name[:1]
            if not element:
                element = atom_name[:1]
            ad4 = AD4_TYPES.get(element, element[:1])
            # PDBQT: cols 1-66 same as PDB, then charge (8.3f), atom type (right-aligned)
            pdb_cols = line[:66].ljust(66)
            # PDBQT: cols 67-70 spaces, 71-76 charge (%6.3f), 77 space, 78-79 atom type
            pdbqt_line = f"{pdb_cols}    {0.000:6.3f} {ad4:<2}\n"
            lines_out.append(pdbqt_line)
    lines_out.append("END\n")
    with open(pdbqt_path, "w") as fh:
        fh.writelines(lines_out)
    print(f"[receptor] Written {len(lines_out)-1} ATOM lines → {pdbqt_path}")


def smiles_to_pdbqt(smiles: str, name: str, out_path: str) -> bool:
    """
    SMILES → 3D conformer (RDKit ETKDGv3 + MMFF94) → Meeko PDBQT
    Returns True on success.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[ligand] ERROR: could not parse SMILES for {name}")
        return False

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result != 0:
        print(f"[ligand] WARNING: ETKDGv3 failed for {name}, trying random coords")
        AllChem.EmbedMolecule(mol, AllChem.ETDG())

    ff_result = AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    if ff_result not in (0, 1):
        print(f"[ligand] WARNING: MMFF did not converge for {name}")

    # Try Meeko
    try:
        from meeko import MoleculePreparation
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        if not mol_setups:
            raise RuntimeError("Meeko returned no setups")
        setup = mol_setups[0]
        # Meeko >=0.6 moved PDBQT serialisation out of the setup object into
        # PDBQTWriterLegacy; <=0.5 exposed setup.write_pdbqt_string(). Support both,
        # newest first, so this does not silently fail into the rigid-ligand path again.
        pdbqt_str = None
        try:
            from meeko import PDBQTWriterLegacy
            res = PDBQTWriterLegacy.write_string(setup)
            if isinstance(res, tuple):          # (string, ok, error) in 0.6/0.7
                pdbqt_str, ok, err = (list(res) + [None, None])[:3]
                if not pdbqt_str or ok is False:
                    raise RuntimeError(f"PDBQTWriterLegacy: {err}")
            else:
                pdbqt_str = res
        except ImportError:
            pdbqt_str = setup.write_pdbqt_string()
        if not pdbqt_str:
            raise RuntimeError("Meeko produced an empty PDBQT string")
        if "TORSDOF 0" in pdbqt_str and "BRANCH" not in pdbqt_str:
            raise RuntimeError("Meeko produced a rigid ligand (TORSDOF 0, no BRANCH)")
        with open(out_path, "w") as fh:
            fh.write(pdbqt_str)
        print(f"[ligand] {name}: Meeko PDBQT written ({len(pdbqt_str)} chars)")
        return True
    except Exception as e:
        print(f"[ligand] Meeko failed ({e})")

    # HARD GATE (Amendment 19, 2026-08-04). The fallback writer emits every atom in
    # a single ROOT block with TORSDOF 0 — a rigid ligand with no torsional search.
    # Under Python 3.7 Meeko failed on 100% of ligands, so an entire benchmark was
    # silently run on frozen molecules and read as evidence that docking does not
    # enrich. A ligand we cannot prepare properly is excluded and counted, never
    # docked rigid. Set BENCH_ALLOW_RIGID_LIGAND=1 only to reproduce the defect.
    if os.environ.get("BENCH_ALLOW_RIGID_LIGAND") == "1":
        print(f"[ligand] {name}: RIGID fallback (TORSDOF 0) — explicitly enabled")
        return _manual_pdbqt(mol, name, out_path)
    print(f"[ligand] {name}: EXCLUDED — no valid torsion tree. "
          f"Run under an environment where Meeko works (see Amendment 19).")
    return False


def _manual_pdbqt(mol, name: str, out_path: str) -> bool:
    """Fallback PDBQT writer using RDKit Gasteiger charges."""
    from rdkit.Chem import rdPartialCharges
    rdPartialCharges.ComputeGasteigerCharges(mol)

    conf = mol.GetConformer()
    lines = [f"REMARK  Name = {name}\n", "ROOT\n"]
    for i, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            continue  # skip hydrogens for brevity (Vina adds implicitly)
        pos   = conf.GetAtomPosition(i)
        chg   = float(atom.GetPropsAsDict().get("_GasteigerCharge", 0.0))
        el    = atom.GetSymbol().upper()
        ad4   = AD4_TYPES.get(el, el[:1])
        # PDBQT ATOM: cols 1-66 standard PDB, 67-70 spaces, 71-76 charge, 77 sp, 78-79 type
        lines.append(
            f"ATOM  {i+1:5d}  {el:<4}LIG A   1    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
            f"  1.00  0.00"
            f"    {chg:6.3f} {ad4:<2}\n"
        )
    lines += ["ENDROOT\n", "TORSDOF 0\n"]
    with open(out_path, "w") as fh:
        fh.writelines(lines)
    print(f"[ligand] {name}: manual PDBQT written")
    return True


def vina_version() -> str:
    """Ask the binary what it is, rather than asserting it.

    This was hardcoded to "1.2.5" and went stale when the binary was upgraded to
    1.2.7 on 27 Jul 2026 — a provenance record that silently stops being true is
    worse than none, because it is quoted into papers and patents. Any published
    artifact is only as trustworthy as the weakest claim in it.
    """
    try:
        out = subprocess.run([VINA_BIN, "--version"], capture_output=True,
                             text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError) as e:
        # narrow: a missing/unrunnable binary is the only thing worth tolerating.
        # A broad `except` here previously swallowed a NameError and reported
        # "unknown", which is the same silent-fallback shape as the rigid-ligand
        # bug — it looked like a provenance answer and was a defect.
        print(f"[vina] WARNING: cannot determine version ({e}); recording 'unknown'")
        return "unknown"
    m = re.search(r"v?(\d+\.\d+\.\d+)", out)
    if not m:
        print(f"[vina] WARNING: unparseable --version output: {out.strip()[:60]!r}")
        return "unknown"
    return m.group(1)


def run_vina(receptor_pdbqt: str, ligand_pdbqt: str, name: str,
             exhaustiveness: int = 16, n_poses: int = 5) -> dict:
    """Run AutoDock Vina and return parsed results."""
    out_pdbqt = str(RESULTS_DIR / f"{name}_out.pdbqt")
    cx, cy, cz = BOX_CENTER
    sx, sy, sz = BOX_SIZE

    cmd = [
        VINA_BIN,
        "--receptor", receptor_pdbqt,
        "--ligand",   ligand_pdbqt,
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x",   str(sx), "--size_y",   str(sy), "--size_z",   str(sz),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(n_poses),
        "--out", out_pdbqt,
    ]

    print(f"\n[vina] Docking {name} ...")
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300
        )
        stdout = proc.stdout + proc.stderr
    except subprocess.TimeoutExpired:
        return {"error": "timeout", "name": name}
    except Exception as e:
        return {"error": str(e), "name": name}

    # Parse affinity table from Vina output
    # Lines like:   1       -7.8      0.000      0.000
    affinities = []
    in_table = False
    for line in stdout.splitlines():
        if "-----" in line:
            in_table = True
            continue
        if in_table:
            parts = line.strip().split()
            if parts and parts[0].isdigit():
                try:
                    affinities.append(float(parts[1]))
                except (IndexError, ValueError):
                    pass

    if not affinities:
        print(f"[vina] WARNING: no affinities parsed for {name}")
        print("Vina output:", stdout[:500])
        return {"name": name, "error": "no_affinities", "raw": stdout[:500]}

    best = affinities[0]
    print(f"[vina] {name}: best pose {best:.2f} kcal/mol "
          f"(poses: {[round(a,2) for a in affinities]})")

    return {
        "name":            name,
        "best_affinity_kcal_mol": best,
        "all_poses_kcal_mol":     affinities,
        "n_poses":         len(affinities),
        "output_pdbqt":    out_pdbqt,
    }
