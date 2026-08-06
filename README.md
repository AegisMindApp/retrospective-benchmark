# Blind Retrospective Benchmark — pre-registered virtual screening study

> **This is a public extract** of a directory inside a private monorepo. Commit SHAs here
> are new and do **not** include the original pre-registration commit `92c6f8bf`
> (2026-07-14) — see `PROVENANCE.md`, which states plainly what can and cannot be verified
> from this repository.
>
> **Archived deposit:** [10.5281/zenodo.21824633](https://doi.org/10.5281/zenodo.21824633)
>
> **The headline finding is not the one this README was originally written for.** The
> question below ("does the platform beat commodity docking?") could not be answered as
> posed: the debias gate excluded all but one of 14 targets. What the study reports instead
> is that **benchmark admissibility is unstable under active-set subsampling** — gate
> verdicts flip in *both* directions with the number of actives sampled, a variable that is
> near-universally capped for compute reasons and rarely reported.

Does the AegisMind platform predict real drug-discovery outcomes **better than commodity
docking, on data it could not have seen**? This harness answers that with a pre-registered,
leakage-controlled, reproducible protocol. It is the one form of validation we can mint
ourselves — no gatekeeper, no wet lab.

**Read `PREREGISTRATION.md` first** — it is the frozen protocol (committed before any result).
Provenance and what is independently verifiable: `PROVENANCE.md`.

## Files

| File | Role | Status |
|---|---|---|
| `PREREGISTRATION.md` | Frozen protocol — the credibility anchor | ✅ committed |
| `metrics.py` | EF@1/5%, BEDROC, AUROC (RDKit-validated) + bootstrap CI + Wilcoxon panel test | ✅ **self-tested** |
| `data_split.py` | ChEMBL temporal split + objective panel freeze → `panel_frozen.json` | ✅ code done; ChEMBL REST **reachability verified** |
| `baselines.py` | chance, fp-similarity (RDKit), Vina wrapper reuse | ✅ chance+fp **verified offline**; Vina needs binary |
| `platform_rank.py` | defines platform = z(−Vina) + z(ensemble); reuses the platform's own ensemble | ✅ compiles; needs vina + API keys to run |
| `run_pilot.py` | orchestrator; enforces blind→unblind order | ✅ compiles |

## Run sequence (each step is a checkpoint in the pre-registration)

```bash
cd analysis/retrospective_benchmark

# 0. FREEZE — already done: PREREGISTRATION.md is committed. Do not edit it post-hoc.

# 1. Build the frozen panel (data-availability only; platform never consulted).
#    Pass candidate ChEMBL target IDs; the qualifying subset is auto-selected.
python3 data_split.py CHEMBL279 CHEMBL240 CHEMBL4282 CHEMBL1824 ...
#    → writes panel_frozen.json ; commit it (freezes the panel).

# 2. Provide receptors.json for the Vina/platform paths:
#    { "CHEMBL279": {"target_name":"VEGFR2","context":"...",
#                    "receptor_pdbqt":"/path/3RXX_receptor.pdbqt",
#                    "box_center":[x,y,z],"box_size":[24,24,24]}, ... }
#    (chance + fp_similarity run without receptors; Vina/platform are skipped without them.)

# 3. Run — scores every method BLIND, then unblinds to metrics + panel test.
python3 run_pilot.py            # → predictions_*.json, report.json
```

## Dependencies

- **Present here:** rdkit, sklearn, scipy, numpy, pandas, requests. ChEMBL REST reachable.
- **Needed to run the docking/platform paths (not in this sandbox):**
  - **Vina binary** at `/tmp/vina` (or edit `VINA_BIN` in
    `vina_tools.py`) + Meeko for ligand prep.
  - **Prepared receptors** per target (PDBQT + binding box) — reuse
    PDB downloads; receptor preparation is in `prepare_receptors.py`.
  - **Provider API keys** (OPENAI/ANTHROPIC/GOOGLE/MISTRAL/XAI) for the ensemble layer.
  - Compute: docking is CPU-cheap; the ensemble layer is the cost driver (bounded by
    design — ranks a docked shortlist, not every compound blindly).

## Honest status of this commit

Built and **verified runnable here:** the metrics (self-tested against perfect/random/worst
rankings), the temporal-split data logic (ChEMBL reachable, fields + year-lookup confirmed),
and the chance/fingerprint baselines. The **Vina and ensemble paths are wired to existing repo
infra and compile**, but need the vina binary + receptors + API keys to execute — i.e. a
machine with the docking stack (the TPU/worker box that ran the AMR/HD screens), not this
sandbox. Nothing has been scored against ground truth yet — by design, the panel isn't frozen
until step 1 runs there.

## Integrity guarantees baked into the code

- Labels are never read during scoring; `run_pilot.py` writes all rankings before unblinding.
- Panel selection uses data-availability only — the platform can't be steered onto easy targets.
- A method missing its inputs is **skipped and reported**, never silently substituted.
- Frozen constants live once (in the modules) and are mirrored in `PREREGISTRATION.md`.


## Obtaining AutoDock Vina

The Vina binary is **not** redistributed here. The study used **v1.2.7** (Linux x86-64):

```
https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7
md5  0c5d02550bdb3e661a77ee0badcfbe78
```

Place it at `bin/vina`, or point `$VINA_BIN` at it. The md5 is recorded because the
protocol identity depends on the exact binary: caches carry a sidecar recording the md5 of
the binary and receptor, and a run refuses to resume against a cache built under a
different one. An earlier version of this code asserted a hardcoded version string that
silently became false when the binary was upgraded — hence digests rather than strings.

## Reproducing

Paths are repo-relative. Set `$OBABEL` if Open Babel is not on `PATH`. The decoy generator
queries the live ChEMBL REST API, so decoy sets are **not** reproducible on demand — a
frozen, md5-stamped compound list is provided for the PD-L1 run for exactly that reason
(`kaggle/frozen_compounds_*.json`, if present).

## Distributed execution

`kaggle/` contains the harness for splitting a run across a local machine and Kaggle CPU
sessions. It ships pre-built ligand PDBQTs, the binary and the receptor rather than
reproducing any of them remotely, and verifies md5s on arrival. Five compounds scored in
both environments agreed to within 0.06 kcal/mol against a pre-committed 0.30 tolerance.
Set `$KAGGLE_USER` and provide `KAGGLE_API_TOKEN` via `$ENV_FILE`.
