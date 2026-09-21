# Blind Retrospective Benchmark — pre-registered virtual screening study

> **This is a public extract** of a directory inside a private monorepo. Commit SHAs here
> are new and do **not** include the original pre-registration commit `92c6f8bf`
> (2026-07-14) — see `PROVENANCE.md`, which states plainly what can and cannot be verified
> from this repository.
>
> **There is no archived deposit.** A Zenodo DOI (`10.5281/zenodo.21824633`) was reserved
> but never minted, and has never resolved; the draft was deleted on 7 Aug 2026. Earlier
> revisions of this file linked it as an "archived deposit", which was wrong. Nothing was
> ever published there and there is nothing to cite. See Amendment 26.
>
> ---
>
> **Read this before the study below: the question this repository was built to answer was
> not answered, and the finding that replaced it was then refuted by our own later work.**
>
> **1. The study could not run as designed.** Enforced without exception, the debias gate
> excluded **every validated target assessed**. Exactly one library passed — and it was
> **not a target**. `CHEMBL612545` is a ChEMBL `UNCHECKED` record with no organism, no
> target components and 2.3 million heterogeneous activities, paired with a genuine PD-L1
> receptor structure by a hand-written label. Its "actives" are HCN1 channel blockers and
> IL-6 release inhibitors sharing no pharmacophore. Real human PD-L1 is `CHEMBL3580522`.
> That error survived twenty amendments and a public release (Amendment 21).
>
> **2. The mechanism we proposed for it is refuted, by us.** Amendment 21 claimed the
> malformed library passed *because* it was incoherent — that a set with no coherent
> property signature gives a debias classifier nothing to separate. **Amendment 23 refuted
> that** (pooled pseudo-targets median 0.6765 against real targets 0.6167, p = 0.773 — the
> pooled median is *higher*). Amendment 24's replacement, distributional identity, **was
> not supported by Amendment 25** (p = 0.711). Per the commitment recorded in A24.4, **no
> third mechanism was proposed.** If you have read the Amendment 21 framing anywhere, it is
> superseded.
>
> **3. What actually survives, and it needs no mechanism.** At the active-set sizes
> routinely used in this literature, **this class of decoy-bias metric does not support
> per-target admissibility judgements.** A library with no target relationship whatsoever
> is statistically indistinguishable from a real target under the gate (0.603 against
> 0.617, p = 0.711) and passes at the same rate (50% against 43%). **Five of seven real
> targets fall inside the null range**, which at n ≈ 30 spans 0.511–0.821 — wider than the
> entire interval between "clean" and "badly biased" as the 0.60 threshold is normally
> read. `CHEMBL612545` scored 0.4832, which against that spread is unremarkable. It did not
> pass because of any special property; verdicts at this sample size are close to
> uninformative.
>
> Secondary finding: gate verdicts depend on how many actives were sampled — SARS-CoV-2
> Mpro qualified at 30 actives (AUROC 0.543) and was excluded at 85 (0.665). Compute-driven
> active caps are near-universal in this literature and seldom reported.
>
> **Caches and results in this repository from the `CHEMBL612545` docking run are retained
> for transparency but are not evidence about anything.** That run was stopped on
> discovery. Note also that `CHEMBL4005` is recorded throughout as "aldose reductase"; it
> is **PI3Kα**.
>
> **Why this repository is public again.** It was set private on 7 Aug 2026 because the
> manuscript it supported was not being submitted (Amendment 26). Amendment 26.4 committed
> to republishing it as a condition of any future use of the work, and it is now cited as a
> worked case in a published audit checklist. The residual limitation stands and is
> restated rather than dropped: an extract cannot carry the original pre-registration
> commit, so the freeze date is attested and reviewable on request rather than
> independently verifiable from this artefact alone.

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
