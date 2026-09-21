# PRE-REGISTRATION — Blind Retrospective Benchmark (Pilot)

**Frozen: 2026-07-14.** This document is the protocol, committed to git *before* any
method is scored against ground truth. **The git commit that adds this file is the
tamper-evident timestamp** — the analysis choices below cannot be revised after
seeing results without an auditable follow-up commit. If any parameter changes, that
is a visible amendment, not a silent tune.

Scope doc: `docs/retrospective_benchmark_scope.md`.

---

## 1. Hypothesis (falsifiable)

On protein–ligand bioactivity records first published **after the model training
cutoff**, the AegisMind full pipeline (`platform` = z(−Vina) + z(ensemble-activity))
ranks true actives **(a) above chance and (b) above Vina docking alone**, by a
statistically significant margin across the pre-registered target panel.

**We commit to reporting the result whichever way it falls, including a null result.**

## 2. Primary endpoint & success criterion (frozen)

- **Primary metric:** `EF@1%` (enrichment factor in the top 1% — early recognition).
- **Primary test:** Wilcoxon signed-rank on paired per-target `EF@1%`, `platform` vs
  `vina`, across the panel.
- **SUCCESS is declared iff BOTH hold:**
  1. `platform` > `vina` on median per-target `EF@1%` **and** Wilcoxon *p* < 0.05; **and**
  2. `platform` > `chance` on the same test (sanity floor).
- **Secondary (reported, not gating):** BEDROC (α=20), EF@5%, AUROC, each with
  bootstrap 95% CIs (seed=0, 2000 resamples).

## 3. Leakage control (the linchpin)

- **Temporal split:** include only compounds whose ChEMBL source document year is
  **strictly greater than `CUTOFF_YEAR = 2025`** (→ 2026+). Rationale: the ensemble's
  training cutoffs are ~2025; post-cutoff records cannot have been memorised.
- The ensemble prompt provides structure + target context and is **never given the
  measured label**.
- Known residual risk (documented, not hidden): a post-cutoff document may re-report a
  long-known compound. Mitigated by preferring new compound–target pairs; not fully
  eliminable. This is why beating Vina — a method with *no* memory — is the real bar.

## 4. Frozen parameters (mirrored in code; code imports these constants)

| Parameter | Value | Where |
|---|---|---|
| `CUTOFF_YEAR` | 2025 (keep year > 2025) | data_split.py |
| Active threshold | pChEMBL ≥ 6.0 | data_split.py |
| Inactive threshold | pChEMBL ≤ 5.0 (5–6 grey zone dropped) | data_split.py |
| Activity types | IC50, Ki, Kd, EC50 | data_split.py |
| Panel eligibility | ≥ 10 actives AND ≥ 10 inactives (post-cutoff) | data_split.py |
| Molecule dedup | max pChEMBL per molecule | data_split.py |
| Primary metric | EF@1% | run_pilot.py / metrics.py |
| BEDROC α | 20 | metrics.py |
| Bootstrap | 2000 resamples, seed 0, 95% CI | metrics.py |
| Panel test | Wilcoxon signed-rank, α = 0.05 | metrics.py |
| Platform score | z(−Vina affinity) + z(ensemble activity 0–1) | platform_rank.py |

## 5. Target panel selection (no cherry-picking)

- Candidate targets are chosen for **data availability and 3D-structure existence only**,
  **before** the platform is run on any of them.
- A candidate **qualifies** purely by the §4 eligibility rule (≥10/≥10 post-cutoff).
- `data_split.build_panel()` freezes the qualifying set + its compounds to
  `panel_frozen.json`. The platform is **never consulted** in selection.
- **Pilot:** first 3 qualifying targets. **Full run:** 10–15 (a separate amendment
  commit once the pilot harness is validated end-to-end).

## 6. Baselines (all scored blind, identically)

1. `chance` — random scores (seed 0), the floor.
2. `fp_similarity` — max Tanimoto to **pre-cutoff** known actives (RDKit Morgan r2/2048).
   Fair: uses only pre-cutoff knowledge.
3. `vina` — −(best AutoDock Vina affinity); the method `platform` must beat.

## 7. Anti-self-deception safeguards (from the MSH3/CHDI post-mortem)

- **Target-identity verification** of every receptor against its structure header
  before docking (the exact class of error that mis-assigned eptifibatide).
- **No shared-oracle pseudo-replication:** the ensemble signal and the Vina signal are
  distinct inputs; we do not dress one oracle up as several.
- **Blind before unblind:** all rankings are written to `predictions_*.json` before any
  join to labels (`run_pilot.py` enforces the order).
- **One pre-registered analysis:** the primary metric and test above are the decision;
  secondary metrics are descriptive, not a menu to pick a flattering winner from.
- **Skipped ≠ substituted:** any method lacking inputs is reported as skipped, never
  back-filled or faked.

## 8. What a null result means

If `platform` does not beat `vina`: the platform's value is **not** point-prediction
accuracy on unseen chemistry. That is a real, reportable finding. The pivot is to
benchmark what the platform *does* do that Vina cannot (hypothesis breadth, cross-domain
reasoning, throughput, adversarial self-critique) — **not** to re-roll the ranking until
it looks good. Re-rolling to beat the metric is prohibited by this pre-registration.

---

# AMENDMENT 1 — 2026-07-14 — property-matched debiased decoys

**Visible amendment.** §1–8 above are unchanged; this section supersedes the specific
parameters it names. The commit adding this section is its tamper-evident timestamp.
Made *before* any registered (post-cutoff) result exists — it is triggered by a dry-run
data-availability finding, not by any platform score.

## A1.1 Why

The harness dry-run (real 2025 ChEMBL data) showed recent bioactivity data is **dominated
by actives** — people publish hits, not misses:

| Target | actives | inactives |
|---|---|---|
| EGFR (CHEMBL203) | 220 | 2 |
| JAK2 (CHEMBL2971) | 16 | 0 |
| CDK2 (CHEMBL301) | 22 | 2 |
| VEGFR2 (CHEMBL279) | 57 | 10 |

An enrichment benchmark needs **rare actives among many inactives** — the opposite of this.
As originally specified it (a) excludes nearly every target for lack of inactives, and (b)
on survivors tests "find actives" in an ~85%-active pool, which is close to meaningless.

## A1.2 What changes

**Negatives are now: any genuine post-cutoff measured inactives (kept) + property-matched,
topologically-dissimilar DECOYS, topped up to a fixed ratio.** Decoys are *presumed*
inactive (no measured activity against the target).

- **Leakage is unaffected where it matters.** The temporal split still governs the
  ACTIVES — the leak-critical class (the "yes, compound X hits target Y" answers an LLM
  could recall). Decoys carry **no target-specific published answer** to leak; they are
  unlabeled background. So decoys do not reintroduce the memorisation risk §3 controls.
- **Panel eligibility (supersedes §4/§5):** a target qualifies on **≥ 15 post-cutoff
  actives** alone (inactives are now supplied). This also un-blocks the targets the dry
  run excluded.

## A1.3 Decoy protocol (frozen)

For each active, draw decoys from a drug-like ChEMBL background (Ro5-compliant), such that
each decoy:
1. **Property-matches** the active within: MW ±25 Da, clogP ±0.5, HBD ±1, HBA ±1,
   rotatable bonds ±1, **formal charge exact**;
2. is **topologically dissimilar** to *every* active: max Tanimoto (Morgan r2, 2048) **< 0.35**
   — avoids latent-actives / analogue leakage (the DUD-E failure mode);
3. has **no recorded ChEMBL activity against the target** (excludes known actives).

## A1.4 Frozen amendment parameters

| Parameter | Value |
|---|---|
| `DECOY_RATIO` | 25 decoys per active |
| `MAX_ACTIVES` | 30 per target (seed-0 subsample if more — bounds docking compute) |
| `MIN_ACTIVES` | 15 (panel eligibility) |
| Property windows | MW ±25, clogP ±0.5, HBD ±1, HBA ±1, RotB ±1, charge exact |
| Decoy↔active max Tanimoto | < 0.35 (Morgan r2, 2048) |
| Background | ChEMBL drug-like (Ro5), sampled seed 0 |
| Ensemble shortlist | ensemble scores the **docking top-150** only; the rest get docking-only platform score (bounds API cost; EF@1%/BEDROC live at the sharp end anyway) |

## A1.5 Debiasing safeguard (the hostile-reviewer gate)

The classic decoy failure is that a **physicochemical-only** classifier separates
actives from decoys, so "enrichment" is a property artifact, not molecular recognition.
Guard: fit a 5-fold-CV logistic regression on [MW, clogP, HBD, HBA, RotB, TPSA, formal
charge] to separate actives from the chosen decoys.

- **Accept the decoy set only if that classifier's AUROC ≤ 0.60.**
- If > 0.60, re-draw decoys (≤ 3 attempts); if still > 0.60, **exclude the target and
  report the exclusion** — never ship a property-separable set.

Also report, per target, the **Bemis–Murcko scaffold count of the actives** (analogue-bias
transparency): a target whose actives are one chemotype is flagged, since "enrichment"
there is really "recognise one series."

## A1.6 Unchanged

Primary metric **EF@1%** (now meaningful under the controlled ratio), the paired
Wilcoxon panel test and success criterion (§2), the platform-score definition
(`z(−Vina) + z(ensemble)`, §4), the leakage temporal split for actives (§3), blind→unblind
ordering, the negative-result commitment, and every §7 anti-self-deception safeguard.

---

# AMENDMENT 2 — 2026-07-27 — ChEMBL data availability: cutoff relaxed to year > 2024

**Visible amendment.** §1–8 and Amendment 1 are unchanged except where explicitly superseded below. This commit is the tamper-evident timestamp; it is added **before** any post-2024 result is scored.

## A2.1 Why

ChEMBL's release cadence lags publication by 6–12 months. As of July 2026 the database (ChEMBL 35) contains essentially no records with `year = 2026`, making the registered `CUTOFF_YEAR = 2025` (keep year > 2025) yield zero qualifying compounds — the benchmark cannot run. This is a data-availability failure of the pre-registration, not a finding about the platform.

## A2.2 What changes

**`CUTOFF_YEAR` is amended from `2025` to `2024`** (keep `year > 2024`, i.e. 2025-dated records). The primary endpoint, target-panel construction, decoy protocol, blind→unblind ordering, and negative-result commitment are unchanged.

## A2.3 Residual leakage risk (documented, not hidden)

The ensemble models' training cutoffs are approximately August 2025. Records with `year = 2025` and published before August 2025 may appear in model training data.

**Bounding the risk:**
- The temporal split acts on ChEMBL source-document publication year, not individual compound date. A document year of 2025 does not mean the compound was in training data — only documents indexed by August 2025 would be.
- The real bar remains **beating Vina** — a method with zero memorisation capability. If the platform beats Vina using 2025 data, either (a) the platform's cross-domain reasoning is genuinely better, or (b) it memorised the answer Vina cannot. We cannot eliminate (b) for H1-2025 data.
- **Guard:** we separately report the fraction of actives whose source document was published pre- vs post-August 2025, and flag any target where >30% of actives are pre-August-2025 documents as `HIGH_LEAKAGE_RISK`.

## A2.4 Frozen amendment parameters

| Parameter | Registered value | Amended value |
|---|---|---|
| `CUTOFF_YEAR` | 2025 (year > 2025) | 2024 (year > 2024) |
| Leakage guard | — | Flag targets with >30% actives from pre-Aug-2025 docs |


---

# Amendment 3 — Subprocess path fix + pilot expanded to all 7 targets (2026-07-27)

## A3.1 Pilot run 1 findings

The first pilot run (2026-07-27) revealed two issues:

1. **Module resolution bug.** `_dock_one` (the function run inside `ProcessPoolExecutor` workers) used `Path(__file__)` to resolve the `run_vina_docking` module path. `__file__` is unreliable in spawned subprocesses; all Vina calls silently returned `No module named 'run_vina_docking'`. Zero valid Vina affinity scores were produced for CDK2. **Fix:** replaced `Path(__file__)` with the module-level `BENCH` variable (absolute path computed at import time). No statistical parameter changes; this is a bug fix.

2. **Debias gate excluded 2 of 3 pilot targets.** VEGFR2 (AUROC=0.701) and EGFR were both excluded. Kinase inhibitor actives cluster tightly in physicochemical property space, so property-matched decoys remain separable from actives even at the registered matching windows. This is correct behaviour of the gate — the test would have been artifactually easy — but it leaves too few targets for a meaningful panel test at PILOT_TARGETS=3.

## A3.2 What changes

- **Bug fix:** `_dock_one` now uses `BENCH.parent / "amr_glass" / "docking"` (absolute) for the Vina module path, and adds the repo root to `sys.path` for `aegismind` imports. No change to statistical logic.
- **`PILOT_TARGETS` expanded from 3 to 7.** The full panel is run; the debias gate still controls which targets contribute to panel tests. The primary endpoint, scoring function, panel test (Wilcoxon signed-rank), and negative-result commitment are unchanged.
- **Receptors prepared for remaining 4 targets:** JAK2 (CHEMBL2971, PDB 3LPB), BRAF (CHEMBL205, PDB 1UWH), Aurora A (CHEMBL325, PDB 2BCJ), BTK (CHEMBL1075104, PDB 3GEN).

## A3.3 Frozen amendment parameters

| Parameter | Registered value | Amended value |
|---|---|---|
| `PILOT_TARGETS` | 3 | 7 (all panel targets) |
| `_dock_one` path | `Path(__file__)` | `BENCH` (module-level absolute) |

---

# Amendment 4 — Switch to diverse non-kinase panel (2026-07-27)

## A4.1 Kinase panel finding

All seven kinase targets in the registered panel fail the debias gate (5/7) or produce EF@1%=0 for all methods including chance (2/7). The root cause is structural: recent ChEMBL kinase inhibitors are predominantly ATP-competitive analogues with tightly clustered physicochemical properties, so property-matched decoys remain separable by logistic regression despite the Tanimoto distance filter. This is a known limitation of standard virtual screening benchmarks for popular kinase targets.

## A4.2 What changes

The pilot panel is replaced with three structurally diverse non-kinase targets, each with post-2024 qualifying actives:

| Target | ChEMBL ID | PDB | Post-2024 actives | Inhibitor classes |
|--------|-----------|-----|-------------------|-------------------|
| SARS-CoV-2 Mpro | CHEMBL4523582 | 7K40 | 85 | Peptidomimetics, α-ketoamides, non-covalent, macrocycles |
| Factor Xa | CHEMBL244 | 2J34 | 44 | Diverse pharmacophore classes (S1 oxyanion/S4 aromatic) |
| BRD4 bromodomain 1 | CHEMBL1163125 | 3MXF | 159 | Thienodiazepines, natural products, diverse scaffolds |

These targets span three distinct binding-site types (cysteine protease, serine protease, bromodomain), providing a genuine test of method generality rather than a kinase-specific benchmark.

## A4.3 Unchanged parameters

All statistical parameters, scoring function, debias gate, primary endpoint (EF@1%), and negative-result commitment are unchanged. The frozen panel (`panel_frozen.json`) is committed before any platform predictions.

---

# Amendment 5 — High-diversity five-target panel (2026-07-27)

## A5.1 Amendment 4 pilot results

The Amendment 4 pilot (Mpro/FXa/BRD4) found:

| Target | Debias AUROC | Gate | Vina EF@1% | Platform EF@1% |
|--------|-------------|------|------------|----------------|
| Mpro (CHEMBL4523582) | 0.586 | PASS | 6.61 | 6.61 |
| FXa (CHEMBL244) | >0.60 | FAIL — excluded | — | — |
| BRD4 (CHEMBL1163125) | >0.60 | FAIL — excluded | — | — |

FXa and BRD4 failed the debias gate (property-separable decoys). Mpro passed and showed EF@1%=6.61 for both Vina and Platform — the ensemble contributed zero signal beyond Vina alone. This is consistent with the haiku-class model lacking binding-mode intuition.

## A5.2 What changes

The panel is expanded to five targets by adding four mechanistically diverse non-kinase targets with post-2024 ChEMBL actives:

| Target | ChEMBL ID | PDB | Post-2024 actives | Binding site type |
|--------|-----------|-----|-------------------|-------------------|
| PD-L1 | CHEMBL612545 | 5J89 | 1209 | PPI inhibitor (β-sheet surface) |
| Acetylcholinesterase | CHEMBL3717 | 4EY6 | 34 | Gorge + catalytic triad |
| Aldose reductase | CHEMBL4005 | 2IKH | 69 | NADPH/anion binding |
| MDM2 | CHEMBL3778 | 4JRG | 39 | PPI hotspot |
| SARS-CoV-2 Mpro | CHEMBL4523582 | 7K40 | 85 | Cysteine protease (anchor from A4) |

Five mechanistically distinct binding modes: PPI-surface, gorge, cofactor-site, p53-hotspot, catalytic dyad. This is the most diverse panel tested.

## A5.3 Unchanged parameters

All statistical parameters (EF@1% primary, Wilcoxon signed-rank, debias gate AUROC ≤ 0.60, 25:1 decoy ratio, temporal split CUTOFF_YEAR=2024), scoring function, and negative-result commitment are unchanged. Targets that fail the debias gate are excluded from the primary analysis exactly as in A4.

---

# Amendment 6 — Drop macromolecule actives before decoy generation (2026-07-28)

## A6.1 Problem found

PD-L1 (CHEMBL612545) has post-2024 actives spanning MW 352–3,187 Da, including antibody fragments and large macrocyclic peptides. These macromolecules:
(a) are not amenable to property-matched small-molecule decoy generation;
(b) cause ChEMBL API calls with extreme parameters (MW±25 around 3,187 Da) to either return zero results or time out after 4×60s retries, stalling the pipeline;
(c) cannot be docked with Vina (molecular weight incompatible with standard AutoDock Vina grid setup).

This is a failure of data curation, not of the benchmark design — virtual screening benchmarks are intrinsically defined over the small-molecule chemical space.

## A6.2 What changes

A drug-like mass cutoff **MW ≤ 800 Da** is applied to the active set (per target, before subsampling to MAX_ACTIVES=30) when building decoys. Actives with MW > 800 Da are dropped from the decoy-building step.

If fewer than MIN_ACTIVES (15) drug-like actives remain after this filter, the target is excluded with status `excluded_macromolecule_dominated`.

## A6.3 Rationale and scope

This filter has no effect on targets with exclusively small-molecule actives (Mpro, MDM2, Aldose reductase). It affects PD-L1 (many macrocycles/peptides in recent datasets) and potentially AChE (some carbamate macrocycles in post-2024 sets). The debias gate still runs after this filter; targets passing MW≤800 but failing the gate are still excluded exactly as before.

## A6.4 Unchanged parameters

All prior parameters are unchanged. The MW≤800 filter is a data curation step, not a statistical parameter change.

---

# Amendment 7 — Add TPSA matching to decoy property windows (2026-07-28)

## A7.1 Problem found

The Amendment 5 diverse panel (run8, 2026-07-28) shows systematic debias gate failure:

| Target | Debias AUROC | Gate |
|--------|-------------|------|
| PD-L1 (CHEMBL612545) | 0.655 | FAIL |
| Acetylcholinesterase (CHEMBL3717) | 0.681 | FAIL |
| Aldose reductase (CHEMBL4005) | 0.717 | FAIL |
| MDM2 (CHEMBL3778) | 0.856 | FAIL |

Root cause: the Amendment 1 property-matching windows (MW ±25, clogP ±0.5, HBD ±1, HBA ±1, RTB ±1, charge exact) exclude TPSA. The debias gate tests TPSA as one of its 7 features. Drug-like actives for a given target cluster in TPSA space (MDM2 inhibitors: hydrophobic pocket → low TPSA ~50–80; kinase inhibitors: ATP-pocket → similar clustering). Random ChEMBL molecules that match on the other 6 properties span a broad TPSA range. The logistic regression then trivially separates actives from decoys on TPSA alone, and the gate correctly rejects the set — but the protocol has no path to produce passing decoys without fixing the root cause.

This is discovered before any enrichment results are observed (only gate exclusion status was available). It was triggered by the observation that 4 of 4 diverse non-kinase targets now fail (not just the kinase series from Amendments 3–4), pointing to a systematic protocol gap rather than a target-class artefact.

**Precedent:** DUD-E (Mysinger et al. 2012) and DUDE-Z (Zhang et al. 2023) both match decoys on PSA/TPSA explicitly. Our protocol was less stringent than the field standard.

## A7.2 What changes

**Add TPSA ±20 Å² to the property-matching windows** for decoy generation, both server-side (ChEMBL API filter `molecule_properties__psa`) and client-side (`property_match()`). The TPSA window ±20 Å² matches the DUDE-Z standard.

No change to: the debias gate threshold (0.60), DECOY_RATIO (25), MAX_ACTIVES (30), Tanimoto cutoff (0.35), primary endpoint, scoring function, or negative-result commitment.

## A7.3 Consequence

Adding a TPSA window makes decoy selection more stringent — fewer candidates qualify. Targets that previously had marginal decoy counts may now fall below the 0.5× floor and be excluded as `excluded_insufficient_decoys`. This is acceptable: a target excluded for insufficient decoys is more honest than a target included with property-separable decoys.

## A7.4 Frozen amendment parameters

| Parameter | Registered value | Amended value |
|---|---|---|
| TPSA window | not matched | ±20 Å² (DUDE-Z standard) |
| Server filter | no PSA filter | `molecule_properties__psa` ±20 |
| `property_match` | 6 features | 7 features (adds TPSA) |

---

# Amendment 8 — Fix receptor PDBQT format (rigid vs. flexible) (2026-07-28)

## A8.1 Problem found

Run8 (2026-07-28) revealed that all Vina docking calls produced `WARNING: no affinities parsed`. Direct inspection showed Vina returning exit code 1 with:

```
PDBQT parsing error: Unknown or inappropriate tag found in rigid receptor.
 > ROOT
```

Root cause: `prepare_receptors.py` used `obabel ... -O receptor.pdbqt` to convert the protein PDB to PDBQT. obabel treats the protein as a *flexible molecule* and inserts `ROOT`/`ENDROOT`/`BRANCH`/`ENDBRANCH`/`TORSDOF` ligand-tree markers (1383 active torsions for a 300-residue protein). AutoDock Vina 1.2's `--receptor` expects a *rigid* PDBQT with no tree markup — any such tag causes a fatal parse error and no poses are produced.

This bug was present in all runs after the initial harness setup. It explains why run8 (and retrospectively the Amendment 4/5 docking steps before) produced zero valid affinity scores. The `EF@1%=6.61` reported for Mpro in the Amendment 4 summary was from a run that either used a pre-existing (correctly prepared) PDBQT or had the bug silently masked by fallback scoring — this run (run9) is the first clean end-to-end test.

## A8.2 What changes

`prepare_receptors.py` now calls `pdb_to_pdbqt_receptor()` from `amr_glass/docking/run_vina_docking.py` (the same manual AD4-atom-type writer already used for the AMR docking work). This produces a plain ATOM-record PDBQT with no tree markup — exactly what AutoDock Vina 1.2 requires.

All receptor PDBQT files regenerated. `receptors.json` updated with absolute paths (using `Path(__file__).resolve().parent`) so subprocess workers inherit correct paths regardless of CWD.

**Validated:** nirmatrelvir docked to Mpro (7K40, box 22 Å centred at [9.9, −31.26, 19.56]) returns best pose −6.68 kcal/mol with 5 poses — consistent with literature values for this inhibitor class.

## A8.3 Unchanged parameters

All statistical parameters, debias gate, primary endpoint, scoring function, and negative-result commitment are unchanged. Only the receptor file preparation method changed.

---

# Amendment 9 — Remove server-side PSA filter (ChEMBL/RDKit inconsistency) (2026-07-28)

## A9.1 Problem found

Run9 (2026-07-28, first clean end-to-end run post Amendment 8 receptor fix) tested
Amendment 7's TPSA matching on the Amendment 5 diverse panel. Results:

| Target | Run8 AUROC (no TPSA) | Run9 AUROC (TPSA Amend.7) | Direction |
|--------|----------------------|---------------------------|-----------|
| PD-L1 | 0.655 | 0.700 | **worse** |
| AChE | 0.681 | 0.819 | **worse** |
| Aldose reductase | 0.717 | 0.625 | better |
| MDM2 | 0.856 | 0.904 | **worse** |

The client-side RDKit TPSA matching in `property_match()` (Amendment 7) should help
by tightening the decoy pool. However, Amendment 7 also added a server-side PSA filter
(`molecule_properties__psa__gte/lte`) to `_fetch_candidates()`. The root cause of the
worsening is the **inconsistency between ChEMBL's server-side PSA and RDKit's
`Descriptors.TPSA`**: these values differ because ChEMBL uses a different partial
charge model and atom-type parameterisation. Adding a ChEMBL-PSA server filter biases
the candidate pool toward molecules where ChEMBL-PSA ≈ active's RDKit-TPSA, but this
cluster happens to be more separable from actives on other physicochemical features,
causing the gate to AUROC to increase rather than decrease.

## A9.2 What changes

Remove the two PSA server-filter lines from `_fetch_candidates()`:
```python
# REMOVED (Amendment 9):
# "molecule_properties__psa__gte": round(a["tpsa"] - TPSA_WIN, 1),
# "molecule_properties__psa__lte": round(a["tpsa"] + TPSA_WIN, 1),
```

The client-side RDKit TPSA check in `property_match()` (Amendment 7) is **retained**:
```python
and abs(a["tpsa"] - c["tpsa"]) <= TPSA_WIN   # consistent RDKit↔RDKit comparison
```

This preserves TPSA matching with a consistent computation on both sides, while
avoiding the server-side bias from the ChEMBL/RDKit PSA inconsistency.

## A9.3 Trigger

Observed before any enrichment results were inspected (only gate exclusion AUROC values
were examined). Discovered systematically from the pattern across 4 targets in run9.

## A9.4 Unchanged parameters

All statistical parameters, debias gate threshold (≤0.60), primary endpoint, scoring
function, and negative-result commitment are unchanged. Only the server-side candidate
fetch filter changes.

---

# Amendment 10 — Debias gate interpretation and future panel expansion (2026-07-28)

## A10.1 Finding

Across all runs (run8, run9, run10), the debias gate persistently excludes PD-L1, AChE,
Aldose reductase, and MDM2. The AUROC range 0.625–0.904 is well above the 0.60 ceiling.

The gate is operating CORRECTLY. These targets have tightly clustered active chemotypes:
- **MDM2** inhibitors: flat hydrophobic PPI disruptors (TPSA 50–80, MW 350–500)
- **AChE** inhibitors: aromatic amine scaffolds with charged nitrogen
- **PD-L1** binders: linear peptide-mimetic biphenyl scaffolds (TPSA 100–130)
- **Aldose reductase** inhibitors: carboxylate/tetrazole warhead scaffolds (TPSA 90–120)

Property-matched decoys drawn from the ChEMBL general drug-like background (which spans
all therapeutic areas) are systematically distinguishable from these target-specific
active sets even after MW/clogP/HBD/HBA/RTB/charge/TPSA matching. This is the
"analogue bias" problem documented in DUDE-Z (Stein et al. 2021).

## A10.2 Implication

Only Mpro (CHEMBL4523582) provides valid enrichment data from the Amendment 5 pilot panel.
Mpro actives are chemotypically diverse (covalent warheads, peptidomimetics, macrocycles,
small fragments) because the COVID pandemic generated a large, chemically diverse set of
post-cutoff data from multiple independent screening campaigns.

## A10.3 Not a protocol failure

The debias gate is working as designed: it prevents reporting spurious enrichment on a
biased benchmark. The correct action is to report Mpro EF@1% as the sole pilot result
and note the four excluded targets. This is scientifically valid.

## A10.4 Future panel amendment (not yet implemented — requires new data_split run)

A future panel should add targets with known chemotypic diversity in post-2025 ChEMBL data:
- Other SARS-CoV-2 proteins (PLpro, Nsp15) — active COVID research pipeline
- Emerging targets with broad scaffold diversity (e.g., recently deorphaned GPCRs)
- Targets from recent phenotypic screens where scaffold diversity is built-in

Selection criterion: debias gate AUROC ≤ 0.60 on a trial decoy build. Only targets that
pass the trial build should enter the registered panel. This "trial gate" prevents the
systematic exclusion seen in the Amendment 5 panel.

## A10.5 Unchanged parameters

All statistical parameters, primary endpoint, and negative-result commitment are unchanged.
This amendment documents the pilot findings and does not alter any methodology for the
current runs.

---

# Amendment 11 — Deterministic ChEMBL candidate ordering (2026-07-28)

## A11.1 Problem found

Between Mpro run8 (EF@1%=6.61) and Mpro run10 (EF@1%=0.00), the same target with the
same SEED produced a completely different decoy set, causing EF@1% to swing from 6.61
to 0.  Root cause: `_fetch_candidates()` calls the ChEMBL REST API with `limit=500` but
no `order_by`, so the server returns an arbitrary 500-molecule subset on each call. Even
with a fixed `rng` seed, different candidate pools lead to different decoy sets.

## A11.2 What changes

Add `"order_by": "molecule_chembl_id"` to the `_fetch_candidates()` API call.
ChEMBL guarantees stable alphabetical ordering of molecule IDs across calls.
This makes the candidate pool deterministic (for the same active property window),
so the decoy set depends only on the rng seed and not on transient API state.

```python
"order_by": "molecule_chembl_id",   # Amendment 11: deterministic ordering
```

## A11.3 Affected runs

- Runs ≤ run10 (Mpro): NOT retroactively fixed. The Mpro run10 result (EF@1%=0.00,
  pre-registered null) stands.
- Run11 (PD-L1, launched before Amendment 11 was committed): NOT retroactively fixed.
  Run11 will use the old non-deterministic code. Its result stands.
- Runs ≥ run12: will use the deterministic ordering.

## A11.4 Unchanged parameters

All statistical parameters, primary endpoint, scoring function, and the pre-registered
negative-result commitment are unchanged.

---

# Amendment 12 — PD-L1 gate failure with deterministic decoys (2026-07-28)

## A12.1 Finding

Run12 (2026-07-28, Amendment 11 deterministic decoys) re-tested PD-L1 (CHEMBL612545).
Gate AUROC = 0.620 → EXCLUDED (above the 0.60 threshold).

In run10-others (non-deterministic candidates): AUROC = 0.586 → passed.
In run12 (deterministic, order_by=molecule_chembl_id): AUROC = 0.620 → excluded.

## A12.2 Interpretation

PD-L1 is borderline. Its gate AUROC is sensitive to which specific ChEMBL molecules
are selected as decoy candidates. With alphabetical ordering (deterministic), the
candidate pool is slightly more property-separable from PD-L1 actives than with a
random pool. AUROC 0.620 vs 0.586 is within sampling error, but the deterministic
result (0.620) is the correct one to use — it is reproducible and not dependent on
API call order.

**Conclusion:** PD-L1 is EXCLUDED from the registered panel. The earlier 0.586
pass was an artifact of random candidate ordering. This is not a protocol failure;
it is the gate working correctly on a borderline target.

## A12.3 What changes

PD-L1 is removed from the qualifying target list. Mpro (CHEMBL4523582) remains as
the only qualifying target from the Amendment 5 panel (gate AUROC ≤ 0.60 confirmed
in run9 at 0.538 and run10 at 0.543).

Future panel candidates: Aldose reductase (CHEMBL4005, run9 gate AUROC 0.625 without
server PSA filter) will be retested with Amendment 11 deterministic code.

## A12.4 Unchanged parameters

All statistical parameters, primary endpoint, and negative-result commitment unchanged.

---

# Amendment 13 — Aldose reductase gate failure (2026-07-28)

## A13.1 Finding

Run12-Aldose (2026-07-28, Amendment 11 deterministic code): Aldose reductase
(CHEMBL4005) gate AUROC = 0.750 → EXCLUDED.

History:
| Code | Gate AUROC | Direction |
|------|-----------|-----------|
| Run8 (MW/clogP/HBD/HBA/RTB only) | 0.717 | baseline |
| Run9 (+ server PSA + client TPSA) | 0.625 | BETTER |
| Run12 (client TPSA only, deterministic) | 0.750 | WORSE than baseline |

The server PSA filter (removed in Amendment 9 because it worsened PD-L1/AChE/MDM2)
was specifically helpful for Aldose reductase. Removing it pushed the gate AUROC
above the baseline (0.717 → 0.750), likely because the alphabetical ordering
(Amendment 11) selects a different subset of ChEMBL candidates that is more
property-separable.

## A13.2 Conclusion: Mpro is the sole qualifying target

After testing 12 targets across 3 panels:
- Kinase panel (8 targets): all excluded (AUROC 0.63–0.93)
- Amendment 4 panel (FXa, BRD4, Mpro): FXa excluded, BRD4 excluded, Mpro QUALIFIED
- Amendment 5 panel (PD-L1, AChE, MDM2, Aldose): all excluded (AUROC 0.62–0.90)

Only SARS-CoV-2 Mpro (CHEMBL4523582, gate AUROC ~0.54) passes the debias gate.
Mpro's post-cutoff active set is chemotypically diverse due to the COVID-19 pandemic
generating compounds from multiple independent discovery campaigns (covalent warheads,
peptidomimetics, macrocycles, fragments). No other tested target achieves this diversity.

## A13.3 Benchmark result: pre-registered null (Mpro pilot, exhaustiveness=4)

With only Mpro qualifying, the Wilcoxon signed-rank panel test cannot be performed
(requires ≥2 targets; minimum n=2 gives non-significant p=0.5). The pilot result for
Mpro (run10, exhaustiveness=4) is:
- Platform EF@1% = 0.00, Vina EF@1% = 0.00 (null — no active in top 1% of 780 cpds)
- Platform AUROC = 0.617 vs Vina AUROC = 0.566 (platform has better ranking, though both
  are only slightly above chance)
- FP-similarity baseline EF@1% = 23.1 (actives ARE chemically distinct from decoys —
  the issue is Vina's insufficient precision at exhaustiveness=4 for early recognition)

**This is a pre-registered null result and is committed to being published.**

## A13.4 Next steps

Higher-exhaustiveness Mpro run: pre-registration contemplated exhaustiveness=16 as the
production run. At exhaustiveness=16 with more docking precision, Vina may achieve non-
zero EF@1% for Mpro. This would not change the panel test conclusion (still single target)
but would provide a more meaningful single-target comparison.

Additional COVID targets: ChEMBL lists only 3 SARS-CoV-2 single-protein targets
(Replicase polyprotein 1ab CHEMBL4523582, Spike CHEMBL4662936, Nucleoprotein CHEMBL5169223).
PLpro (NSP3) and NSP15 do not have independent single-protein target entries. Spike and
Nucleoprotein have 0 post-2024 small-molecule actives (pChEMBL ≥ 6.0).

Additional non-COVID diverse targets: extended search for targets with wide post-2024
property distributions — see Amendments 14 and beyond.

---

# Amendment 14 — TCPTP/PTPN2 gate failure (2026-07-28)

## A14.1 Finding

TCPTP (tyrosine-protein phosphatase non-receptor type 2, CHEMBL3807) was identified
as a potentially qualifying target due to its wide post-2024 property profile:
- 77 total post-2024 actives (pChEMBL ≥ 6.0); 75 unique small molecules (MW ≤ 800)
- clogP range: -0.6 to 7.6 (sd=2.0), MW range: 257-676 (sd=98), scaffolds: 33/46 (72%)
- Property distribution superficially similar to Mpro (wide clogP, wide MW)

Gate result: **EXCLUDED, AUROC = 0.712** (threshold ≤ 0.60).

Despite the wide property profile, TCPTP inhibitors have a distinctive pharmacophore
signature (charged phosphate mimetics for active-site inhibitors; more lipophilic compounds
for allosteric inhibitors; the bimodal distribution is itself a separable signal). The
72% scaffold diversity is high, but scaffold diversity is not equivalent to
physicochemical property overlap with ChEMBL background.

## A14.2 Interpretation

The debias gate correctly identifies that TCPTP actives are property-separable from
background drug-like molecules. The gate is functioning as designed. The wide clogP
range is dominated by two subpopulations (polar phosphomimetics and lipophilic allosteric
compounds), and a logistic regression reliably identifies this bimodal pattern.

## A14.3 Updated target count

After 13 targets across 4 panels (including TCPTP): only Mpro qualifies. Gate correctly
identifies structural specificity in all tested non-Mpro targets.

---

# Amendment 15 — KRAS gate failure (2026-07-28)

## A15.1 Finding

KRAS (GTPase KRas, CHEMBL2189121) was evaluated as a candidate with 119 post-2024 actives.
After deduplication: 79 unique molecules; after subsampling to MAX_ACTIVES=30: 29 actives
with 20 distinct scaffolds (69% scaffold diversity ratio — comparable to Mpro).

Property profile: MW 431-733 (mean 629, sd=67), clogP 3.8-8.2 (sd=1.0).
High clogP clustering: KRAS inhibitors (G12C covalent/noncovalent, G12D inhibitors) are
systematically more lipophilic than the ChEMBL drug-like background at similar MW.

Gate result: **EXCLUDED, AUROC = 0.664** (threshold ≤ 0.60).

## A15.2 Interpretation

KRAS inhibitors cluster in high-clogP, high-MW space that is distinguishable from
property-matched background. Despite having 20 diverse scaffolds, the lipophilicity
signature is systematic enough for the logistic regression to achieve > 0.60 AUROC.
This is consistent with KRAS inhibitors being a pharmacologically specific, lipophilic
compound class that is rare in the general ChEMBL background at these MWs.

## A15.3 Final target count

After 14 targets tested (TCPTP + KRAS extending Amendment 13 panel to 14 targets):
only Mpro qualifies. The debias gate's 0.60 AUROC threshold enforces a stringent
requirement: post-cutoff active sets must have property distributions genuinely
overlapping with the ChEMBL background — not just wide ranges within a specific
pharmacological class.

Mpro passes because COVID drug discovery engaged orthogonal medicinal chemistry
strategies (covalent warheads, peptidomimetics, macrocycles, fragments, repurposed drugs)
that collectively span a property space indistinguishable from background ChEMBL molecules.

---

# Amendment 16 — production re-run of the qualifying target (2026-08-04)

**Registered before execution. No results from this configuration have been observed at
the time of writing.**

## A16.1 What changes

| Parameter | Registered / run10 | Amendment 16 | File |
|---|---|---|---|
| `EXHAUSTIVENESS` | 4 (run10; registered default 16) | **32** | `run_pilot_fast.py` |
| `MAX_ACTIVES` | 30 (frozen, Amendment 1 §A1.4) | **100** | `decoys.py` |

`DECOY_RATIO` (25:1), all property windows, `MAX_SIM`, `DEBIAS_AUROC_MAX` (0.60), the
temporal split and the primary endpoint are **unchanged**.

## A16.2 Why

Both changed constants were compute-bounding choices of ours, not properties of the
methods under test, and between them they made the pre-registered primary endpoint unable
to discriminate.

**Exhaustiveness.** run10 used 4 against a registered default of 16, for throughput. Any
null attributed to docking at exhaustiveness 4 is open to the objection that it reflects
the setting rather than the method. 32 is above the registered default, so this cannot be
read as tuning down to obtain a null.

**Active cap.** `MAX_ACTIVES = 30` discarded 55 of the 85 qualifying Mpro actives. Because
`DECOY_RATIO` is fixed, the cap does **not** change prevalence (3.78% → 3.82%); what it
changes is the absolute expected count of actives in the top slice:

| Configuration | N | EF@1% slice | Expected actives by chance @1% | @5% |
|---|---|---|---|---|
| run10 (30 actives) | 793 | 7.93 | **0.30** | 1.50 |
| Amendment 16 (85 actives) | 2223 | 22.23 | **0.85** | 4.25 |

**This does not fully repair EF@1%.** The expected count by chance remains below 1, so a
zero is still a plausible chance outcome and EF@1% remains weakly powered. EF@5%, at 4.25
expected, is the endpoint with genuine power under this configuration. We state this now,
before seeing the result, so that the interpretation cannot be chosen afterwards.

## A16.3 Known asymmetry left unchanged

LLM activity scoring covers only the top `top_n = 150` compounds by Vina rank; the
remainder receive a constant 0.5. At N = 793 that was the top 18.9%; at N = 2223 it is the
top 6.7%. The reported endpoints (top 1% = 22, top 5% = 111 compounds) both fall inside
the scored region in both runs, so EF is unaffected, but **platform AUROC is computed over
a library in which most compounds carry a constant score and is not comparable between
runs.** `top_n` is left unchanged so the platform definition is stable; AUROC is
interpreted accordingly.

## A16.4 Status of run10

run10 is **not superseded and not withdrawn.** It is reported alongside this run as the
low-exhaustiveness, capped-active configuration. Both appear in any publication.

## A16.5 Pre-commitment

The outcome is published regardless of direction. Specifically, and stated in advance:

- If enrichment appears at exhaustiveness 32 that was absent at 4, the conclusion is that
  the run10 null reflected the exhaustiveness setting, and we say so.
- If the platform again fails to exceed Vina at EF@5%, that stands as a second observation
  of the same effect.
- If the fingerprint baseline again dominates both, that is reported as the headline.

Nothing in the analysis plan is contingent on which of these occurs.

---

# Amendment 17 — the qualifying target does not qualify at full active count (2026-08-04)

**Discovered while executing Amendment 16, before any docking was run.**

## A17.1 Result

Rebuilding decoys for SARS-CoV-2 Mpro (CHEMBL4523582) with the Amendment 16 cap
(`MAX_ACTIVES = 100`, i.e. all 85 qualifying actives):

| Configuration | n actives | Scaffolds | Gate AUROC | Registered threshold | Outcome |
|---|---|---|---|---|---|
| run10 (Amendment 1 cap) | 30 (subsampled, seed 0) | 20 | **0.543** | ≤ 0.60 | Qualified |
| Amendment 16 (full set) | 85 | 44 | **0.665** | ≤ 0.60 | **Excluded** |

Status returned: `excluded_property_separable`.

## A17.2 Consequence

**The panel now has no qualifying target. 14 of 14 tested targets fail the debias gate.**

Mpro qualified in run10 only as an artefact of subsampling. The 30-active subsample drawn
under seed 0 happened to be property-matchable against the ChEMBL background; the full
85-active set is not. The full set is the more faithful characterisation of the target,
since subsampling discards information rather than adding it.

## A17.3 What is NOT being done

The gate threshold is 0.60 and was registered on 2026-07-14. Mpro at full active count
scores 0.665 and is excluded. We are **not**:

- lowering the threshold,
- searching seeds for a subsample that passes,
- reporting the n = 30 result as the target's gate value without this context, or
- selecting an intermediate cap that restores qualification.

Any of these would convert a pre-registered gate into a post-hoc filter, which is the
practice this benchmark exists to test. The Amendment 16 docking re-run is **cancelled**:
there is no qualifying target to run it on.

## A17.4 Standing of the run10 enrichment results

run10's enrichment figures are retained and reported, now correctly qualified: they were
measured on a target that passes the gate **only under a 30-active subsample**, and are
therefore not evidence about performance on an unbiased library. They are reported as a
description of what was run, not as a validation result.

## A17.5 Scientific finding

Two findings replace the previous framing:

1. **Property-matched decoy generation failed to produce an unbiased library for any of
   14 targets** spanning kinases, proteases, phosphatases, bromodomains, coagulation
   factors, immune checkpoints and oncogenes.
2. **Debias-gate outcomes are sensitive to active-set subsampling.** The same target
   passes at n = 30 and fails at n = 85. Any benchmark that subsamples actives for compute
   reasons may qualify targets that a full active set would exclude — and subsampling is
   near-universal in this literature.

Finding 2 was not sought. It emerged from a protocol change registered for an unrelated
reason (restoring power to the primary endpoint) and is reported as encountered.

---

# Amendment 18 — PD-L1 qualifies at full active count; docking proceeds (2026-08-04)

**Registered before execution. No enrichment result for PD-L1 has been observed.**

## A18.1 Re-gate result at full active count

Amendment 16 raised `MAX_ACTIVES` to 100. Re-running the registered gate at that cap for
every target with local panel data:

| Target | n actives | Gate AUROC | Verdict at full count | Verdict at n ≤ 30 |
|---|---|---|---|---|
| CHEMBL3778 | 39 | 0.9026 | Excluded | Excluded |
| CHEMBL3717 | 34 | 0.7608 | Excluded | Excluded |
| CHEMBL244 | 44 | 0.7081 | Excluded | Excluded |
| Mpro (CHEMBL4523582) | 85 | 0.6653 | Excluded | **Qualified (0.543)** |
| CHEMBL1163125 | 87 | 0.6641 | Excluded | Excluded |
| Aldose reductase (CHEMBL4005) | 48 | 0.6510 | Excluded | Excluded (0.750) |
| **PD-L1 (CHEMBL612545)** | **82** | **0.5850** | **QUALIFIED** | **Excluded (0.620, Amend. 12)** |

## A18.2 Correction to Amendment 17

Amendment 17 stated that no target qualifies and that the 14-of-14 result was
*conservative*, reasoning from Mpro that subsampling makes qualification easier. **Both
claims are withdrawn.**

Subsampling moves gate verdicts in **both** directions:

- Mpro: qualified at n = 30 (0.543), excluded at n = 85 (0.665).
- PD-L1: excluded at n = 30 (0.620), **qualified** at n = 82 (0.585).

The bias is not systematic and its direction cannot be anticipated. This is a stronger
methodological finding than either previous framing: whether a target is admissible to a
benchmark depends on a compute-driven sampling choice in an unpredictable direction.

## A18.3 Consequence

The panel is **not** empty. PD-L1 qualifies at full active count under the registered gate,
so the enrichment arm cancelled in Amendment 17 is reinstated on PD-L1.

## A18.4 Docking run, registered in advance

- Target: PD-L1 (CHEMBL612545), 82 actives, 25:1 synthetic decoys, plus real measured
  inactives from the same temporal split.
- `EXHAUSTIVENESS = 32` (Amendment 16; above the registered default of 16, so this cannot
  be read as tuning down to obtain a null).
- Arms unchanged: chance, AutoDock Vina, platform (Vina + LLM activity likelihood), 2D
  fingerprint similarity.
- Endpoints unchanged: EF@1% primary; EF@5%, BEDROC, AUROC secondary.

**Power, stated before the result.** At 82 actives the top 1% slice is ~25 compounds and
the expected count by chance is ~0.82 — still below 1, so EF@1% remains weakly powered and
a zero remains a plausible chance outcome. EF@5% (expected ~4.1) is the endpoint with
genuine power. This is stated now so the interpretation cannot be selected afterwards.

## A18.5 Pre-commitment

Published regardless of direction, specifically:

- If the platform exceeds Vina at EF@5%, that is reported as a positive result.
- If the platform again lands at or below Vina, that is a second independent observation
  of the effect first seen on Mpro.
- If the 2D fingerprint baseline again dominates both structure-based arms, that is
  reported as the headline.
- If Vina at exhaustiveness 32 enriches where exhaustiveness 4 did not, the earlier null is
  attributed to the setting and we say so plainly.

---

# Amendment 19 — all prior docking was rigid-ligand; both targets re-run (2026-08-04)

**Registered before execution. This amendment invalidates prior enrichment results.**

## A19.1 The defect

Ligand preparation calls Meeko and, on failure, falls back to `_manual_pdbqt()` in
`analysis/amr_glass/docking/run_vina_docking.py`. That writer emits every atom inside a
single `ROOT` block, no `BRANCH`/`ENDBRANCH` records, and `TORSDOF 0`.

In the environment used for all runs to date (Python 3.7), **Meeko fails on every ligand**
with `'type' object is not subscriptable` — a Python-version incompatibility. Observed
failure rate in the PD-L1 run: **87 of 87 ligands, 100%**. Every molecule was therefore
docked as a rigid body, frozen in a single arbitrary RDKit conformer, with zero torsional
degrees of freedom.

The PD-L1 ligand set has a median of 5 rotatable bonds (mean 9, max 99); **100% have at
least one.** All of that conformational freedom was discarded before docking began.

## A19.2 What this invalidates

- **All enrichment results reported to date**, on Mpro (run10) and PD-L1 alike. Both used
  the same code path.
- **The Amendment 16 rationale.** Exhaustiveness governs the thoroughness of the
  conformational and positional search. With `TORSDOF 0` there is no conformational space
  to search, so raising exhaustiveness from 4 to 32 could not have had the effect the
  amendment was written to test.
- **Any inference that docking does not enrich.** Vina at AUROC 0.566 against a 2D
  fingerprint baseline at 0.890 is the expected signature of a structure-based arm unable to
  sample ligand conformations. It is not evidence about docking. We have not yet tested
  docking on these panels.

The debias-gate results (§3.1) and the subsample-instability finding (§3.2, Amendments
16–18) are **unaffected** — neither involves docking.

## A19.3 Fix

Execute under the `docking` conda environment: Python 3.10.20, Meeko 0.7.1, RDKit 2026.03.3,
Vina 1.2.7. Verified on a 7-rotatable-bond probe: Meeko succeeds and emits 7 `BRANCH`
records with `TORSDOF 7`, against `TORSDOF 0` under Python 3.7.

**Ligand preparation is now a hard gate, not a fallback.** Any ligand for which Meeko fails
is excluded and counted, rather than silently docked rigid. The fallback writer is retained
only for receptors, where rigidity is correct.

## A19.4 Re-run

Both qualifying-at-some-active-count targets, exhaustiveness 32, arms and endpoints
unchanged:

- **PD-L1 (CHEMBL612545)** — qualifies at full active count, gate AUROC 0.585.
- **Mpro (CHEMBL4523582)** — qualifies only under 30-active subsampling; re-run so the
  subsample-qualified case is also reported on correctly prepared ligands.

## A19.5 Pre-commitment

Unchanged and restated: the outcome is published in whichever direction it falls. If
enrichment appears with flexible ligands where it was absent with rigid ones, the earlier
null is attributed to ligand preparation and we say so plainly. If it does not appear, that
is then a result about docking rather than about our pipeline.

## A19.6 Process note

This defect was visible in the run log from the first ligand — `Meeko failed ... falling
back to manual PDBQT writer` — and was not acted on until five hours of compute had been
spent. Amendment 16 was written and committed on a rationale that could not have held.
Recorded because a pre-registration that omits its own process failures is worth less than
one that includes them.

---

# Amendment 20 — distributing the Amendment 19 re-run across environments (2026-08-05)

**Written and committed before any shard is dispatched.**

## A20.1 Why

The Amendment 19 re-run is executing on a dual-core i7-6500U at 8.85 compounds/hour.
PD-L1 alone needs ~10.8 days and Mpro is queued behind it — ~22 days serial. That
delay has no scientific content; it is a hardware limit. Kaggle permits two
concurrent CPU sessions, so the remaining compounds are split between the local
process and two remote workers.

**Nothing about the protocol changes.** Exhaustiveness stays 32, the box, receptor,
scoring function and endpoints are unchanged. The only change is which machine
executes a given `vina` invocation.

## A20.2 The risk this creates, and how it is closed

Distributing work across environments risks exactly the failure Amendment 19
documents: a difference in ligand preparation that is invisible in the results.
Meeko and RDKit build the torsion tree, and an RDKit upgrade alone has already moved
gate AUROC in this project (0.585 → 0.598).

The design removes the risk rather than guarding it:

- **Ligand PDBQTs are prepared locally and shipped.** The remote side has no RDKit,
  no Meeko, no OpenBabel. It runs `vina` and nothing else.
- **The Vina binary and receptor are shipped byte-for-byte** and verified by md5 on
  arrival.
- **The protocol record stores md5s, not version strings.** A version string is a
  claim about a binary; an md5 is a fact about it. This was not academic:
  `run_vina_docking.py` recorded `"vina_version": "1.2.5"` while the binary in use
  was **1.2.7**.

## A20.3 The compound list is frozen, not regenerated

Decoys are drawn live from the ChEMBL REST API. Amendment 11 made those responses
deterministically *ordered*, which is not the same as *stable* — ChEMBL is a living
database. A remote worker regenerating the list could legitimately receive a
different decoy set, and the two environments would then be docking two different
benchmarks.

The list is therefore frozen once (`kaggle/freeze_panel.py`), md5-stamped, and
shipped. For PD-L1, which was already part-way through docking when this amendment
was written, the frozen list is **verified against the live run**: every compound
already docked must appear in it. If any is absent, ChEMBL has moved and the target
is not shardable — the script exits non-zero and no shard is dispatched.

## A20.4 Cross-environment calibration

Each shard re-docks five compounds already scored locally. Vina is run without
`--seed`, so these will not agree exactly; the acceptance criterion is
**|local − remote| ≤ 0.3 kcal/mol** on every calibration compound. A wider spread
indicates the environments differ by more than the sampler's RNG, and the shard is
rejected in full rather than partially merged.

This is a pre-committed threshold. If the calibration fails, the distributed plan is
abandoned and the run completes locally, however long that takes.

## A20.5 Stopping rule

The local process predates the resume fix and will work through all 2,526 compounds
regardless of what the remote workers return. Local takes the head of the frozen
list, remote workers the tail. The run is complete when the **union** covers the
frozen list, not when either side individually does; `kaggle/check_coverage.py`
evaluates this and is the only authority for stopping the local process. Compounds
docked twice are discarded in favour of the local score.

## A20.6 What is not claimed

Distributing compute is a logistics change and is reported as such. It does not
strengthen the result, and no enrichment claim rests on it. Recorded here only
because the benchmark's value depends on every deviation being visible.

## A20.7 Calibration result (6 Aug 2026)

The pre-committed cross-environment check in A20.4 **passes**. Five compounds already
scored locally were re-docked on Kaggle from the shipped bundle:

| Compound | Local | Kaggle | \|Δ\| |
|---|---|---|---|
| CHEMBL1169441 | −6.06 | −6.05 | 0.01 |
| CHEMBL3236427 | −5.28 | −5.34 | 0.06 |
| CHEMBL108628 | −4.92 | −4.91 | 0.01 |
| CHEMBL105503 | −5.48 | −5.48 | 0.00 |
| CHEMBL6159690 | −4.78 | −4.82 | 0.04 |

Worst 0.06, mean 0.02, against the pre-committed tolerance of 0.30 kcal/mol.

This is tighter than anticipated. Vina is run without `--seed`, so scatter of 0.1–0.2
from the sampler alone would have been unremarkable; agreement at 0.00–0.06 indicates
the two environments are performing near-identical searches. That is the intended
consequence of shipping the binary, receptor and pre-built ligand PDBQTs rather than
reproducing any of them remotely — the version surface was eliminated, not narrowed.

Two limits on what this establishes. It demonstrates **environment equivalence only**,
and says nothing about whether the docking enriches — that remains the open question the
benchmark exists to answer. And the five controls span −4.8 to −6.1 kcal/mol; agreement
there does not strictly extend to strong binders, though no mechanism is proposed by
which it would differ.

Remote results are therefore admissible. Calibration replants are discarded in favour of
the local score for the same compound. First merge: 130 remote results accepted,
bringing the PD-L1 cache from 416 to 546 of 2,526.

**Process note.** The controls were initially appended to the END of each shard, where
they would have been docked after roughly 60 hours — defeating the purpose of an early
warning. `merge_shards.py` refused both shards with "no calibration overlap", which is
the gate behaving correctly, but it should not have been reachable. Ordering corrected
before any remote result was merged.

---

# Amendment 21 — the qualifying target is not a target (2026-08-07)

**Discovered while writing §3.1/§3.2 of the manuscript, after the repository had been made
public and 1,145 of 2,526 compounds had been docked. All runs stopped on discovery.**

## A21.1 What was found

`prepare_receptors.py` pairs PDB **5J89** — a genuine PD-L1 structure — with ChEMBL target
**CHEMBL612545**, labelled `target_name: "PD-L1"`. The receptor is PD-L1. The ChEMBL
identifier is not.

CHEMBL612545 is a ChEMBL **`UNCHECKED`** record: `pref_name` "Unchecked", no organism,
**zero target components**, and **2,317,536 activities** spanning unrelated assay types.
It is a catch-all bucket for bioactivities with no validated target assignment.

Human PD-L1 is **CHEMBL3580522**.

## A21.2 Evidence

The panel's stored `pchembl` values match records in the unchecked bucket exactly:

| Panel record | Assay actually matched in CHEMBL612545 |
|---|---|
| CHEMBL6195028, 6.51, IC50 | Inhibition of human HCN1/PEX5L with TRIP8b, HEK293 |
| CHEMBL5569030, 7.54, IC50 | Inhibition of LPS-stimulated IL-6 release, human PBMC |
| CHEMBL6192612, 7.78, IC50 | Inhibition of LPS-stimulated IL-6 release, human PBMC |

Of 12 sampled "actives", **1** has any recorded activity against CHEMBL3580522. The set
is composed largely of HCN1 channel blockers and IL-6 release inhibitors.

An initial reading of this attributed one active to a metabolic-stability assay; that came
from an unfiltered query and is corrected here.

## A21.3 What it invalidates

- PD-L1's gate results at both active counts (0.620 at n=21; 0.585 at n=82).
- The claim that gate verdicts flip **in both directions**. Only the Mpro reversal
  survives — qualified at 30 actives (0.543), excluded at 85 (0.665) — and it is a real
  target, unaffected by this.
- The claim that **one** target qualifies at full active count. The correct figure is
  **zero**; the qualification was an artefact.
- Amendment 18, and the entire §3.3 enrichment run authorised by it.

## A21.4 Why it passed the gate — the finding that replaces it

The debias gate excludes a target when its actives are separable from property-matched
decoys. Thirteen coherent target sets were excluded because their actives cluster in
property space in target-specific ways.

An incoherent set — HCN1 blockers, IL-6 inhibitors and unrelated chemistry sharing no
pharmacophore — has **no coherent property signature**, so the logistic regression cannot
separate it from property-matched decoys. It passed the gate *because* it is not a target.

**A debias gate cannot distinguish "property-diverse actives" from "not a target at all",
and the second case passes most easily.** A target-coherence check therefore belongs
upstream of the debias gate. This is a stronger and more transferable result than the
finding it replaces, and it was produced by the failure rather than despite it.

## A21.5 Process failure

The pre-flight gates adopted after the MSH3/MSH2 chain error include target-identity
verification. That check is written for the **receptor** and passed correctly: 5J89 is
PD-L1. The error was in the **ligand set**, which no gate examined. Target identity must be
verified on both sides of a docking experiment, and the identifier used to fetch actives
must be confirmed to be a validated single-protein target — not merely to exist.

The label "PD-L1" was carried in a hand-written dictionary and never checked against the
ChEMBL record it was paired with. It survived twenty amendments and a public release.

## A21.6 Actions

1. All docking stopped (local, harvest, Mpro queue, relaunch watcher). Caches retained:
   1,145 of 2,526 in the merged cache; shard caches 233 and 363.
2. Zenodo deposit **10.5281/zenodo.21824633 not published** — it was still a draft.
3. `AegisMindApp/retrospective-benchmark` is public and carries the error; correction to
   follow rather than deletion, since the repository is now citable.
4. `CHEMBL4005` is additionally mislabelled: it is **PI3Kα**, not aldose reductase.

---

# Amendment 22 — does an incoherent library pass the gate more easily? (2026-08-07)

**Written and committed before the experiment is run.**

## A22.1 Why

Amendment 21 reports that the one library admitted by the debias gate was not a target,
and proposes a mechanism: a set with no shared pharmacophore has no coherent property
signature, so a property-based classifier cannot separate it from property-matched decoys,
and it passes more easily than any real target.

As it stands that mechanism is supported by **one accident**. A single malformed library
scoring 0.585 is an anecdote, and it is legitimate to object that we made a data-handling
error and are presenting it as a finding. This amendment converts the claim into a
measurement, or refutes it.

## A22.2 Design

The debias gate (`decoys.debias_auroc`) is a logistic regression on seven descriptors. It
requires no docking, so this experiment costs hours rather than weeks.

- **Group A — real targets.** Each target with retained panel actives, subsampled to a
  fixed *n* so that set size cannot confound the comparison.
- **Group B — pooled pseudo-targets.** Sets of the same size *n*, drawn at random from the
  union of actives across ≥ 3 different targets. These are incoherent by construction and
  are exactly what a malformed identifier produces.
- **Group C — positive control.** `CHEMBL612545`, the malformed library discovered in
  Amendment 21, at the same *n*.

Every set passes through the identical decoy protocol and the identical gate.

## A22.3 Pre-committed prediction

**Group B scores lower gate AUROC than Group A** — that is, incoherent libraries pass the
gate more easily than real targets. Group C sits within the Group B distribution rather
than the Group A distribution.

Primary comparison: Mann–Whitney U on Group A versus Group B AUROCs, two-sided, α = 0.05.
Reported alongside the medians and the full distributions.

## A22.4 Falsification, and what it costs us

**If Group B does not score lower than Group A, the mechanism proposed in Amendment 21 is
wrong.** The explanation for CHEMBL612545 passing would then be unknown — plausibly chance —
the paper would have no central finding, and we would not publish it. The repository and
this pre-registration would remain as the record.

We state this before running because the result decides whether a manuscript exists, and a
prediction registered afterwards would be worthless.

## A22.5 What this cannot establish

One decoy generator, one descriptor set, one threshold. A positive result shows that *this*
gate admits incoherent libraries preferentially; it does not establish the behaviour of
every debias procedure. Pooled sets are also an imperfect model of a malformed identifier —
they are drawn from real targets and so retain more internal structure than a genuine
catch-all bucket. If anything that biases **against** the prediction, which we regard as
the conservative direction.

---

# Amendment 23 — the Amendment 21 mechanism is REFUTED (2026-08-07)

**Result of the experiment pre-registered in Amendment 22. Reported as registered,
including that it goes against us.**

## A23.1 Outcome

| Group | n | Median gate AUROC | Range |
|---|---|---|---|
| A — real single targets (n ≈ 30) | 7 | **0.6167** | 0.5084 – 0.8873 |
| B — pooled pseudo-targets (n ≈ 30, ≥ 3 sources) | 12 | **0.6765** | 0.4773 – 0.8172 |
| C — `CHEMBL612545` malformed control | 1 | **0.4832** | — |

Mann–Whitney U = 38.0, two-sided **p = 0.773**. The pooled median is **higher** than the
real median, not lower.

**The prediction registered in A22.3 is refuted.** Incoherent libraries do not pass the
debias gate more easily. The mechanism proposed in Amendment 21 — that a set with no shared
pharmacophore has no property signature to separate — does not survive its own test, and
the claim is withdrawn.

Per A22.4 this means the manuscript has no central finding as framed, and the reframed
paper of 7 Aug is not publishable in that form.

## A23.2 What still stands, and the corrected hypothesis

The *observation* is unaffected: `CHEMBL612545` scores 0.4832 here and 0.585 at full count,
passing comfortably while real targets do not. What is wrong is our explanation.

Pooling actives from three or four real targets does not produce what a malformed
identifier produces. A pooled set is a union of a few internally coherent clusters — still
structured, still distinguishable from a broad drug-like background. The malformed bucket
holds 2.3 million compounds spanning all of ChEMBL, so a sample from it is effectively a
**random draw from drug-like chemical space** — which is the same space the decoys are
drawn from.

The corrected hypothesis is therefore not *incoherence* but **distributional identity**: the
pseudo-target passes because its actives and its decoys are sampled from the same
distribution, so no classifier can separate them. That is a sharper claim and it makes a
different prediction — a set of random drug-like ChEMBL compounds used as pseudo-actives
should score ≈ 0.5 and pass trivially.

It is **not tested here**, and is recorded as a hypothesis, not a result. Testing it
requires its own pre-registration.

## A23.3 An unregistered but robust observation

Group A was run at n ≈ 30 while the Amendment 16/17 re-gate used full active counts. The
same targets, same protocol, differing only in active count:

| Target | AUROC at n ≈ 30 | Verdict | AUROC at full count | Verdict |
|---|---|---|---|---|
| PI3Kα | 0.5084 | **qualifies** | 0.6510 | excluded |
| Mpro | 0.5949 | **qualifies** | 0.6653 | excluded |
| MET | 0.5204 | **qualifies** | 0.7608 | excluded |
| BRD4 | 0.6167 | excluded | 0.6641 | excluded |
| FXa | 0.6936 | excluded | 0.7081 | excluded |
| IRAK4 | 0.8873 | excluded | 0.9026 | excluded |
| VEGFR2 | 0.7856 | excluded | — | — |

**3 of 7 targets qualify at n ≈ 30; 0 of 6 qualify at full active count.** Every target
moves in the same direction — subsampling lowers the measured AUROC and makes qualification
easier — and no target moves the other way.

This is a far stronger version of §3.2, which rested on a single reversal (Mpro) plus one
that turned out to belong to the pseudo-target. It is now seven targets moving consistently,
with three verdict changes. It was not pre-registered as a hypothesis and is reported as an
observation arising from the control group of another experiment.

## A23.4 Consequence

The paper's centre of gravity moves from "a debias gate admits non-targets preferentially"
(refuted) to **"debias gate verdicts are an artefact of active-set size"** (measured, seven
targets, consistent direction). The malformed-target episode becomes a documented failure of
target validation rather than the central result.

Whether that is publishable is a separate judgement, taken with the numbers in hand rather
than in advance of them.

---

# Amendment 24 — does a library with no target relationship pass the gate? (2026-08-07)

**Written and committed before the experiment is run.**

## A24.1 Why

Amendment 23 refuted the "incoherence" mechanism and proposed a replacement: the malformed
library passes not because its actives are incoherent but because they are drawn from the
**same distribution as the decoys**. Decoys are selected from drug-like ChEMBL space; a
sample from a 2.3-million-compound catch-all is effectively a random draw from that same
space, so the two are indistinguishable by construction.

That hypothesis is recorded in A23.2 and untested. This amendment tests it.

## A24.2 Design

- **Group D — null libraries.** Sets of *n* ≈ 30 molecules sampled at random from drug-like
  ChEMBL space, used as "actives". These have **no target relationship of any kind**. Each
  passes through the identical decoy protocol and the identical gate.
- Compared against Group A (real targets, median 0.6167) and Group B (pooled pseudo-targets,
  median 0.6765) already measured under Amendment 22, and Group C (`CHEMBL612545`, 0.4832).

## A24.3 Pre-committed prediction

**Group D scores ≈ 0.5 and passes the 0.60 gate in the large majority of replicates**, and
scores **lower than Group A**. `CHEMBL612545` (Group C, 0.4832) sits with Group D.

Primary comparison: Mann–Whitney U on Group A versus Group D, two-sided, α = 0.05, plus the
proportion of Group D replicates scoring below the 0.60 threshold.

## A24.4 Falsification

**If Group D scores like Group A (≥ 0.60, no significant difference), the distributional
hypothesis is also wrong.** We would then have no explanation for the malformed library's
behaviour, and — per the standard set in A22.4 — the malformed-target line would be dropped
entirely. What would remain of the study is the active-set-size finding of A23.3 and
nothing else.

Two mechanisms will then have been proposed and refuted. We record in advance that a third
will not be attempted: continuing to generate explanations until one survives is the
re-rolling prohibited by §8 of the original registration.

## A24.5 Anticipated objection

A positive result is close to tautological: a property-matching procedure will of course
fail to separate two samples from the same distribution. That is the point rather than a
defect — the gate is *used* as a check that a benchmark is unbiased, and a check passed
most easily by a library with no signal cannot support that use. We record the objection
here so that it is not mistaken for an unanticipated criticism later, and we make no claim
that the observation is novel until the literature has been checked.

---

# Amendment 25 — the null distribution is the finding (2026-08-07)

**Result of the experiment pre-registered in Amendment 24.**

## A25.1 Outcome

| Group | n | Median | Range | Pass (< 0.60) |
|---|---|---|---|---|
| A — real targets (n ≈ 30) | 7 | 0.6167 | 0.508 – 0.887 | **3/7 (43%)** |
| D — null libraries, no target relationship | 12 | 0.6030 | 0.511 – 0.821 | **6/12 (50%)** |
| C — `CHEMBL612545` malformed | 1 | 0.4832 | — | passes |

Mann–Whitney U = 47.0, two-sided **p = 0.711**.

The A24.3 prediction was that Group D would score **≈ 0.5** and pass in the **large
majority** of replicates. It scored 0.603 and passed in half. **The prediction is not
supported on its quantitative content**, and the distributional hypothesis is not
established.

## A25.2 Both proposed mechanisms are now dropped

Amendment 21 proposed incoherence; Amendment 23 refuted it. Amendment 23 proposed
distributional identity; this amendment does not support it. Per the commitment recorded in
A24.4, **no third mechanism is proposed.**

The correct reading requires no mechanism at all. `CHEMBL612545` scored 0.4832, which sits
below but adjacent to a null distribution spanning 0.511–0.821. Given that spread, a single
library scoring 0.48 is **unremarkable**. The malformed library did not pass because of any
special property. It passed because verdicts at this active-set size are close to
uninformative — which is the actual finding, arrived at by two failed explanations.

## A25.3 What the experiment did establish

**A library with no target relationship whatsoever is statistically indistinguishable from
a real target under this gate** (0.603 vs 0.617, p = 0.711), and passes at the same rate
(50% vs 43%).

**Five of seven real targets fall inside the null range.** At n ≈ 30 the metric's null
distribution spans 0.511–0.821 — wider than the entire interval between "clean" and
"badly biased" as the 0.60 threshold is normally interpreted.

This is a statement about the metric, not about any target. A bias metric whose null
distribution is that wide cannot support a per-target verdict at this sample size. It also
explains the instability of A23.3 mechanically: small active sets produce high-variance
estimates, so verdicts flip when *n* changes.

## A25.4 The study's remaining claim

Two measured results survive, and they are the same result seen from two directions:

1. **Verdicts shift systematically with active-set size.** 3 of 7 targets qualify at
   n ≈ 30; 0 of 6 at full active count; every target moves in the same direction (A23.3).
2. **At n ≈ 30 the null distribution spans 0.511–0.821**, and real targets are not
   separable from null libraries (p = 0.711). Five of seven real targets lie inside the
   null range.

Together: **a decoy-bias metric of this form, applied to active sets of the size routinely
used, does not support the per-target admissibility judgements that are made with it.**
Compute-driven active caps are near-universal and rarely reported.

That is narrower than anything this study set out to claim, and it is the first claim in it
that has survived an experiment designed to kill it. Nothing here concerns whether the
AegisMind platform outperforms docking; that question (§1 of the original registration) was
never executable and is not addressed.

## A25.5 Literature position, stated honestly

Training a classifier to separate actives from property-matched decoys, scored by AUROC, is
established practice as a *bias diagnostic* — it is used, for example, in the DeepCoy work
(Imrie et al., *Bioinformatics* 2021) and in Sieg et al. (*JCIM* 2019). Our use of it as a
prospective inclusion criterion, and the observations above about its null distribution and
sample-size dependence, were not located in two targeted literature searches. **Two searches
are not a systematic review**, and no novelty claim should be made until a proper search is
done.

---

# Amendment 26 — publication artefacts withdrawn (2026-08-07)

Recorded so that the actions taken on this study are as fully logged as its results.

## A26.1 What was withdrawn

- **Zenodo deposit `10.5281/zenodo.21824633` deleted.** It was never published: the DOI was
  reserved but never minted, and `https://doi.org/10.5281/zenodo.21824633` has never
  resolved. Nothing was released and there is nothing to retract.
- **`AegisMindApp/retrospective-benchmark` set private.** It had been public for
  approximately one day, with zero stars, forks and watchers. Verified inaccessible
  anonymously via API, raw content and web page.

## A26.2 Why

The manuscript's central claim did not survive its own tests. Amendment 21's mechanism was
refuted by Amendment 23; the replacement was not supported by Amendment 25. What remains is
a narrow, measured claim about the null distribution and sample-size dependence of a
decoy-bias metric (A25.4) — worth a short methods note at most, and not currently being
submitted.

Publishing a repository and minting a permanent DOI in support of a manuscript that is not
being submitted, and whose framing changed twice in the two days the artefacts existed, has
no purpose. The work is retained in full privately.

## A26.3 What this is not

This is **not** a retraction, and it is not the withdrawal of a public claim: no claim was
ever published. Nothing external cited either artefact. Had the DOI been minted, or had the
repository been cited, the correct action would have been correction in place rather than
removal — which is what was done on 7 Aug when the repository was public and carried the
Amendment 21 error.

## A26.4 Condition on any future use

If any part of this study is submitted anywhere, the repository and this pre-registration
must be made public again as a condition of submission. The residual limitation stands and
must be stated: an extract cannot carry the original pre-registration commit `92c6f8bf`, so
the freeze date is attested and reviewable on request rather than independently verifiable
from the artefact.

## A26.5 Record of what was learned

Retained locally and unaffected by the withdrawal: twenty-six amendments; the gate results
for every target assessed; the pooled-target and null-library experiments and their data;
the docking caches; and the process failures — silent rigid-ligand preparation, a target
identifier that was not a target, a stale hardcoded version string, and two proposed
mechanisms that did not survive testing.

---

# Amendment 27 — corrections to the public extract (2026-09-21)

Recorded on republication. Neither item changes any result; both correct the **artefact** rather
than the analysis, and both are the kind of defect this study exists to document.

## A27.1 The published `make_bundle.py` was the pre-fix version

The extract taken on 7 Aug 2026 shipped the bundling script as it stood on that date. It resolved
the receptor by the naming convention `receptors/{target_id}_receptor.pdbqt` rather than from the
path recorded in `receptors.json`. On PD-L1 the convention path held **5J89**, which had failed its
redocking gate, while the configuration pointed at **5J8O**, which passed. A run would have docked
455 ligands against the failed structure using the passing structure's box — a different crystal
frame, so the box sat in empty solvent — and returned a well-formed AUROC with no exception and
nothing in any log. It was caught by comparing receptor md5s (`b83fcd22` against `4ca3a9d9`)
minutes before compute, and fixed in the private tree on 23 Aug.

**That fix had not reached this repository**, so anyone reusing the published script inherited the
defect. The corrected version is now here: it reads the receptor path from `receptors.json` and
fails loudly when the file is absent.

The general rule, which is cheaper than the md5 check that caught it: **a path built by naming
convention will eventually disagree with the path in the configuration.** Prefer the configuration
field.

## A27.2 The configured Mpro receptor was 7K40; the study used 7VU6

Amendments 19 and 20 list **7K40** as the SARS-CoV-2 Mpro receptor, and it appears in a validation
note. Those entries are correct as a record of what was configured at the time and are left
unedited, because this log is append-only. But 7K40 was subsequently **rejected**, for two reasons
recorded in the Mpro result:

- its co-crystal ligand is **covalent to Cys145**, which a non-covalent scoring function cannot
  represent; and
- the configured box sat **6.3 Å off the ligand centroid**.

Four alternatives were tested and **7VU6 was adopted**, passing the redocking gate at
exhaustiveness 32 — the value the screen uses, not the 16 the gate was first validated at — with
the top pose at **1.01 Å at rank 1**, the next pose 5.88 Å and 0.38 kcal/mol behind.

A reader following 7K40 out of this document would be following a receptor the study discarded.

## A27.3 What is still not in this repository

Stated so the omission is not mistaken for absence of the work. The extract carries the
pre-registration, the amendment log and the code. It does **not** carry the per-target result
documents — Mpro, PD-L1, Factor Xa, the target-selection record, or the model-versus-docking
comparison. Those postdate the 7 Aug extract in part (Factor Xa was run 25–29 Aug) and have not
been exported. Until they are, results quoted from this study elsewhere cannot be checked against
this repository, and should be treated as attested rather than verifiable — the same standing
limitation `PROVENANCE.md` states about the freeze commit.

---

# Amendment 28 — result documents exported (2026-09-21)

Amendment 27.3 recorded that this repository carried the design and the code but none of the
per-target results, and that results quoted from this study elsewhere therefore could not be
checked against it. That is now closed for the results themselves.

## A28.1 What was added

`results/` — Mpro, PD-L1 and Factor Xa result documents, the two pre-analyses fixing their reading
rules before the data existed, the target-3 selection record, the model-versus-docking comparison
(which carries its own correction: its panel is `CHEMBL612545`, not a target), and the assessment
of whether modern deep-learning docking had been tested.

## A28.2 Each carries a status banner, and the public copies therefore differ from the private ones

Three headline numbers moved after these documents were written, so publishing them unmarked would
leave a reader with a superseded figure and no signal. Each file gains a banner pointing at
`RESULTS_STATUS.md`; the bodies are otherwise byte-identical to the private originals and are **not**
edited to match the present, because a record edited to agree with today cannot show what was
believed when.

The changes, none of which flipped a verdict:

- **Mpro** 0.427 -> **0.4530**, and **Factor Xa** 0.6657 -> **0.6775**, both from repairing
  receptors that were missing polar hydrogens. Below chance and above chance respectively, before
  and after.
- **PD-L1's conclusion was replaced by a stronger one.** "Not above chance" becomes "the panel is
  property-separable, so no docking result on it bears on binding" — seven free descriptors reach
  0.9145 on it.

## A28.3 The claim the study actually supports

Stated here because no single result document states it: every one of them compares docking against
chance, and chance is the wrong comparator. Against seven free physicochemical descriptors, docking
adds **-0.0023 [-0.0044, -0.0002]** on Mpro and **+0.0150 [-0.0033, +0.0339]** on Factor Xa. On
neither does it add anything demonstrable; on Mpro it significantly subtracts. The general claim
"docking is below chance" is **withdrawn** — true of Mpro, false of Factor Xa.

## A28.4 What remains attested rather than verifiable

The repaired-receptor and descriptor-baseline figures come from a separate analysis not exported
here. They are quoted so a reader is not left with a superseded number, but they cannot be checked
from this repository. The limitation is the same one `PROVENANCE.md` states about the freeze commit,
and it is named rather than left for a reader to discover.

---

# Amendment 29 — A28.4 was wrong: the correcting figures are checkable (2026-09-21)

Amendment 28.4, written earlier the same day, said the repaired-receptor and descriptor-baseline
figures "cannot be checked from this repository" and were "attested, not verifiable". The first
half is true and the second is **false**, and the distinction is the entire point of saying it.

Those figures are public. They are in the companion repository
`github.com/AegisMindApp/screening-decomposition` — `analysis/docking_value/MARGINAL_VALUE.md`
carries the marginal-value table verbatim, `analysis/receptor_prep/` carries the donor defect that
caused the numbers to move, `analysis/docking_audit/RESULT.md` carries PD-L1's 0.9145, and
`analysis/boltz2/` carries the Boltz-2 margins. `RESULTS_STATUS.md` now maps each figure to its
file.

**How the error happened, since it is the same shape as several the study documents.** The question
asked was "are these figures in *this* repository", the answer was no, and that was recorded as
"unverifiable" — a property of the claim — when it was a property of where I had looked. An absence
in one place read as an absence everywhere. That is the failure this project has logged repeatedly
under other names, committed here in a paragraph whose purpose was to be scrupulous about exactly
this.

**What was NOT done, deliberately.** The obvious response is to copy those directories into this
repository. Two copies of one result in two repositories drift, and a reader then cannot tell which
is current. One location, cross-referenced from the other, is correct. A28.4 is superseded by this
amendment rather than edited, per the append-only rule.

The freeze-date limitation in `PROVENANCE.md` is unaffected and still stands.
