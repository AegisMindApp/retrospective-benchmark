# Factor Xa: docking is above chance — the Mpro inversion does not replicate

> **Status — read `../RESULTS_STATUS.md` first.** This document is the record of the analysis on the date it carries. Its 0.6657 is on the 888-compound panel; on the repaired-receptor panel of 854 the figure is **0.6775**. Both are above chance and the verdict is unchanged.

**Run 25–29 August 2026, Kaggle CPU, AutoDock Vina 1.2.5, exhaustiveness 32, receptor 2W26.**
Analysis fixed in [`FXA_PREANALYSIS.md`](FXA_PREANALYSIS.md) before the bundle was built and not
revised. Raw verdict in `fxa_result.json`.

## Result

**AUROC 0.6657, 95% CI [0.6298, 0.7016]** on 888 compounds (405 active / 483 inactive).

The pre-registered rule keys on the CI alone, three exclusive cases, no `else`:

> **CI entirely above 0.5 → docking is above chance here. Consistent with PD-L1; the Mpro
> inversion stays a single-target result.**

That is the branch taken. This is the third target in the programme and the first where docking
demonstrably discriminates on a panel that cleared every pre-flight gate before compute.

| target | AUROC | 95% CI | reading |
|---|---|---|---|
| Mpro | 0.427 | below 0.5 | ranks non-inhibitors above inhibitors |
| PD-L1 | 0.5948 | above 0.5 | above chance, but on an 80%-active panel with MW-AUROC 0.1228 |
| **Factor Xa** | **0.6657** | **[0.630, 0.702]** | **above chance, on a panel that passed every gate** |

## The five pre-registered checks

1. **Two independent AUROC computations.** Rank-sum 0.6657, direct pairwise count 0.6657 —
   agree to 4 dp. The script aborts if they do not; an orientation slip inverted this
   measurement once on Mpro (0.573 vs 0.427).
2. **Orientation.** Vina affinity is kcal/mol, more negative = stronger, so the score is the
   negated affinity and AUROC > 0.5 means docking works.
3. **MW confound on the panel as actually docked: 0.6037**, against 0.6011 at screening — a shift
   of 0.0026. The confound is unchanged by which compounds happened to dock, so the reading is
   not qualified.
4. **Within-MW-quartile AUROC**, reported only where the minority class holds ≥ 20. All four
   qualified, and they are not uniform:

   | quartile | MW range | n | actives/inactives | AUROC | 95% CI |
   |---|---|---|---|---|---|
   | Q1 | < 438 Da | 221 | 67 / 154 | 0.5152 | [0.432, 0.598] — **spans 0.5** |
   | Q2 | 438–500 Da | 223 | 107 / 116 | **0.7205** | [0.653, 0.788] |
   | Q3 | 500–547 Da | 222 | 108 / 114 | **0.7537** | [0.690, 0.818] |
   | Q4 | > 547 Da | 222 | 123 / 99 | 0.5646 | [0.489, 0.640] — **spans 0.5** |

5. **Analysable n, with failures separated.** 1,006 frozen → 112 ligand-prep failures → 894
   prepared → 6 dock failures → **888 analysable**, 99.3% of prepared.

## The quartile structure is the most informative part

Docking's signal is **concentrated in the middle of the mass range and absent at both ends**.
Within Q2 and Q3 the molecular-weight confound is held nearly constant, so AUROCs of 0.72 and
0.75 there are discrimination docking earns beyond size — the strongest evidence in this
programme that it can. In Q1 and Q4 the CIs span 0.5 and it is at chance.

So "docking works on Factor Xa" is true in aggregate and misleading as a summary. It works on
438–547 Da compounds and does not, measurably, outside that band.

## The six dock failures, which are a protocol limit rather than noise

Five are 2-hour timeouts at exhaustiveness 32 on the most flexible ligands in the panel
(28–86 rotatable bonds; the median docked compound has 8). One, CHEMBL2373367, contains a
**boron atom**, for which AutoDock Vina has no parameters — permanent under any budget. Both
categories are recorded with reasons rather than silently dropped.

The 112 prep failures are a larger exclusion and predate this run.

## What this does not license

**A general claim about docking.** Three targets: Mpro below chance, PD-L1 above chance on a
compromised panel, Factor Xa above chance on a clean one. Three measurements, one compromised,
are not a characterisation of a method. In particular this does **not** retract the Mpro result —
it establishes that Mpro's inversion is not universal, which is what a second and third target
were run to find out.

## Reproducing

```
python3 fxa_analysis.py     # needs rdkit; reads the shard caches
```
