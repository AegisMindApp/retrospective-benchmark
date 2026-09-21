# Mpro: docking ranks measured inactives above measured actives

> **Status — read `../RESULTS_STATUS.md` first.** This document is the record of the analysis on the date it carries. Its AUROC of 0.427 was measured on a receptor later found to be missing its polar hydrogens; repaired, the figure is **0.4530** on 751 compounds. The direction is unchanged and the reasoning here stands.

**SARS-CoV-2 main protease, 753 compounds, 97.7% coverage. Completed 19 Aug 2026.**

| | |
|---|---|
| **AUROC** | **0.427**, 95% CI **[0.381, 0.471]** |
| chance | 0.500 — **outside the interval** |
| Mann–Whitney | p = 1.01 × 10⁻³ |
| actives (n = 257) | mean **−7.392** kcal/mol |
| measured inactives (n = 496) | mean **−7.582** kcal/mol |

Docking scores the compounds that *do not* inhibit Mpro as the stronger binders, and the effect is
statistically significant rather than a wash.

## Why this one counts when the previous attempt did not

The withdrawn CHEMBL612545 result failed because nothing about it had been checked. Every input
here was verified before the run, and two of those checks changed the experiment:

| check | outcome |
|---|---|
| Target identity | 10/10 sampled actives confirmed against CHEMBL4523582 (`SINGLE PROTEIN`, SARS-CoV-2). The pseudo-target scored 1/12 on the same test. |
| Panel construction | Mpro **failed** the study's own decoy debias gate (`excluded_property_separable`), as does every real target. Switched to **measured inactives** (pchembl < 5) — no synthetic decoys, so no debias gate applies. |
| Receptor | Configured 7K40 **rejected**: its ligand is covalent to Cys145, and the box sat 6.3 Å off the ligand centroid. Four alternatives tested; **7VU6 adopted**, redocking gate PASS. |
| Protocol match | Gate re-run at exhaustiveness **32**, the value the screen uses, not the 16 it was first validated at. Top pose 1.01 Å at rank 1, next pose 5.88 Å and 0.38 kcal/mol behind. |
| Bundle integrity | Receptor md5 in the shipped bundle byte-identical to the gated file. |

So the pose prediction on this receptor is demonstrably correct, and the affinity ranking is
demonstrably inverted. Those are separable failures, and only the second one is present here.

## Stability

The estimate did not move as coverage grew and the class skew resolved:

| coverage | actives in subset | AUROC |
|---|---|---|
| 46% | 70.9% | 0.416 |
| 67% | 49.4% | 0.417 |
| **97.7%** | **34.1%** (panel: 35.1%) | **0.427** |

## An error worth recording

The first computation returned **0.573 and "not distinguishable from chance"** — the opposite
conclusion. Cause: ranking ascending on raw Vina score, where a *more negative* value is the
stronger prediction, computes P(active scored weaker), the complement. The published figure is now
cross-checked by two independent methods (rank formula and direct pairwise count), which agree to
four decimal places. Any future rerun should keep that cross-check: the failure is silent, and it
inverts the finding.

## What this supports, and what it does not

**Supports:** on SARS-CoV-2 Mpro, against a receptor whose redocking control passes at the exact
screening protocol, Vina's affinity ranking is significantly worse than chance at separating
measured actives from measured inactives.

**Does not support:** a general claim about docking. This is one target. Measured inactives are
compounds someone thought worth testing, so they resemble the actives more than a screening library
would — a harder contrast than a prospective screen, which is the conservative direction but still
not the same experiment.

**Not reported as enrichment.** EF@1% is 2.09 over 7 compounds at a 34% active rate. That is not
screening enrichment and should not be quoted as such; AUROC is the metric here.

## Corroboration

Three independent lines now point the same way: the AcrB and CTSS decoy controls (score 70%
explained by molecular weight, same four large molecules topping unrelated proteins), the redocking
failures (1STP top pose 6 Å with the near-native pose ranked 6th), and this inverted ranking on a
gate-passing receptor.

## Three-way comparison (752 compounds scored by all predictors)

| predictor | AUROC | 95% CI | p vs chance |
|---|---|---|---|
| **Seven descriptors** (5-fold CV logistic) | **0.763** | [0.723, 0.800] | — |
| chance | 0.500 | — | — |
| Vina docking | **0.427** | [0.382, 0.472] | 1.0 × 10⁻³ |
| Phase 2 PDBbind GNN | **0.417** | [0.372, 0.461] | 1.8 × 10⁻⁴ |

Both structure-based predictors sit significantly below chance. They are **not distinguishable
from each other**: ΔAUROC = +0.010, 95% CI [−0.036, +0.054]. Spearman between them is +0.474, so
they are making correlated errors rather than independent ones.

This also disposes of the only positive claim from the withdrawn CHEMBL612545 analysis. Phase 2
scored 0.677 there and looked like it beat docking; on a real target it scores 0.417 and is
indistinguishable from it. The 0.677 was enrichment of a set with no true actives.

## The size confounder runs in docking's FAVOUR, and it still fails

Our AcrB and CTSS controls showed docking score tracks molecular size, so the obvious objection is
that the inversion is a size artefact. It is not — and checking makes the result stronger:

| | |
|---|---|
| actives mean MW | **435.8** Da |
| inactives mean MW | 400.3 Da (Mann–Whitney p = 6.3 × 10⁻⁸) |
| Spearman(MW, Vina prediction) | **+0.503** |
| **AUROC of molecular weight alone** | **0.620** |

On this panel the actives are the *larger* molecules and docking prefers larger molecules, so the
size bias should have carried it to roughly 0.62. A predictor that read nothing but molecular
weight would have beaten Vina by ~0.19 AUROC. Docking's structure-specific contribution is
therefore not merely uninformative; it is negative enough to destroy a signal that MW alone
captures.

Within-band check, holding size roughly constant (MW quartiles): Vina AUROC 0.414, 0.413, 0.240,
0.451 — inverted in every quartile.

## Verification performed

1. Two independent AUROC computations (rank formula, direct pairwise count) agree to 4 dp for both predictors.
2. Predictions are non-degenerate: Vina 656/752 unique values, Phase 2 752/752.
3. Phase 2 outputs a plausible pKd range (2.60–9.52, median 6.51) for a model with val_rmse 1.449.
4. Checkpoint loaded with 0 missing / 0 unexpected keys.
5. Score orientation verified on both arms — Vina negated, Phase 2 pKd raw. An orientation error had already produced 0.573 and the opposite conclusion once.
6. Bootstrap CIs (5,000 resamples) and Mann–Whitney tests.
7. Paired bootstrap on the AUROC difference.
8. Size confounder tested and excluded.
9. Property baseline computed from this panel by cross-validated logistic regression, not inherited — the 0.585 in the codebase belongs to the CHEMBL612545 decoy generator and does not transfer.

## Outstanding

Training-set leakage is excluded: Mpro has **0** complexes in the PDBbind refined set (HIV-1
protease, with 74, confirms the test detects leakage when present).
