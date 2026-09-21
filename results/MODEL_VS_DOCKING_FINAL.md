# Phase 2 (PDBbind GNN) vs Vina docking — CHEMBL612545 panel

> **Status — read `../RESULTS_STATUS.md` first.** This document is the record of the analysis on the date it carries. Carries its own correction banner: the panel is `CHEMBL612545`, which is not a target. Retained because removing it would hide the error.

> **CORRECTION, 17 Aug 2026 — the target is not PD-L1, and this is not evidence about docking.**
>
> This panel's actives were fetched with ChEMBL id `CHEMBL612545`, verified against the ChEMBL API
> as `pref_name: "Unchecked"`, `target_type: UNCHECKED`, no organism, **0 target components** — a
> catch-all bucket of ~2.3M bioactivity records with no validated target assignment. Real human
> PD-L1 is `CHEMBL3580522`. The receptor docked against (PDB 5J89) IS PD-L1; the ligand set is not.
> Sampled "actives" are HCN1 channel blockers and IL-6 release inhibitors. The benchmark paper
> documents this as its own §3.3 finding; I ran this comparison without reading it.
>
> **What must be withdrawn:** "docking underperforms a property-only baseline". There were no true
> actives to enrich, so the comparison says nothing about docking's ability to rank real binders.
>
> **What survives, reframed:** this is an unintended *negative control* — a real receptor paired
> with a ligand set carrying no true signal. Read that way it is informative, and not in docking's
> disfavour:
>
> | predictor | AUROC | reading |
> |---|---|---|
> | chance | 0.500 | — |
> | Vina docking | **0.533** | at chance, which is the CORRECT behaviour on a null set |
> | Phase 2 PDBbind GNN | **0.677** | **above chance on a set with no true actives** |
>
> Vina behaves properly here. The learned model finds 0.677 of apparent enrichment where no
> enrichment exists, which is evidence it reads property//scaffold structure rather than binding —
> the same concern Limitation 9 raises about PDBbind provenance, now with a measurement behind it.
>
> The comparison still needs running on a genuine target before anything can be claimed about
> either predictor. SARS-CoV-2 Mpro (`CHEMBL4523582`, verified SINGLE PROTEIN, 85 actives) is the
> obvious candidate.


Run 16 Aug 2026 at **2,495/2,526 docked (98.8%)**. 31 compounds never scored and are not
recoverable by retrying; they are the high-torsion tail, so their absence is not neutral, but at
1.2% it cannot move the comparison.

Paired comparison on the 2,490 compounds both predictors scored (the only fair one):

| predictor | n | actives | AUROC | EF@1% | EF@5% |
|---|---|---|---|---|---|
| chance | — | — | 0.500 | 1.00 | 1.00 |
| **decoy property-only baseline** | — | — | **0.585** | — | — |
| Vina docking (exh 32) | 2490 | 77 | **0.533** | **0.00** | 1.30 |
| Phase 2 (PDBbind-trained GNN) | 2490 | 77 | **0.677** | **5.17** | 3.13 |

## Two findings

**1. Docking underperforms a property baseline.** Vina's AUROC of 0.533 is barely above chance
and sits *below* the 0.585 achievable from molecular descriptors alone. On this target, seven
descriptors beat the physics. EF@1% is 0.00 — in the top 25 compounds Vina finds no actives.

This is the same phenomenon the AcrB and CTSS decoy controls isolated directly: docking score
tracks molecular size (r² = 0.70 against MW on the CTSS panel), and the same four large
molecules top unrelated targets. If score is largely a size read-out, a size baseline should beat
it, and here it does.

**2. The learned model beats both.** AUROC 0.677 clears chance and clears the property baseline,
and EF@1% = 5.17 against Vina's 0.00. At the 1% cutoff (~25 compounds) it recovers roughly four
actives where docking recovers none. Small absolute numbers — treat EF@1% as directional.

## The limitation that governs how far this can be stated

Phase 2 is trained on PDBbind. **PD-L1 complexes may be in that training set**, and we have not
excluded them. The model's advantage over docking is therefore not cleanly attributable to
generalisation; some part may be memorisation of this target's chemistry. This is Limitation 9 in
the benchmark paper and it is not a formality — it is the difference between

- defensible: *on this panel, the PDBbind-trained model ranks actives better than Vina, with a
  training-set overlap we have not ruled out*; and
- not defensible: *our model outperforms docking*.

Establishing the second needs a target absent from PDBbind, or a checkpoint retrained with PD-L1
complexes held out.

## Status of the claim

This is a measurement of two predictors against measured labels. It is **not evidence about any
compound**, and nothing here licenses a wet-lab approach on any molecule in the panel.


## Training-set leakage: measured, not assumed (18 Aug 2026)

Earlier notes here treated "the PDBbind-trained model's margin may be memorisation" as an
unresolvable caveat requiring a held-out retrain. It is resolvable in minutes, and the answer is
that there is no target leakage to worry about.

Phase 2 fine-tunes on the PDBbind v2020 **refined** set — 5,316 complexes (`pdbbind_data.py` loads
`INDEX_refined_data.2020`; the module docstring's "general set (~19K)" is aspirational and wrong).
Resolving each evaluation target's UniProt to its PDB entries via the RCSB search API and
intersecting with those 5,316:

| target | UniProt | PDB entries | in PDBbind refined set |
|---|---|---|---|
| SARS-CoV-2 Mpro | P0DTD1 | 3,653 | **0** |
| human PD-L1 | Q9NZQ7 | 78 | **0** |
| HIV-1 protease (control) | P03366 | 452 | 74 |

The HIV-1 protease row is there to show the method detects leakage when it exists — 74 of its
complexes ARE in the set, so a zero is a real zero and not a broken query.

**Consequence:** the model's performance on Mpro cannot be memorisation of Mpro complexes, because
it has never seen one. The same holds for the CHEMBL612545/PD-L1 receptor. The claim "our model
outperforms docking" is still not established — that needs the Mpro comparison to finish — but the
specific objection that it is reading training data is answered.

**Residual, weaker channel: homology.** Two SARS-CoV (2003) replicase structures (2z94, 3r24) are
in the refined set, and SARS-CoV 3CL protease is a close homolog of SARS-CoV-2 Mpro. MERS and
HCoV-229E contribute none. So the model has seen two related-protease complexes out of 5,316. That
is a far weaker channel than target leakage and should be stated rather than ignored, but it does
not support a memorisation explanation.

A held-out retrain is therefore NOT needed for this comparison. The 658 MB refined set and the
index remain in gs://aegismind-tpu-results/phase2_data/ if a homology-excluded retrain is ever
wanted.
