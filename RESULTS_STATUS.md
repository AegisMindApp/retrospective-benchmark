# Where each result now stands

The documents in `results/` are records of an analysis on the date each one carries. Several
headline numbers have moved since, and one target's conclusion was replaced by a stronger one.
This file states the current position for each, so that a reader arriving at any single document
by a direct link knows whether it is the latest word.

Nothing here is a retraction of a result document. Each is left as written, because a record that
is edited to match the present cannot show what was believed when.

## The current position, in one table

| target | as recorded | now | changed by |
|---|---|---|---|
| SARS-CoV-2 Mpro | AUROC **0.427** [0.381, 0.471], 753 compounds | **0.4530**, 751 compounds | receptor repair |
| PD-L1 | AUROC **0.5661** [0.4876, 0.6419], not above chance | **conclusion replaced** — the panel is property-separable, so no docking result on it bears on binding | descriptor audit |
| Factor Xa | AUROC **0.6657** [0.6298, 0.7016], 888 compounds | **0.6775**, 854 compounds | receptor repair |
| `CHEMBL612545` | withdrawn | withdrawn — not a target | Amendment 21 |

## The three changes, and why each happened

**1. Receptor repair moved every docking number.** The receptors used for the original runs were
missing their polar hydrogens. Repairing them changed Mpro from 0.427 to 0.4530 and Factor Xa from
0.6657 to 0.6775, on slightly different compound counts. **No verdict flipped.** Mpro is below
chance before and after; Factor Xa is above chance before and after. The documents' reasoning is
unaffected; only the values moved.

**2. PD-L1's conclusion was replaced by a stronger one.** `PDL1_RESULT.md` reports "not above
chance" and then, to its credit, says the number should not be quoted as a clean negative, because
the 38 missing compounds are class-selected. A later audit settled it in the other direction:
**seven free physicochemical descriptors reach AUROC 0.9145 on a PD-L1 panel.** The panel is very
nearly separable from ligand properties alone, so nothing docked against it is evidence about
binding at any sample size. PD-L1 is not a weak result — it is not a result.

**3. The headline claim of the whole study is narrower than any single document states.** Each
result document compares docking against chance. That is the wrong comparator. Against **seven
free descriptors** — molecular weight, clogP, donors, acceptors, rotatable bonds, TPSA, formal
charge, costing microseconds and needing no structure — docking adds:

| target | Vina | seven descriptors | **docking adds** | 95% CI |
|---|---|---|---|---|
| Mpro (751, 257 act) | 0.4530 | 0.7653 | **−0.0023** | [−0.0044, −0.0002] **sig.** |
| Factor Xa (854, 388 act) | 0.6775 | 0.7129 | **+0.0150** | [−0.0033, +0.0339] **n.s.** |

**On neither target does docking add anything demonstrable over a baseline that is free**, and on
Mpro it significantly subtracts. This is the claim the study supports. "Docking is below chance",
which an earlier framing of this work asserted generally, is **withdrawn** — it was true of Mpro
and false of Factor Xa.

What does clear a pre-registered bar is Boltz-2's binary affinity probability added to those
descriptors: **+0.0934 [+0.0616, +0.1253]** on Mpro and **+0.0751 [+0.051, +0.099]** on Factor Xa,
against a +0.04 margin fixed in advance. Boltz-2 alone does not clear it. So the structure-based
information that helps is not coming from docking.

## Where the correcting figures can be checked

They are not in this repository, and they do not need to be copied here: they live in the
companion repository **`github.com/AegisMindApp/screening-decomposition`**, which is public. Each
figure quoted above resolves to a file there.

| figure | what it is | file in `screening-decomposition` |
|---|---|---|
| **0.4530** | Mpro Vina, repaired receptor | `analysis/pose_ensemble/PROTONATED_RESULT.md`, `analysis/boltz2/RESULT.md` |
| **0.6775** | Factor Xa Vina, repaired receptor | `analysis/receptor_prep/FXA_RESULT.md`, `analysis/docking_audit/RESULT.md` |
| **0.7653 / 0.7129** | seven-descriptor baselines | `analysis/docking_value/MARGINAL_VALUE.md` |
| **−0.0023 / +0.0150** | what docking adds over them | `analysis/docking_value/MARGINAL_VALUE.md` |
| **0.9145** | PD-L1 descriptor separability | `analysis/docking_audit/RESULT.md` |
| **+0.0934 / +0.0751** | descriptors + Boltz-2 margins | `analysis/boltz2/RESULT.md`, `analysis/boltz2/FXA_RESULT.md` |
| the donor defect itself | why repair changed the numbers | `analysis/receptor_prep/DONOR_DEFECT.md` |

**Copying them here was considered and rejected.** Two copies of the same result in two
repositories drift, and a reader then has no way to tell which is current — that is the defect this
project keeps finding in other people's work and has now found four times in its own. One location,
cited from the other, is the correct arrangement.

## What is genuinely still not checkable

The **freeze date** of the original pre-registration. `PROVENANCE.md` states it: this repository is
an extract, its commit SHAs are new, and the original commit `92c6f8bf` (2026-07-14) is in a
private tree. That limitation stands and is not closed by anything above.
