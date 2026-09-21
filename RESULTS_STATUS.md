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

## What is still not checkable from this repository

The repaired-receptor and descriptor-baseline numbers quoted above were produced in a separate
analysis that has not been exported here. They are stated so that a reader is not left with a
superseded figure, but within this repository they are **attested, not verifiable** — the same
standing limitation `PROVENANCE.md` records about the pre-registration freeze commit. The result
documents in `results/` are verifiable here; the corrections to them are not yet.

Saying so is cheaper than the alternative, which is a reader discovering it.
