# Target 3: Factor Xa, chosen by a gate PD-L1 would have failed

> **Status — read `../RESULTS_STATUS.md` first.** This document is the record of the analysis on the date it carries. The selection stands; Factor Xa was run and reported.

**25 August 2026.** All screening local — ChEMBL queries and RDKit descriptors, no docking, no
Kaggle. That is the point: PD-L1 cost three Kaggle sessions and two resume cycles to produce a
panel a two-minute local check would have rejected.

## The gate, in order

| step | check | threshold |
|---|---|---|
| 1 | Identity | SINGLE PROTEIN, real organism, ≥1 component |
| 2 | Balance | ≥200 in the minority class (PD-L1 had 84) |
| 3 | **Separability** | **\|MW-alone AUROC − 0.5\| ≤ 0.15** |
| 4 | Electrophile load | low fraction of covalent warheads among actives |
| 5 | Receptor redocking gate at exhaustiveness 32 | *next step, not yet run* |

Step 3 is new, and calibrated on the only two panels we have run: Mpro at 0.620 (0.120 from
chance) **passes** — its size confound ran in docking's favour and docking failed anyway, which
made the result stronger; PD-L1 at 0.1228 (0.377 from chance) **fails**.

## Result

14 candidates screened, 9 reached stage 2, 7 passed the MW gate.

| target | MW AUROC | dev | 7-desc | actives | inactives | electrophilic |
|---|---|---|---|---|---|---|
| EGFR | 0.5338 | 0.034 | 0.646 | 345 | 987 | **23.4%** |
| VEGFR2 | 0.4655 | 0.035 | 0.766 | 428 | 581 | **17.3%** |
| SRC | 0.5472 | 0.047 | 0.663 | 416 | 443 | 13.5% |
| AKT1 | 0.5542 | 0.054 | 0.648 | 397 | 258 | 1.2% |
| COX-2 | 0.5727 | 0.073 | **0.844** | 423 | 726 | 0.0% |
| **Factor Xa** | **0.6011** | **0.101** | **0.652** | **467** | **539** | **0.0%** |
| JAK2 | 0.3977 | 0.102 | 0.660 | 408 | 325 | 0.9% |

**Rejected by the MW gate:** HDAC2 (0.6681, dev 0.168) and **mTOR (0.8166, dev 0.317)** — mTOR is
almost as badly confounded as PD-L1 was, in the opposite direction. The gate is not decorative; it
removed two of nine.

**Rejected earlier:** MET, on 147 inactives.
**Unmeasured, not rejected:** PI3Kα, ABL1, HDAC1 — persistent ChEMBL HTTP errors. Recorded as
unmeasured. AChE returned 4,518/1,511 on one attempt and errored on later ones, so it is a live
candidate whose stage-2 numbers were never obtained.

## Why Factor Xa over the higher-ranked ones

EGFR and VEGFR2 top the MW ranking but carry **23.4%** and **17.3%** acrylamide/Michael-acceptor
actives. Docking scores a non-covalent pose; a covalent binder's potency reflects bond formation
the scoring function does not model, so those compounds are unrankable by construction. That is
the CTSS defect — a non-covalent screen against a covalent-inhibitor site — arriving through the
ligand set instead of the receptor. Nearly a quarter of EGFR's actives is too much.

COX-2 is clean on warheads but its seven-descriptor baseline is **0.844**: properties alone
classify that panel very well, leaving docking little room to be informative and making any
apparent performance hard to attribute.

**Factor Xa is the only candidate clean on every axis at once:** 0.0% electrophilic, a
seven-descriptor baseline of 0.652, and the best class balance in the set — 467 actives against
539 inactives, a 46% active fraction, where PD-L1 was 80% active on 84 inactives.

One further property decided it. Its MW AUROC of **0.6011 runs in the same direction as Mpro's
0.620**: higher molecular weight predicts *active*, and Vina prefers larger molecules, so the size
confound runs **in docking's favour**. That is the asymmetry that made the Mpro result strong — if
docking fails while the confound is helping it, the failure is real. It also makes Factor Xa
directly comparable to the one result we still trust, rather than a third incommensurable panel.

## Corrections made during the screen

**The covalent check over-flagged and would have rejected a good candidate.** The first version
summed all warhead classes and reported SRC at 46.7% "covalent" — almost entirely plain aryl
nitriles, which are ordinary non-covalent substituents. Its docstring claimed the patterns
"undercount"; for nitriles the opposite is true. Split into ELECTROPHILIC (acrylamides, vinyl
sulfones, chloroacetamides — decision-relevant) and AMBIGUOUS (nitriles, aldehydes, boronic acids
— reported, not acted on). SRC's real electrophilic load is 13.5%.

**Rate-limit errors were nearly mistaken for absent candidates.** 8 of 14 targets failed stage 1
on ChEMBL HTTP errors, and acetylcholinesterase had *succeeded* on an earlier attempt with the
same query. Retried with backoff, SRC and JAK2 both passed — two viable candidates that a
single-pass screen would have discarded silently, with a log that looked complete. Failures are
labelled "count failed", never "reject".

## Next

Receptor selection and the redocking gate at exhaustiveness 32, reusing `gate_pdl1_receptors.py`.
Factor Xa has abundant non-covalent co-crystals, so the covalency trap that made 3N4C ungateable
should not arise — but the gate decides that, not this note. **No compute is committed until a
receptor passes.**
