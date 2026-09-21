# PD-L1: what this run can and cannot conclude, fixed before the result

> **Status — read `../RESULTS_STATUS.md` first.** This document is the record of the analysis on the date it carries. The reading rule fixed here was honoured; see `PDL1_RESULT.md` and the note above it.

**Written 23 August 2026, while both shards are still docking.** Nothing here may be revised after
the AUROC is known. The point of writing it now is that "underpowered" is an argument anyone can
reach for once they dislike a number, and it is only worth anything stated in advance.

## The question

Does the inverted docking ranking measured on SARS-CoV-2 Mpro — AUROC 0.427, 95% CI
[0.381, 0.471], significantly below chance — generalise to a second verified target?

This matters beyond the paper. `solver.press/about` currently states that docking ranks
non-inhibitors above inhibitors, and rests that claim on Mpro alone.

## What the panel can detect

461 compounds fetched, 455 prepared, so the analysable panel is roughly **357 measured actives and
98 measured inactives** (pchembl ≥ 6 vs < 5). PD-L1 is inactive-poor in ChEMBL: 187 activity
records below pchembl 5, ~99 unique molecules, and relaxing the cutoff to 5.5 adds only ~40 more
while blurring the class boundary. 98 negatives is the ceiling, not a choice.

Hanley–McNeil standard errors at that panel size:

| true AUROC | SE | 95% CI | what we could say |
|---|---|---|---|
| 0.40 | 0.033 | [0.335, 0.465] | inversion detected |
| **0.427** (Mpro's value) | 0.033 | [0.361, 0.493] | **inversion detected — no margin to spare** |
| 0.44 | 0.033 | [0.375, 0.505] | nothing |
| 0.50 | 0.033 | [0.435, 0.565] | nothing |
| 0.55 | 0.032 | [0.487, 0.613] | nothing |
| 0.60 | 0.031 | [0.539, 0.661] | above chance |

**The detection window is `≤ 0.43` or `≥ 0.60`. Everything between is a blind zone**, and that
blind zone covers most of the interesting range — including every outcome in which docking is
mildly bad, useless, or mildly useful.

## The decision rule

Stated before unblinding, and binding on how the result is written up.

1. **AUROC ≤ ~0.43 with the CI excluding 0.5** → the inversion replicates. Two verified targets,
   two gate-passing receptors. The site copy stands and gets stronger.

2. **AUROC in ~0.44–0.59, CI spanning 0.5** → **report as a power limit, not as a failure to
   generalise.** This panel cannot distinguish "docking works fine on PD-L1" from "docking is
   inverted here too, by less than Mpro's margin". The correct sentence is *the PD-L1 panel was
   too small to test the question*, and the correct action is to say so and stop — not to soften
   the Mpro claim, and not to present it as replication either.

3. **AUROC ≥ ~0.60 with the CI excluding 0.5** → docking beats chance on PD-L1. The inversion does
   **not** generalise, and the site copy must be narrowed to Mpro that day.

Note the asymmetry, and that it is not special pleading: outcome 3 is a real finding and outcome 2
is not, because the panel has power against a large positive effect and no power against a small
negative one. The rule follows from the table above, which was computed before the data existed.

## Why this is written down

The Mpro result was itself reversed once by a silent orientation bug — ranking ascending on raw
Vina score gave 0.573 and "not distinguishable from chance", the exact opposite conclusion. A
number in the 0.44–0.59 band will look like a clean negative result and will be indistinguishable
from an underpowered one. That is this project's most frequent failure mode
(`feedback_silent_failure_looks_like_data`), and the defence is a rule fixed in advance.

## Also required before reading any number

- Merge both shards and confirm coverage; report the analysable n, not the fetched n.
- Compute AUROC by **two independent methods** (rank formula and direct pairwise count) and
  require agreement to 4 dp. This is the check that catches the orientation inversion.
- Confirm score orientation explicitly: Vina is negated, more negative is a stronger prediction.
- Report the molecular-weight-alone AUROC on this panel. On Mpro the size confounder ran in
  docking's favour and it still failed; if that does not hold here, say so.
