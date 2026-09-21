# PD-L1 full set: not above chance, and the coverage gap is class-selected

> **Status — read `../RESULTS_STATUS.md` first.** This document is the record of the analysis on the date it carries. Superseded in the stronger direction: a later audit found seven free descriptors reach **0.9145** on a PD-L1 panel, so it is property-separable and **no docking conclusion on PD-L1 bears on binding**. This document's own caveats already point that way.

**3 September 2026.** Final 98 compounds completed, merged with the existing protonated cache.

    n = 417   actives 334 (80.1%)
    AUROC 0.5661   95% CI [0.4876, 0.6419]   -> CI includes 0.5, NOT above chance

## This resolves the gate flip

The benchmark record noted PD-L1 flipping gate admissibility under active subsampling —
failed at a 30-active cap, qualified at 82. With the full set assembled, **PD-L1 does not
qualify**: the interval includes chance. The earlier "qualified" reading was an artefact of
subsampling, not a property of the target.

## Second finding: the missing 38 are not random

    scored 417 / 455    missing 38
    active fraction   missing 65.8%   scored 80.1%   difference -14.3%

The coverage gap is **class-selected** — the same failure shape as arm 2's 149 size-selected
dropouts. Only one compound is a recorded Vina failure, so these are largely compounds never
dispatched rather than compounds that failed to dock. The AUROC above is therefore untrustworthy
independently of its interval, and should not be quoted as a clean negative.

A third issue worth recording: the single recorded failure, `bench_CHEMBL6048352`, logs
`rc=0 no-score: 5 -10.1 1.505 2.355` — Vina **did** produce a score (−10.1) and the parser
missed it. That is a silent-failure pattern, not a docking failure.

## Standing caveat

At 80.1% actives this set is not a realistic screening scenario; AUROC is unstable at that
imbalance regardless of the above.
