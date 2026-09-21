# Factor Xa: what this run can conclude, fixed before any compute

> **Status — read `../RESULTS_STATUS.md` first.** This document is the record of the analysis on the date it carries. The reading rule fixed here was honoured and not revised after the fact.

**25 August 2026, written before the bundle was built.** Nothing here may be revised once a
number exists. Target 3 of the docking-evaluation programme.

## Why this target

It is the first to clear every pre-flight gate before compute was committed, which is the process
change PD-L1 forced:

| gate | Factor Xa | PD-L1, for contrast |
|---|---|---|
| Identity | `CHEMBL244`, SINGLE PROTEIN, *Homo sapiens*, 1 component | passed |
| Balance | 467 actives / 539 inactives — **46.4% active** | 80% active on 84 inactives |
| **MW-alone AUROC** | **0.6011** (0.101 from chance) | **0.1228** (0.377 from chance) |
| Electrophilic actives | **0.0%** | not checked at the time |
| 7-descriptor baseline | 0.652 | — |
| Receptor redocking gate @ exh 32 | **2W26, 0.88 Å at rank 1** | both 5J89 pockets failed, 4.13 / 4.02 Å |

The MW confound runs in the **same direction as Mpro's** (0.6011 vs 0.620: higher weight predicts
active, and Vina prefers heavier molecules). That is deliberate. It means the size confound helps
docking here, so a failure would be meaningful in exactly the way Mpro's was — and it makes this
panel directly comparable to the one result still standing.

## The decision rule, and the flaw it fixes

**PD-L1's rule had a gap and I nearly fell into it.** It keyed branches to point-estimate
thresholds (≤0.43, ≥0.60) *and* CI conditions, so a result at 0.5948 with a CI excluding 0.5 fell
through to an `else` labelled "power limit" — a label its own definition ("CI spanning 0.5")
contradicted. The thresholds were only ever proxies for what the CI does.

So this rule keys on **the CI alone**. Three cases, mutually exclusive and exhaustive, with no
`else` to fall into:

| the 95% CI | conclusion |
|---|---|
| **entirely below 0.5** | The Mpro inversion **replicates**. Docking ranks non-inhibitors above inhibitors on a second gate-passing target. |
| **entirely above 0.5** | Docking is **above chance** here. Consistent with PD-L1; the inversion stays a single-target result. |
| **spans 0.5** | **Not distinguishable from chance.** Report the effect size the panel could have detected. Not evidence for either side. |

No point-estimate threshold appears anywhere in it.

## What the panel can detect

Hanley–McNeil at 467/539:

| true AUROC | SE | 95% CI | reading |
|---|---|---|---|
| 0.43 (Mpro's value) | 0.0180 | [0.395, 0.465] | inversion detected |
| 0.46 | 0.0182 | [0.424, 0.496] | inversion detected |
| 0.47 | 0.0182 | [0.434, 0.506] | spans 0.5 |
| 0.53 | 0.0182 | [0.494, 0.566] | spans 0.5 |
| 0.55 | 0.0182 | [0.514, 0.586] | above chance |

**The blind zone is 0.47–0.53** — about a tenth of the range. PD-L1's was 0.44–0.59, roughly
three times wider. At ~92% dock coverage (PD-L1 achieved 91.9%) this degrades only to
[0.393, 0.467] at 0.43, so the conclusion does not depend on full coverage.

## Checks that must run before the number is read

1. **Two independent AUROC computations** — rank formula and direct pairwise count — agreeing to
   4 dp. An orientation slip already inverted this measurement once on Mpro (0.573 vs 0.427).
2. **Orientation stated explicitly**: Vina affinity negated, more negative = stronger prediction.
3. **MW-alone AUROC on the panel as actually docked.** It was 0.6011 at screening. If the docked
   subset differs materially, the confound has shifted and the reading is qualified.
4. **Within-MW-quartile AUROCs — reported only where each quartile holds ≥20 of the minority
   class.** On PD-L1, quartiles 1–3 held 8, 1 and 4 inactives and produced AUROCs of 0.893, 0.990
   and 0.720 that meant nothing. A minimum-count rule prevents quoting them.
5. **Analysable n, not fetched n**, with prep failures and dock failures stated separately.

## What no outcome licenses

A general claim about docking. This would be three targets: Mpro below chance, PD-L1 above chance
on a confounded panel, and Factor Xa. Three points, one of them compromised, is not a
characterisation of a method — it is three measurements, each reportable on its own terms.

## Compute

Full panel, no cap. 1,006 compounds across 2 Kaggle CPU shards at ~12–20 compounds/hour, so
roughly 3 sessions per shard with resume cycles. The capped alternative (250/300, CI [0.382,
0.478] at 0.43) was rejected in favour of the tighter blind zone the full panel buys.

---

## Run log: why the analysis did not run at 802 of 894 (27 Aug 2026)

Both Kaggle shards hit the 12h session cap and were cancelled mid-run — shard-0 at
110/123 remaining-queue, shard-1 at 90/138 — leaving **802 of 894 ligands docked (89.7%)**.

It was tempting to analyse there. The composition of the unfinished tail says not to:

| | n | actives | inactives | active rate |
|---|---|---|---|---|
| cached | 802 | 382 | 420 | **47.6%** |
| remaining | 92 | 24 | 68 | **26.1%** |
| full set | 894 | 406 | 488 | 45.4% |

**The missing 92 are decoy-enriched.** The workers process ligands in a fixed order and
were cut inside a decoy-rich block — visible as an anomalously low 9th decile in both
shards (13% and 18% active against a ~45% base rate). So the completed subset is not a
random sample of the panel; it is active-enriched by construction.

This project has already measured how much that matters: gate admissibility **flipped in
both directions** under active subsampling on the earlier targets (Mpro qualified at 30
actives and failed at 85; PD-L1 failed at 30 and qualified at 82). An AUROC computed on a
non-randomly truncated panel is exactly the quantity that was shown to move.

Both shards were therefore resumed from their published caches (427 and 375, no shrink)
rather than analysed. Remaining: 24 on shard-0 (~2.5h), 68 on shard-1 (~8.5h) — both
inside one 12h session, so this should be the final cycle. The pre-registered rule above
is unchanged and will be applied to the complete 894.
