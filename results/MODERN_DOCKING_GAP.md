# Have we tested DiffDock or any modern deep-learning docking? No — and it has been named as the fix twice

> **Status — read `../RESULTS_STATUS.md` first.** This document is the record of the analysis on the date it carries. Answered since: Boltz-2 was tested, and combined with seven descriptors it clears the pre-registered margin on two targets (+0.0934 and +0.0751) where docking does not.

**Assessed 31 August 2026.** Prompted by the question directly; recorded because the answer is
"no" and the reason it stayed "no" is worth having on record.

## What the repo actually contains

`grep -rli diffdock|equibind|gnina|boltz|neuralplexer` returns **`paper_draft.md` and nothing
else**. No code, no weights, no install, no results. It appears three times, always as a
recommendation:

> *"requires a protein-aware oracle (e.g., DiffDock, EquiBind) to draw valid conclusions"* — §4.3
> *"A protein-aware surrogate (e.g., DiffDock, EquiBind, or a pocket-conditioned GNN) would be
> required to discriminate between targets."* — §4.7

Named as the necessary next step across two phases, never executed.

## What we HAVE tested, and why it does not count

A PDBbind-trained **ligand-only** GNN, on the CHEMBL612545 panel — and that result is
**withdrawn**. The "target" was a catch-all ChEMBL bucket with 0 target components, so there were
no true actives to enrich. The GNN's AUROC 0.677 on a set with no signal is evidence it reads
property and scaffold structure rather than binding. That is the opposite of validating a learned
docking method.

So: no structure-based deep-learning docking method has ever been run here.

## Where that leaves us on Mpro (753 compounds, measured actives and inactives)

| predictor | AUROC |
|---|---|
| Vina, exhaustiveness 32 | **0.427** |
| Vina, exhaustiveness 4 | 0.408 |
| chance | 0.500 |
| **seven physicochemical descriptors, 5-fold CV logistic** | **0.763** |

**Docking is 0.34 AUROC worse than seven descriptors and a logistic regression.** Any new method
has to clear **0.763**, not 0.5. That is the bar, and it is already in the repo.

## The obvious explanation is wrong — checked, not assumed

Docking score tracks molecular size, so the natural read is that actives here are small and the
size bias ranks them last. Measured on the same 753:

| | |
|---|---|
| actives vs inactives, heavy atoms | **+2.3** (actives are *bigger*) |
| corr(heavy atoms, Vina score) | −0.588 (bigger scores "better") |
| corr(heavy atoms, is_active) | +0.134 |
| **corr(Vina score, is_active)** | **+0.104** (worse score → more likely active) |

Size should push docking to rank actives **better** on both legs, and it still ranks them worse.
This is a genuine anti-signal on this target, not a size confound. Which makes the method question
sharper, not softer.

## But DiffDock is probably not the fix, and the literature says why

2025–26 benchmarks are consistent: **DiffDock's own confidence score is near-random for virtual
screening** — AUROC 0.538 and 0.56, enrichment factors below 1.0 at 1–10%. It is a pose predictor,
and its confidence is not an affinity estimate. Swapping our 0.427 for a method whose screening
signal is ~0.54 does not clear 0.763.

What *does* work in those benchmarks is **decoupling pose generation from scoring**: DiffDock-L
poses rescored with a physics function beat Vina's own pose sampling (BEDROC 0.22 vs 0.10, EF@1%
10.92 vs 3.89), and DiffDock-L + gnina reaches BEDROC 0.33–0.36.

## The experiment that would actually be worth running

Not "DiffDock instead of Vina". Three arms on the same 753, against the 0.763 bar:

1. **Vina pose + Vina score** — done, 0.427.
2. **DiffDock-L pose + Vina rescore** — tests whether our failure is in *sampling* or in *scoring*.
3. **gnina rescore** of both pose sets — a CNN scoring function built for screening, unlike a
   diffusion confidence head.

If arm 2 moves and arm 3 moves more, the failure is scoring. If neither moves, the failure is the
target or the benchmark, and no docking method fixes it. Either answer is worth having, and arm 1
already exists.

Needs a GPU (Kaggle, 9h cap). DiffDock's dependency stack (torch-geometric, ESM) is the fiddly
part, and is the reason to scope this as its own run rather than bolt it onto an existing one.

## Sources

- [Benchmarking Single-Pose Docking, Consensus Rescoring, and Supervised ML on LIT-PCBA](https://arxiv.org/pdf/2605.01681)
- [Physics beats diffusion: agentic AI-driven virtual screening on a GPCR target](https://www.researchsquare.com/article/rs-9142847/v1)
- [Integrating ML-Based Pose Sampling with Established Scoring Functions (JCIM)](https://pubs.acs.org/doi/10.1021/acs.jcim.5c00380)
