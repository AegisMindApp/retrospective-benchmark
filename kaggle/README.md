# Kaggle shard harness

Splits a docking run across the local machine and two Kaggle CPU sessions.
Authorised by **Amendment 20**; read that first — it records the pre-committed
calibration threshold and the stopping rule.

## The one design rule

**Nothing that affects the protocol is reproduced remotely.**

Ligand PDBQTs are built locally with this machine's RDKit 2026.03.3 and Meeko 0.7.1
and shipped. The Vina 1.2.7 binary and the receptor are shipped byte-for-byte and
md5-verified on arrival. Kaggle runs `vina` and nothing else — no RDKit, no Meeko,
no version negotiation.

That is not caution for its own sake. Meeko and RDKit generate the torsion tree; an
RDKit bump alone moved gate AUROC 0.585 → 0.598 here, and Amendment 19 exists
because ligand prep silently differed from what was assumed. Pinning versions would
be a claim. Shipping artifacts is a fact.

## Sequence

```bash
cd analysis/retrospective_benchmark/kaggle
D=python

# 1. Freeze the compound list and verify it against the live run.
#    Exits non-zero if ChEMBL has moved — then the target is NOT shardable.
$D freeze_panel.py ../panel_pdl1_only.json CHEMBL612545 \
     --verify-cache ../vina_cache_CHEMBL612545.json

# 2. Build two bundles from the TAIL of the remaining compounds.
$D make_bundle.py --tid CHEMBL612545 --shards 2

# 3. Upload bundle_CHEMBL612545/shard0 and shard1 as two Kaggle Datasets.
#    In each notebook (CPU, internet off is fine — nothing is downloaded):
#       !python /kaggle/input/<ds>/kaggle_shard_worker.py \
#            --bundle /kaggle/input/<ds> --workers 4
#    Save/commit the notebook so /kaggle/working persists.

# 4. Resume a session: attach the previous run's output as an input dataset and
#    re-run the same cell. The worker picks up its own cache.

# 5. Merge — dry run first. Refuses on protocol, list or calibration failure.
$D merge_shards.py --tid CHEMBL612545 shard0_cache.json shard1_cache.json
$D merge_shards.py --tid CHEMBL612545 shard*.json --apply

# 6. Should the local run stop? This is the ONLY authority for that decision.
$D check_coverage.py --tid CHEMBL612545 shard*.json --log ../pdl1_exh32_flex.log
```

## Relaunching between sessions

A Kaggle session dies at 12h, so each shard takes several runs. One command per
cycle:

```bash
$D relaunch.py --shard 0 --status   # look, change nothing
$D relaunch.py --shard 0            # fetch results, publish, re-push
```

It refuses to act while the kernel is still running, and refuses to publish an
empty or *shrinking* cache — a progress dataset that went backwards would erase
completed work at the next resume.

Each cycle: download the finished session's cache → **write a local copy first**
(that is what `merge_shards.py` gates on, and it survives anything happening to the
Kaggle side) → publish it as `pdl1-shard-N-progress` → re-push with bundle +
progress attached. The worker merges every cache it finds under `/kaggle/input`
rather than taking the first, because with more than one prior dataset attached,
first-match would silently drop a session's work.

## Why shards come from the tail

The local process runs pre-fix code and will dock all 2,526 compounds whatever the
remote workers return. Local works the head, Kaggle the tail; they meet in the
middle. Duplicated compounds are waste, not error — the local score wins. The run is
done when the **union** covers the frozen list, which is what `check_coverage.py`
computes and no one should be computing by hand at 3am.

## The gate that actually matters

`merge_shards.py` runs three checks — protocol md5s, compound-list md5, and
cross-environment calibration. The third is the one that would catch a genuinely
different remote environment: five compounds already scored locally are re-docked on
Kaggle and must agree within **0.3 kcal/mol**. They will not agree exactly, because
Vina runs without `--seed`.

If calibration fails, the shard is rejected **in full**. Per Amendment 20 the
distributed plan is then abandoned and the run finishes locally. Do not lower the
tolerance to make a shard pass.

## Expected arithmetic

| | rate | notes |
|---|---|---|
| local | 8.85 cmpd/h | 3 workers, i7-6500U, load-bound |
| Kaggle × 2 | ~30 cmpd/h | 4 cores each, 12h sessions |

~2,291 remaining at a combined ~40/h is roughly 2.5 days of session time versus
10.8 days local-only — but Kaggle sessions need restarting every 12h, so wall-clock
depends on how often they are relaunched.

## Files

| File | Role |
|---|---|
| `freeze_panel.py` | reproduce and freeze the exact compound list; verify vs live run |
| `make_bundle.py` | prepare ligands locally, assemble shard bundles + calibration |
| `kaggle_shard_worker.py` | remote worker; dependency-free, resumable, budget-capped |
| `merge_shards.py` | three-gate merge; dry run by default |
| `check_coverage.py` | union coverage; the stopping rule |
