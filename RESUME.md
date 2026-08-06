# Resuming an interrupted docking run

A flexible-ligand run at exhaustiveness 32 takes ~11 days. Losing it to a reboot,
an OOM kill or a stray `pkill` would be expensive, so resume is now wired and
tested rather than assumed.

## How to restart after a glitch

```bash
cd analysis/retrospective_benchmark
conda activate docking
nohup python run_pilot_fast.py panel_pdl1_only.json \
      --out-prefix pdl1_exh32_flex --workers 3 >> pdl1_exh32_flex.log 2>&1 &
```

Use the **same `--out-prefix`**. That is what names the log and the cache; change
it and the run starts from zero. Append (`>>`), never truncate (`>`) — the log is
a recovery source in its own right.

On start you should see:

```
  [resume] reusing 235 cached results; docking 2291
```

If that line is missing, or the number is 0, stop and diagnose before letting it
run — it means the cache was not found and you are about to repeat 26 hours.

## What protects the work

Two independent layers, either of which is sufficient:

1. **`vina_cache_<target>.json`** — every completed dock is written the moment it
   lands, keyed `bench_<CHEMBL_ID>`.
2. **The log** — `recover_from_log()` re-parses `[vina] NAME: best pose N kcal/mol`
   lines.

On start the two are merged **unconditionally** (`{**from_log, **from_cache}`),
not only when the cache is empty. A cache that exists but *lags* the log — harvest
died at entry 100 while the log reached 235 — would otherwise be taken at face
value and the extra 135 re-docked. `load()` returns `{}` on unreadable JSON, so a
deleted or corrupted cache is rebuilt from the log rather than crashing.

Verified on the live PD-L1 run: log and cache agreed on all 235 entries with zero
mismatches. The loss ceiling is whatever is in flight — at most `--workers`
compounds, so ~20 minutes.

### Do not stop `harvest_cache.sh` while PD-L1 is running

The live PD-L1 process (started 4 Aug) is executing **pre-fix code that never
writes the cache itself**. Its 235 entries exist *only* because harvest scrapes
them out of the log. Killing harvest removes the sole thing checkpointing the run
that has 26 hours in it.

Harvest becomes genuinely redundant only for runs started with the current code,
which record directly. It stays harmless either way: writes are atomic via
`os.replace`, and since it rebuilds from the append-only log it can only restore
entries, never lose them.

**Do not restart PD-L1 to pick up the fix.** There is no gain — the fix only
matters at restart, and log recovery already covers that — against a real chance
of losing in-flight work.

## Protocol guard

The cache filename is keyed on **target only**: nothing in `vina_cache_CHEMBL612545.json`
distinguishes exhaustiveness 8 from 32, or rigid ligands from flexible. Re-running a
target under a changed protocol would silently reuse the old scores and blend two
protocols into one result set — which is precisely the Amendment 19 failure.

A sidecar `vina_cache_<target>.protocol.json` now records exhaustiveness, receptor,
box and the flexible-ligand flag. On start the run compares and **aborts on
mismatch** rather than resuming. Failing loudly beats quietly re-docking from zero.

The sidecar is deliberately a separate file: `harvest_cache.sh` and
`queue_mpro_after_pdl1.sh` merge flat `{name: score}` dicts through the
`vina_cache.py` CLI, and the Mpro one fires unattended in ~11 days. Restructuring
the scores JSON under a running shell loop was not worth the risk.

The live PD-L1 cache has been stamped with its true protocol (exh 32, flexible,
`CHEMBL612545_receptor.pdbqt`, box 22³ at 13.9/6.07/183.23), taken from the running
vina command line — so a restart sees `ok`, not `legacy`.

## The bug this replaced

Until 5 Aug 2026 the cache was loaded and a `[resume]` line was printed, but
`parallel_vina_scores` was still called with the full compound list. Nothing was
skipped and the cached scores were discarded — the run *reported* resuming and
then silently re-docked everything from the start.

Guard against regression:

```bash
python test_resume.py   # 9 assertions; exits non-zero on regression
```

The test monkeypatches `_dock_one` and checks that cached compounds are not
re-docked, that exactly the uncached ones are, that cached scores land at the
correct positions in the original ordering, that new results reach the cache, and
that the protocol guard refuses an exhaustiveness change, refuses rigid-vs-flexible,
and reports an unstamped cache as `legacy` rather than `ok`. Note that
`ProcessPoolExecutor` forks, so the ledger of docked compounds must be written to
a file — an in-memory list in the child is invisible to the parent and will make
the test look like it passed when it proved nothing.
