#!/bin/bash
# Wait for each Kaggle shard to finish, then run one relaunch cycle and merge.
#
# Sessions end ~11h apart and not together, so each shard is handled independently
# the moment it completes rather than waiting for both.
#
# Runs up to MAX_CYCLES per shard so coverage is continuous overnight and through
# the following day: a session finishing at 05:00 is relaunched immediately rather
# than idling until someone is awake. Bounded, not open-ended — it stops after
# MAX_CYCLES or when a shard's assignment is complete.
#
# Safety comes from the tools it calls, not from this script:
#   relaunch.py      refuses a running kernel; refuses an empty or SHRINKING cache
#   merge_shards.py  refuses on protocol mismatch, list mismatch, or calibration
#                    outside 0.3 kcal/mol — so an auto-merge cannot admit results
#                    from an environment that has drifted
#
# Logs everything to watch_relaunch.log for review.

cd "$(dirname "$0")" || exit 1
SP=${VENV_DIR:-./.venv}
PY=$SP/kgvenv/bin/python
DOCK=python
LOG=watch_relaunch.log
TID=CHEMBL612545

say() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a $LOG; }

MAX_CYCLES=${MAX_CYCLES:-4}
c0=0; c1=0
say "watcher started; up to $MAX_CYCLES cycles per shard"

for i in $(seq 1 200); do           # 200 * 15min = 50h ceiling
  for s in 0 1; do
    eval "c=\$c$s"
    [ "$c" -ge "$MAX_CYCLES" ] && continue

    st=$(env -u KAGGLE_USERNAME -u KAGGLE_KEY $PY relaunch.py --shard $s --status 2>/dev/null \
         | grep -o 'status = [a-zA-Z]*' | cut -d' ' -f3)
    if [ -z "$st" ]; then
      say "shard$s: status query failed, will retry"
      continue
    fi
    if [ "$st" = "running" ]; then
      continue
    fi

    say "shard$s: status=$st -> running relaunch cycle"
    env -u KAGGLE_USERNAME -u KAGGLE_KEY $PY relaunch.py --shard $s >> $LOG 2>&1
    rc=$?
    say "shard$s: relaunch exit=$rc"

    if [ $rc -eq 0 ]; then
      say "shard$s: merging (calibration gate enforced by merge_shards)"
      $DOCK merge_shards.py --tid $TID vina_cache_${TID}_shard${s}.json --apply >> $LOG 2>&1
      say "shard$s: merge exit=$?"
      n=$($DOCK -c "import json;print(len(json.load(open('../vina_cache_${TID}.json'))))" 2>/dev/null)
      say "shard$s: main cache now $n/2526"
    else
      say "shard$s: relaunch did not succeed — NOT merging. Inspect $LOG."
    fi
    eval "c$s=\$((c$s + 1))"
    eval "c=\$c$s"
    say "shard$s: cycle $c of $MAX_CYCLES complete"

    # stop early if this shard's assignment is finished
    remain=$($DOCK -c "
import json
m=json.load(open('bundle_${TID}/shard${s}/manifest.json'))
c=json.load(open('vina_cache_${TID}_shard${s}.json'))
print(len(m['compounds'])-len(c))" 2>/dev/null)
    if [ -n "$remain" ] && [ "$remain" -le 0 ]; then
      say "shard$s: assignment COMPLETE, no further cycles"
      eval "c$s=$MAX_CYCLES"
    fi
  done

  if [ "$c0" -ge "$MAX_CYCLES" ] && [ "$c1" -ge "$MAX_CYCLES" ]; then
    say "both shards reached their cycle limit; watcher exiting"; exit 0
  fi
  sleep 900
done

say "watcher hit its time ceiling; cycles done: shard0=$c0 shard1=$c1"
