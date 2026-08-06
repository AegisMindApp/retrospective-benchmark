#!/bin/bash
# Queue the Mpro re-run to start once PD-L1 finishes.
#
# Both targets must be re-docked under Amendment 19 (all prior docking was rigid-
# ligand). Running them concurrently on 4 cores would roughly halve each, so they
# are serialised. Mpro is the subsample-qualified case: it passes the gate only at
# MAX_ACTIVES=30 and fails at 85, so it is reported alongside PD-L1 rather than
# instead of it.
cd "$(dirname "$0")"
D=python

echo "[queue] waiting for PD-L1 to finish..."
while pgrep -f panel_pdl1_only >/dev/null; do sleep 300; done
echo "[queue] PD-L1 ended $(date -Is)"

# Final harvest of the PD-L1 cache before moving on.
python3 vina_cache.py pdl1_exh32_flex.log vina_cache_CHEMBL612545.json || true

echo "[queue] starting Mpro re-run $(date -Is)"
BENCH_EXHAUSTIVENESS=32 nohup $D run_pilot_fast.py panel_mpro_only.json \
  --out-prefix mpro_exh32_flex --workers 3 > mpro_exh32_flex.log 2>&1 &
MPID=$!
echo "[queue] Mpro pid $MPID"

# Harvest Mpro's log too, so it inherits the same 30-minute loss ceiling.
while kill -0 $MPID 2>/dev/null; do
  python3 vina_cache.py mpro_exh32_flex.log vina_cache_CHEMBL4523582.json >/dev/null 2>&1
  sleep 1800
done
python3 vina_cache.py mpro_exh32_flex.log vina_cache_CHEMBL4523582.json
echo "[queue] Mpro finished $(date -Is)"
