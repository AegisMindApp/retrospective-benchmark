#!/bin/bash
# Continuously harvest completed docking results from the live log into the cache,
# so an interruption at any point loses at most the last 30 minutes rather than
# the whole run. Safe to run alongside the docking process — read-only on the log,
# atomic write on the cache.
cd "$(dirname "$0")"
while pgrep -f panel_pdl1_only >/dev/null; do
  python3 vina_cache.py pdl1_exh32_flex.log vina_cache_CHEMBL612545.json >/dev/null 2>&1
  sleep 1800
done
python3 vina_cache.py pdl1_exh32_flex.log vina_cache_CHEMBL612545.json
