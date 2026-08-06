"""Checkpoint / resume for the docking arm.

run_pilot_fast.py docks every compound in one blocking call and writes nothing
until it finishes. At ~42 min/compound for flexible ligands, a 2,526-compound
target is a ~25-day run with a single point of failure at the end. A power cut at
day 24 loses everything.

Two capabilities:

  recover_from_log(log)   parse completed "[vina] NAME: best pose X" lines out of
                          a run log into a cache dict. Lets an interrupted run's
                          work be reclaimed after the fact.

  load / save / merge     a {name: score} JSON cache the runner can consult before
                          docking and append to after each result.

Usage after an interrupted run:
    python3 vina_cache.py pdl1_exh32_flex.log vina_cache_CHEMBL612545.json
"""
from __future__ import annotations
import hashlib, json, os, re, sys

LINE = re.compile(r"\[vina\]\s+(\S+):\s+best pose\s+(-?\d+\.?\d*)\s+kcal/mol")


def recover_from_log(log_path: str) -> dict[str, float]:
    """Extract every completed docking result from a run log."""
    out: dict[str, float] = {}
    if not os.path.exists(log_path):
        return out
    with open(log_path, errors="replace") as fh:
        for line in fh:
            m = LINE.search(line)
            if m:
                out[m.group(1)] = float(m.group(2))
    return out


def load(cache_path: str) -> dict[str, float]:
    try:
        with open(cache_path) as fh:
            return json.load(fh)
    except Exception:
        return {}


def save(cache: dict[str, float], cache_path: str) -> None:
    tmp = cache_path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(cache, fh)
    os.replace(tmp, cache_path)          # atomic; survives a kill mid-write


def record(name: str, score: float, cache_path: str) -> None:
    """Append one result and flush. Called per compound, so nothing is ever lost."""
    c = load(cache_path)
    c[name] = score
    save(c, cache_path)


# ── protocol guard ─────────────────────────────────────────────────────────────
# The cache filename is keyed on TARGET only, so nothing in it distinguishes
# exhaustiveness 8 from 32, or rigid ligands from flexible. Re-running a target
# under a changed protocol would silently reuse the old scores and blend two
# protocols into one result set — which is exactly the Amendment 19 failure.
#
# The protocol is recorded in a SIDECAR file rather than inside the scores JSON,
# so that harvest_cache.sh / queue_mpro_after_pdl1.sh (which merge flat dicts via
# the CLI below) keep working untouched.

def file_md5(path: str) -> str:
    """Identity of an artifact, as a fact rather than a claimed version string."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def protocol_path(cache_path: str) -> str:
    return cache_path[:-5] + ".protocol.json" if cache_path.endswith(".json") \
        else cache_path + ".protocol.json"


def check_protocol(cache_path: str, protocol: dict) -> tuple[str, dict | None]:
    """Compare `protocol` against the one this cache was built under.

    Returns (status, stored):
      'new'      no cache yet -> sidecar written, safe to proceed
      'ok'       sidecar matches -> safe to resume
      'legacy'   scores exist with no sidecar -> caller must decide, then stamp()
      'mismatch' sidecar disagrees -> caller MUST abort, not silently re-dock
    """
    sp = protocol_path(cache_path)
    stored = load(sp) or None
    if stored is None:
        if load(cache_path):
            return "legacy", None
        stamp(cache_path, protocol)
        return "new", protocol
    return ("ok" if stored == protocol else "mismatch"), stored


def stamp(cache_path: str, protocol: dict) -> None:
    """Record the protocol a cache was built under."""
    save(protocol, protocol_path(cache_path))


if __name__ == "__main__":
    log = sys.argv[1] if len(sys.argv) > 1 else "pdl1_exh32_flex.log"
    out = sys.argv[2] if len(sys.argv) > 2 else "vina_cache_recovered.json"
    rec = recover_from_log(log)
    existing = load(out)
    merged = {**existing, **rec}
    save(merged, out)
    print(f"recovered {len(rec)} results from {log}")
    print(f"cache now holds {len(merged)} -> {out}")
