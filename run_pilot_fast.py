#!/usr/bin/env python3
"""
run_pilot_fast.py — parallelised pilot: first 3 qualifying targets, exhaustiveness=4.

Identical to run_pilot.py in its statistical logic (same metrics, same blind→unblind
order, same pre-registration). Differences from the registered workflow are practical,
not analytical: lower Vina exhaustiveness (4 vs. default 16) and concurrent workers.
Both vina and platform use the same exhaustiveness so the comparison is still fair.

NOT a separate pre-registration — this is the pilot (§5 of PREREGISTRATION.md: first 3
qualifying targets). Full-exhaustiveness production run is a follow-up amendment.
"""
from __future__ import annotations

import asyncio
import json
import os
import random
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH.parent.parent))   # repo root (for aegismind imports)

os.environ.setdefault("OBABEL", "obabel")

from metrics import score_ranking, bootstrap_ci, paired_panel_test, EF_FRACTIONS
from baselines import chance_scores, fp_similarity_scores
import decoys

SHUFFLE_SEED    = 0
PRIMARY_METRIC  = "EF@1%"
PILOT_TARGETS   = 5           # high-diversity panel: PD-L1, AChE, Aldose reductase, MDM2 + Mpro anchor
EXHAUSTIVENESS  = int(os.environ.get('BENCH_EXHAUSTIVENESS', 32))  # Amendment 16: 32 (was 4; registered default 16)
N_WORKERS       = int(os.environ.get("VINA_WORKERS", max(1, os.cpu_count() - 1)))  # override: VINA_WORKERS=1


# ── parallel Vina ─────────────────────────────────────────────────────────────
def _dock_one(args):
    smi, name, receptor_pdbqt, box_center, box_size, exhaustiveness = args
    import sys, os
    # Path(__file__) is unreliable in ProcessPoolExecutor subprocesses; use
    # module-level BENCH (absolute, resolved at import time) instead.
    _root = str(BENCH.parent.parent)
    for p in (_amr, _root):
        if p not in sys.path:
            sys.path.insert(0, p)
    os.environ.setdefault("OBABEL", "obabel")
    import vina_tools as rvd
    rvd.BOX_CENTER = tuple(box_center)
    rvd.BOX_SIZE   = tuple(box_size)
    lig = str(rvd.RESULTS_DIR / f"bench_{name}.pdbqt")
    if not rvd.smiles_to_pdbqt(smi, f"bench_{name}", lig):
        return float("nan")
    res = rvd.run_vina(receptor_pdbqt, lig, f"bench_{name}", exhaustiveness=exhaustiveness)
    aff = res.get("best_affinity_kcal_mol")
    return -float(aff) if aff is not None else float("nan")


def _cache_path(out_prefix: str, tid: str) -> str:
    return f"vina_cache_{tid}.json"


def parallel_vina_scores(smiles_list: Sequence[str], receptor_pdbqt: str,
                         box_center, box_size,
                         names: Optional[Sequence[str]] = None,
                         exhaustiveness: int = EXHAUSTIVENESS,
                         cache_path: Optional[str] = None) -> List[float]:
    """Dock every compound, reusing any already present in `cache_path`.

    Resume (added 5 Aug 2026). A flexible-ligand run is ~10 days; without this a
    single interruption discards all of it. Cached compounds are skipped entirely
    and every new result is written the moment it lands, so the worst case is the
    loss of whatever was in flight rather than the whole run.

    Note the compound key is the ligand NAME, matching the '[vina] NAME: best pose'
    line the log recovery parses, so a cache built from a log and a cache built
    live are interchangeable.
    """
    import vina_cache as _vc
    names = names or [f"cmpd_{i}" for i in range(len(smiles_list))]
    cached = _vc.load(cache_path) if cache_path else {}

    scores = [float("nan")] * len(smiles_list)
    todo = []
    for i, (smi, nm) in enumerate(zip(smiles_list, names)):
        key = f"bench_{nm}"
        if key in cached:
            scores[i] = cached[key]                      # reuse, do not re-dock
        elif nm in cached:
            scores[i] = cached[nm]
        else:
            todo.append((i, (smi, nm, receptor_pdbqt, list(box_center),
                             list(box_size), exhaustiveness)))

    reused = len(smiles_list) - len(todo)
    if reused:
        print(f"  [resume] reusing {reused} cached results; docking {len(todo)}")

    if not todo:
        return scores

    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(_dock_one, a): i for i, a in todo}
        done = 0
        for fut in as_completed(futs):
            i = futs[fut]
            try:
                scores[i] = fut.result()
                if cache_path and scores[i] == scores[i]:      # not NaN
                    _vc.record(f"bench_{names[i]}", scores[i], cache_path)
            except Exception as e:
                print(f"  [dock] {names[i]}: error {e}")
            done += 1
            if done % 20 == 0:
                print(f"  [dock] {done}/{len(todo)} done ({reused} reused)")
    return scores


# ── ensemble scoring ───────────────────────────────────────────────────────────
def _zscore(x):
    a = np.asarray(x, dtype=float)
    m, s = np.nanmean(a), np.nanstd(a)
    a = np.where(np.isnan(a), m, a)
    return (a - m) / s if s > 0 else np.zeros_like(a)


def _load_api_key():
    """Read ANTHROPIC_API_KEY from env or .env file."""
    import os
    key = os.environ.get("ANTHROPIC_API_KEY")
    if key:
        return key
    env_file = BENCH.parent.parent / ".env"
    if env_file.exists():
        for line in env_file.read_text().splitlines():
            line = line.strip()
            if line.startswith("ANTHROPIC_API_KEY="):
                return line.split("=", 1)[1].strip().strip('"').strip("'")
    return None


def _score_one_sync(args):
    """Single Anthropic API call (runs in thread pool). Returns float 0–1."""
    import re
    import anthropic as _ant
    api_key, smi, target, context = args
    client = _ant.Anthropic(api_key=api_key)
    prompt = (
        "Rate the probability (0.0–1.0) that this compound is a potent sub-µM binder "
        f"of {target}. Context: {context}. SMILES: {smi}. "
        "Reply with ONLY a number then a 5-word reason."
    )
    try:
        resp = client.messages.create(
            model="claude-haiku-4-5-20251001",
            max_tokens=80,
            messages=[{"role": "user", "content": prompt}],
        )
        text = resp.content[0].text if resp.content else ""
        m = re.search(r"(0?\.\d+|[01]\.?\d*)", text)
        if m:
            v = float(m.group(1))
            if 0.0 <= v <= 1.0:
                return v
    except Exception:
        pass
    return 0.5


async def ensemble_activity_scores(smiles_list, target, context, top_n=150):
    """Score top_n compounds via direct Anthropic API (claude-haiku), 8 concurrent threads."""
    import asyncio
    from concurrent.futures import ThreadPoolExecutor

    api_key = _load_api_key()
    scores = [0.5] * len(smiles_list)
    if not api_key:
        print("  [ensemble] No ANTHROPIC_API_KEY found — returning 0.5 defaults")
        return scores

    loop = asyncio.get_event_loop()
    subset = smiles_list[:top_n]
    args = [(api_key, smi, target, context) for smi in subset]

    print(f"  [ensemble] scoring {len(subset)} compounds via Anthropic API (8 threads)...")
    with ThreadPoolExecutor(max_workers=8) as pool:
        futs = [loop.run_in_executor(pool, _score_one_sync, a) for a in args]
        results = await asyncio.gather(*futs)

    for i, v in enumerate(results):
        scores[i] = v
    n_scored = sum(1 for v in results if v != 0.5)
    print(f"  [ensemble] {n_scored}/{len(subset)} returned non-default scores")
    return scores


# ── main ───────────────────────────────────────────────────────────────────────
async def run(panel_path="panel_frozen.json", receptors_path="receptors.json",
              out_prefix="pilot_fast"):
    panel_all = json.load(open(panel_path))["targets"]
    receptors = json.load(open(receptors_path)) if Path(receptors_path).exists() else {}

    targets = [t for t in panel_all if t in receptors][:PILOT_TARGETS]
    print(f"Pilot targets ({len(targets)}): {targets}")
    print(f"Workers: {N_WORKERS}, exhaustiveness: {EXHAUSTIVENESS}\n")

    per_target_metrics: Dict[str, Dict] = {}
    decoy_report: Dict[str, dict] = {}
    skipped: List[str] = []

    for tid in targets:
        comps = panel_all[tid]
        rec   = receptors[tid]
        print(f"\n{'='*60}\n{tid} — {rec['target_name']}")

        actives      = [c for c in comps if c["label"] == 1]
        real_inact   = [c for c in comps if c["label"] == 0]

        print(f"  Building decoys ({len(actives)} actives → max {decoys.MAX_ACTIVES})...")
        dres = decoys.build_decoys(
            [{"molecule_chembl_id": a["molecule_chembl_id"], "smiles": a["smiles"]} for a in actives], tid)
        decoy_report[tid] = {k: dres.get(k) for k in
                             ("status","n_actives","n_decoys","debias_auroc","scaffold_count")}
        print(f"  Decoys: {dres['status']}, n={dres.get('n_decoys',0)}, auroc={dres.get('debias_auroc'):.3f}" if dres.get('debias_auroc') else f"  Decoys: {dres['status']}")
        if dres["status"] != "ok":
            skipped.append(f"{tid}:decoys({dres['status']})")
            continue

        scored = (
            [{"smiles": a["smiles"], "molecule_chembl_id": a["molecule_chembl_id"], "label": 1}
             for a in dres["actives"]]
            + [{"smiles": d["smiles"], "molecule_chembl_id": d["molecule_chembl_id"], "label": 0}
               for d in dres["decoys"]]
            + [{"smiles": c["smiles"], "molecule_chembl_id": c["molecule_chembl_id"], "label": 0}
               for c in real_inact]
        )
        random.Random(SHUFFLE_SEED).shuffle(scored)
        smiles = [c["smiles"] for c in scored]
        labels = [c["label"] for c in scored]
        names  = [c["molecule_chembl_id"] for c in scored]

        print(f"  Scoring {len(smiles)} compounds (chance + fp)...")
        methods: Dict[str, List[float]] = {}
        methods["chance"] = chance_scores(len(smiles), seed=0)
        try:
            methods["fp_similarity"] = fp_similarity_scores(smiles, tid)
        except Exception as e:
            skipped.append(f"{tid}:fp({e})")

        print(f"  Vina docking {len(smiles)} compounds with {N_WORKERS} workers...")
        try:
            # Checkpoint/resume (added 5 Aug 2026). Docking is a multi-day blocking
            # call that previously wrote nothing until it finished, so any
            # interruption discarded the whole run. Consult the cache first and
            # write each result as it lands.
            import vina_cache as _vc
            import vina_tools as _rvd
            _cache_file = _cache_path(out_prefix, tid)

            # Refuse to resume across a protocol change. The cache is keyed on
            # target alone, so without this an exhaustiveness-8 or rigid-ligand
            # cache would be reused wholesale under a different protocol.
            _proto = {"exhaustiveness": EXHAUSTIVENESS,
                      "num_modes": 5,
                      "receptor": os.path.basename(rec["receptor_pdbqt"]),
                      # md5s, not version strings. A version string is a claim about
                      # the binary; an md5 is a fact about it. vina_tools.py
                      # hardcoded "vina_version": "1.2.5" while the binary was 1.2.7.
                      "vina_md5": _vc.file_md5(_rvd.VINA_BIN),
                      "receptor_md5": _vc.file_md5(rec["receptor_pdbqt"]),
                      # rounded: a recomputed centroid can differ in the last
                      # float digit between runs, and a false "mismatch" that
                      # blocks a legitimate 3am resume is worse than no guard
                      "box_center": [round(float(v), 2) for v in rec["box_center"]],
                      "box_size": [round(float(v), 2) for v in rec["box_size"]],
                      "flexible_ligand": True}
            _status, _stored = _vc.check_protocol(_cache_file, _proto)
            if _status == "mismatch":
                raise RuntimeError(
                    f"cache {_cache_file} was built under a DIFFERENT protocol "
                    f"({_stored}) than this run ({_proto}). Refusing to mix. "
                    f"Move the cache aside to start clean.")
            if _status == "legacy":
                print(f"  [resume] cache has no protocol record; adopting current "
                      f"protocol. Verify it really was exh={EXHAUSTIVENESS} flexible.")
                _vc.stamp(_cache_file, _proto)

            # Merge log and cache unconditionally, not only when the cache is
            # empty: a cache that exists but lags the log (harvest died mid-run)
            # would otherwise be used as-is and the extra results re-docked.
            _cached = {**_vc.recover_from_log(f"{out_prefix}.log"),
                       **_vc.load(_cache_file)}
            if _cached:
                _vc.save(_cached, _cache_file)
                print(f"  [resume] {len(_cached)} of {len(smiles)} already docked; "
                      f"{len(smiles)-len(_cached)} remaining")

            vina = parallel_vina_scores(
                smiles, rec["receptor_pdbqt"], rec["box_center"], rec["box_size"], names,
                cache_path=_cache_file)
            methods["vina"] = vina
        except Exception as e:
            skipped.append(f"{tid}:vina({e})")
            json.dump({m: s for m, s in methods.items()}, open(f"{out_prefix}_predictions_{tid}.json","w"))
            continue

        # ensemble: reorder by vina rank, score top 150
        order = np.argsort(-np.where(np.isnan(vina), -999, vina))
        smiles_ranked = [smiles[i] for i in order]
        try:
            ens_ranked = await ensemble_activity_scores(smiles_ranked, rec["target_name"], rec["context"])
            ens = [0.5] * len(smiles)
            for rank_i, orig_i in enumerate(order):
                ens[orig_i] = ens_ranked[rank_i]
            z_vina = _zscore(vina)
            z_ens  = _zscore(ens)
            methods["platform"] = list(z_vina + z_ens)
        except Exception as e:
            skipped.append(f"{tid}:ensemble({e})")

        # blind predictions persisted before unblinding
        json.dump({m: s for m, s in methods.items()},
                  open(f"{out_prefix}_predictions_{tid}.json", "w"))

        per_target_metrics[tid] = {}
        for m, s in methods.items():
            met = score_ranking(s, labels)
            p, lo, hi = bootstrap_ci(s, labels, PRIMARY_METRIC)
            met[f"{PRIMARY_METRIC}_ci"] = [round(p,3), round(lo,3), round(hi,3)]
            per_target_metrics[tid][m] = met
            print(f"  {m:15s} EF@1%={met['EF@1%']:.2f}  AUROC={met['AUROC']:.3f}")

    # panel test
    report = {"primary_metric": PRIMARY_METRIC, "exhaustiveness": EXHAUSTIVENESS,
              "n_workers": N_WORKERS, "per_target": per_target_metrics,
              "decoy_report": decoy_report, "skipped": skipped, "panel_tests": {}}
    def col(method):
        return {t: per_target_metrics[t][method][PRIMARY_METRIC]
                for t in per_target_metrics if method in per_target_metrics[t]}
    if any("platform" in per_target_metrics[t] for t in per_target_metrics):
        report["panel_tests"]["platform_vs_vina"]   = paired_panel_test(col("platform"), col("vina"))
        report["panel_tests"]["platform_vs_chance"] = paired_panel_test(col("platform"), col("chance"))
        report["panel_tests"]["vina_vs_chance"]     = paired_panel_test(col("vina"),     col("chance"))

    out_path = f"{out_prefix}_report.json"
    json.dump(report, open(out_path, "w"), indent=2)
    print(f"\n{'='*60}")
    print(f"Report written: {out_path}")
    print(json.dumps({"skipped": skipped, "panel_tests": report["panel_tests"]}, indent=2))
    return report


if __name__ == "__main__":
    import argparse as _ap
    _p = _ap.ArgumentParser()
    _p.add_argument("panel", nargs="?", default="panel_frozen.json")
    _p.add_argument("--out-prefix", default="pilot_fast")
    _p.add_argument("--workers", type=int, default=None,
                    help="Override N_WORKERS (default: cpu_count-1)")
    _args = _p.parse_args()
    if _args.workers is not None:
        N_WORKERS = _args.workers
    asyncio.run(run(panel_path=_args.panel, out_prefix=_args.out_prefix))
