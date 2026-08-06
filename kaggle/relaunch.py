"""One command to run the next Kaggle session for a shard.

Kaggle caps a session at 12h, so a shard takes several runs. Each cycle is:

    finished kernel -> download its cache -> keep a local copy (for merging)
                    -> publish it as a small "progress" dataset
                    -> re-push the kernel with bundle + progress attached
                    -> worker resumes from the cache and continues

Doing that by hand is roughly ten dataset wirings across both shards, each an
opportunity to attach the wrong version and silently lose a session's work.

The local copy is the point as much as the relaunch: it is what merge_shards.py
gates on, and it means a shard's results survive even if the Kaggle side is lost.

Safety: refuses to relaunch a kernel that is still running, and refuses to publish
an empty or shrinking cache — a progress dataset that went backwards would erase
completed work on the next resume.

Usage:
    python relaunch.py --shard 0            # one cycle
    python relaunch.py --shard 0 --status   # look, change nothing
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import sys
import tempfile
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(PARENT))
TID = "CHEMBL612545"
# Kaggle account that owns the shard datasets and kernels.
KAGGLE_USER = os.environ.get("KAGGLE_USER", "your-kaggle-username")


def token() -> str:
    for line in open(os.path.join(REPO, ".env")):
        if line.strip().startswith("KAGGLE_API_TOKEN="):
            return line.split("=", 1)[1].strip().strip('"').strip("'")
    raise SystemExit("KAGGLE_API_TOKEN not in .env")


def api_get(path: str, tok: str):
    req = urllib.request.Request(f"https://www.kaggle.com/api/v1/{path}",
                                 headers={"Authorization": f"Bearer {tok}"})
    return json.load(urllib.request.urlopen(req, timeout=120))


def api_post(path: str, payload: dict, tok: str):
    req = urllib.request.Request(
        f"https://www.kaggle.com/api/v1/{path}", data=json.dumps(payload).encode(),
        headers={"Authorization": f"Bearer {tok}", "Content-Type": "application/json"})
    try:
        return json.load(urllib.request.urlopen(req, timeout=180))
    except urllib.error.HTTPError as e:
        return {"_error": e.code, "_body": e.read(600).decode(errors="replace")}


def kernel_status(shard: int, tok: str) -> str:
    d = api_get(f"kernels/status?userName={KAGGLE_USER}&kernelSlug=pdl1-dock-shard-{shard}", tok)
    return d.get("status", "unknown")


def fetch_cache(shard: int, tok: str) -> dict | None:
    """Pull the shard cache out of the finished kernel's output."""
    import kagglehub  # noqa: F401  (only needed for the download path)
    d = api_get(f"kernels/output?userName={USER}&kernelSlug=pdl1-dock-shard-{shard}", tok)
    want = f"vina_cache_{TID}_shard{shard}.json"
    for f in d.get("files", []):
        if f.get("fileName", "").endswith(want) or want in f.get("fileName", ""):
            url = f.get("url")
            if not url:
                continue
            req = urllib.request.Request(url, headers={"Authorization": f"Bearer {tok}"})
            return json.loads(urllib.request.urlopen(req, timeout=300).read())
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--status", action="store_true", help="report only, change nothing")
    ap.add_argument("--force", action="store_true",
                    help="relaunch even if the kernel is still running")
    ap.add_argument("--repush-only", action="store_true",
                    help="skip fetch/publish; just restart the kernel against the "
                         "current bundle and the already-published progress dataset. "
                         "Use when the bundle changed mid-session.")
    a = ap.parse_args()
    tok = token()
    s = a.shard

    st = kernel_status(s, tok)
    local_p = os.path.join(HERE, f"vina_cache_{TID}_shard{s}.json")
    have = json.load(open(local_p)) if os.path.exists(local_p) else {}
    print(f"shard{s}: kernel status = {st}; local copy holds {len(have)} results")

    if a.status:
        return 0

    if a.repush_only:
        merged = have          # whatever is already published stands
        print(f"re-push only: reusing published progress ({len(have)} results)")
    elif st == "running" and not a.force:
        print("still running — nothing to do. Re-run when it completes.")
        return 1

    # 1. retrieve this session's results
    got = None if a.repush_only else fetch_cache(s, tok)
    if a.repush_only:
        got = {}
    if got is None:
        print("no cache in the kernel output. If the run failed, read the log first —\n"
              "relaunching over a failure just repeats it.")
        return 1
    if not a.repush_only:
        print(f"fetched {len(got)} results from the finished session")

    merged = {**have, **got}
    if len(merged) < len(have):
        print(f"REFUSING: merged cache ({len(merged)}) is smaller than the local copy "
              f"({len(have)}). Investigate before publishing — a shrinking progress "
              f"dataset would erase completed work on the next resume.")
        return 1
    if not merged:
        print("REFUSING: empty cache, nothing to publish.")
        return 1

    # 2. local copy first — it is what merge_shards.py gates on, and it survives
    #    anything that happens to the Kaggle side.
    json.dump(merged, open(local_p, "w"), indent=1)
    print(f"local copy updated: {len(merged)} results -> {local_p}")

    # 3. publish as a small progress dataset the next session can attach
    handle = f"{KAGGLE_USER}/pdl1-shard-{s}-progress"
    if a.repush_only:
        print(f"reusing existing {handle} (nothing re-published)")
    else:
        import kagglehub
        from kagglehub.config import set_kaggle_api_token
        set_kaggle_api_token(tok)
        tmp = tempfile.mkdtemp()
        json.dump(merged, open(os.path.join(tmp, f"vina_cache_{TID}_shard{s}.json"), "w"))
        kagglehub.dataset_upload(handle, tmp, version_notes=f"{len(merged)} results")
        print(f"published {handle}")

    # 4. re-push with bundle + progress attached
    src = open(os.path.join(HERE, "kernel_launcher.py")).read()
    r = api_post("kernels/push", {
        "slug": f"{KAGGLE_USER}/pdl1-dock-shard-{s}",
        "newTitle": f"pdl1-dock-shard-{s}",
        "text": src, "language": "python", "kernelType": "script",
        "isPrivate": True, "enableGpu": False, "enableTpu": False,
        "enableInternet": False,
        "datasetDataSources": [f"{KAGGLE_USER}/pdl1-docking-shard-{s}", handle],
        "competitionDataSources": [], "kernelDataSources": [],
        "modelDataSources": [], "categoryIds": [],
    }, tok)
    if r.get("_error") or r.get("error"):
        print(f"push FAILED: {json.dumps(r)[:400]}")
        return 1
    print(f"relaunched v{r.get('versionNumber')}  "
          f"https://kaggle.com/code/{KAGGLE_USER}/pdl1-dock-shard-{s}")
    man = json.load(open(os.path.join(HERE, f"bundle_{TID}", f"shard{s}", "manifest.json")))
    assigned = len(man["compounds"])
    print(f"shard{s}: {len(merged)}/{assigned} done, {assigned - len(merged)} remaining")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
