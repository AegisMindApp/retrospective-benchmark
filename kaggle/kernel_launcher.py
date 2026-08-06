import glob, os, subprocess, sys, zipfile

# Do not assume a layout. The first push hardcoded
# /kaggle/input/<ds>/kaggle_shard_worker.py and died before the worker's own
# detection could run. Kaggle may present the kagglehub upload extracted, nested,
# or still as archive.zip -- handle all three, and if none match, print the tree
# so the log explains itself instead of just "No such file".
print("inputs:", sorted(glob.glob("/kaggle/input/*")), flush=True)

def find_worker():
    hits = glob.glob("/kaggle/input/**/kaggle_shard_worker.py", recursive=True)
    return hits[0] if hits else None

worker = find_worker()

if worker is None:
    zips = glob.glob("/kaggle/input/**/*.zip", recursive=True)
    print("no worker at input root; zips found:", zips, flush=True)
    if zips:
        dest = "/kaggle/working/bundle"
        os.makedirs(dest, exist_ok=True)
        with zipfile.ZipFile(zips[0]) as z:
            z.extractall(dest)
        hits = glob.glob(f"{dest}/**/kaggle_shard_worker.py", recursive=True)
        worker = hits[0] if hits else None
        print("after extract, worker:", worker, flush=True)

if worker is None:
    print("WORKER NOT FOUND. tree:", flush=True)
    for r in sorted(glob.glob("/kaggle/input/*")):
        for dp, dn, fn in os.walk(r):
            if dp.count(os.sep) > 5:
                continue
            print(" ", dp, f"({len(fn)} files)", sorted(fn)[:8], flush=True)
    sys.exit(2)

bundle = os.path.dirname(worker)
print("worker:", worker, flush=True)
print("bundle:", bundle, sorted(os.listdir(bundle))[:8], flush=True)

subprocess.run([sys.executable, worker, "--bundle", bundle, "--workers", "4"],
               check=True)
