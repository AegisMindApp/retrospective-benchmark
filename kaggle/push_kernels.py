import json, os, urllib.request, urllib.error

KAGGLE_USER = os.environ.get("KAGGLE_USER", "your-kaggle-username")
tok=[l.split('=',1)[1].strip().strip('"').strip("'") for l in
     open(os.environ.get('ENV_FILE', '.env'))
     if l.strip().startswith('KAGGLE_API_TOKEN=')][0]
src_t = open(os.environ.get('KERNEL_SRC', 'kernel_launcher.py')).read()
for s in (0,1):
    body = {"slug": f"{KAGGLE_USER}/pdl1-dock-shard-{s}",
            "newTitle": f"pdl1-dock-shard-{s}",
            "text": src_t, "language": "python", "kernelType": "script",
            "isPrivate": True, "enableGpu": False, "enableTpu": False,
            "enableInternet": False,
            "datasetDataSources": [f"{KAGGLE_USER}/pdl1-docking-shard-{s}"],
            "competitionDataSources": [], "kernelDataSources": [],
            "modelDataSources": [], "categoryIds": []}
    req=urllib.request.Request("https://www.kaggle.com/api/v1/kernels/push",
        data=json.dumps(body).encode(),
        headers={"Authorization": f"Bearer {tok}", "Content-Type":"application/json"})
    try:
        r=json.load(urllib.request.urlopen(req, timeout=120))
        print(f"shard{s}: v{r.get('versionNumber')} err={r.get('error')!r} badds={r.get('invalidDatasetSources')}")
    except urllib.error.HTTPError as e:
        print(f"shard{s}: HTTP {e.code} {e.read(300).decode()[:200]}")
