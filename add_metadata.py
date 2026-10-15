import json
import time
from pathlib import Path

import pandas as pd
import requests

URL = "https://www.encodeproject.org/experiments/{}/?format=json&frame=embedded"

try:
    here = Path(__file__).resolve().parent
except NameError:                                   # running in a notebook
    here = Path.cwd()
manifest_dir = next((d / "data/histone_manifests" for d in (here, here.parent, Path.cwd(), Path.cwd().parent)
                     if (d / "data/histone_manifests").is_dir()), None)
if manifest_dir is None:
    raise SystemExit("data/histone_manifests not found: put this script in the project folder or r_work")
cache_fn = manifest_dir / "encode_metadata_cache.json"   # delete to refresh from ENCODE


def fetch(accession):
    r = requests.get(URL.format(accession), headers={"Accept": "application/json"}, timeout=30)
    r.raise_for_status()
    return r.json()


def describe_treatment(t):
    if not isinstance(t, dict):
        return str(t)
    amount = f"{t['amount']} {t.get('amount_units', '')}".strip() if "amount" in t else ""
    duration = f"for {t['duration']} {t.get('duration_units', '')}".strip() if "duration" in t else ""
    return " ".join(x for x in (t.get("treatment_term_name", "unknown"), amount, duration) if x)


def parse(exp):
    lab = exp.get("lab")
    lab = lab.get("title") if isinstance(lab, dict) else lab
    treatments, seen = set(), False
    for rep in exp.get("replicates", []):
        library = rep.get("library") if isinstance(rep, dict) else None
        biosample = library.get("biosample") if isinstance(library, dict) else None
        if isinstance(biosample, dict):
            seen = True
            treatments.update(describe_treatment(t) for t in biosample.get("treatments", []))
    if not seen:                                    # biosample details missing from the response
        return {"lab": lab, "treatment": "unknown"}
    return {"lab": lab, "treatment": "; ".join(sorted(treatments)) or "none"}


def main():
    cache = json.loads(cache_fn.read_text()) if cache_fn.exists() else {}
    manifests = sorted(manifest_dir.glob("histone_experiments_*.csv"))
    ids = sorted({i for m in manifests for i in pd.read_csv(m)["experiment_id"]})

    failed = []
    for acc in ids:
        if acc in cache:
            continue
        try:
            cache[acc] = parse(fetch(acc))
        except Exception as e:
            failed.append((acc, str(e)))
        time.sleep(0.2)
        cache_fn.write_text(json.dumps(cache, indent=1))

    for m in manifests:
        df = pd.read_csv(m)
        backup = m.with_name(m.name + ".bak")
        if not backup.exists():
            df.to_csv(backup, index=False)
        df = df.drop(columns=[c for c in ("lab", "treatment") if c in df.columns])
        df["lab"] = df["experiment_id"].map(lambda a: cache.get(a, {}).get("lab"))
        df["treatment"] = df["experiment_id"].map(lambda a: cache.get(a, {}).get("treatment"))
        df.to_csv(m, index=False)
        print(m.name)
        print(df.groupby(["mark", "lab", "treatment"], dropna=False).size().to_string(), "\n")

    if failed:
        print("could not fetch:", *failed, sep="\n  ")


if __name__ == "__main__":
    main()
