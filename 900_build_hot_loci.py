"""
Downloads ENCODE TF ChIP-seq peaks for a list of cell lines, builds
blacklist-filtered HOT loci for each, and separately downloads a
specified histone-mark panel per cell line.

Layout produced under BASE_DIR:
    data/
      peaks/<cell_line>/<tf_name>_<experiment_id>.bed.gz      (transient - deleted unless KEEP_RAW_PEAK_FILES)
      hot_loci/merged_tf_counts_<cell_line>.bed
      hot_loci/hot_regions_<cell_line>.bed
      hot_loci/hot_regions_<cell_line>.clean.bed              (final HOT loci, blacklist-filtered)
      histone_peaks/<cell_line>/<mark>_<experiment_id>.bed.gz (kept - these feed downstream analysis)
      histone_manifests/histone_experiments_<cell_line>.csv
      hg38-blacklist.v2.bed.gz                                (downloaded once, shared across cell lines)

Cleanup:
    - per-TF raw peak downloads are deleted immediately after being tagged
      into the merged BedTool (set KEEP_RAW_PEAK_FILES = True to retain them)
    - the intermediate decompressed .bed (from each .bed.gz) is removed as
      soon as it's been loaded into pybedtools
    - pybedtools.cleanup(remove_all=True) is called after every cell line to
      purge the temp files that merge()/sort()/intersect()/each() create
      internally, so they don't accumulate across the whole run
    - pybedtools's tempdir is pointed at a subfolder of BASE_DIR so all of
      its scratch files are easy to find and are guaranteed removed by cleanup
"""

import os
import re
import gzip
import shutil
import requests
import pybedtools
import pandas as pd
import time

# ---------------------------------------------------------------------------
# Configuration - edit this block per run
# ---------------------------------------------------------------------------
CELL_LINES = ["HepG2"]
HISTONE_MARKS = ["H3K27me3", "H3K9me3", "H3K4me3", "H3K27ac"]

ASSEMBLY = "GRCh38"
HOT_QUANTILE = 0.99
KEEP_RAW_PEAK_FILES = False  # per-TF peak downloads are only an intermediate; set True to retain them

BASE_DIR = "data"
PEAKS_DIR = os.path.join(BASE_DIR, "peaks")
HOT_DIR = os.path.join(BASE_DIR, "hot_loci")
HISTONE_PEAKS_DIR = os.path.join(BASE_DIR, "histone_peaks")
HISTONE_MANIFEST_DIR = os.path.join(BASE_DIR, "histone_manifests")
BLACKLIST_PATH = os.path.join(BASE_DIR, f"{ASSEMBLY.lower()}-blacklist.v2.bed.gz")
PYBEDTOOLS_TMPDIR = os.path.join(BASE_DIR, "tmp_pybedtools")

HEADERS = {"accept": "application/json"}
ENCODE_SEARCH_URL = "https://www.encodeproject.org/search/"
ENCODE_EXPERIMENT_URL = "https://www.encodeproject.org/experiments/{}/?format=json"
BLACKLIST_URL = "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"

PEAK_PRIORITY = [
    "optimal idr thresholded peaks",
    "IDR thresholded peaks",
    "conservative idr thresholded peaks",
    "conservative IDR thresholded peaks",
    "replicated peaks",
]

for d in (PEAKS_DIR, HOT_DIR, HISTONE_PEAKS_DIR, HISTONE_MANIFEST_DIR, PYBEDTOOLS_TMPDIR):
    os.makedirs(d, exist_ok=True)

pybedtools.set_tempdir(PYBEDTOOLS_TMPDIR)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def safe_name(name):
    """make a string filesystem-safe for use in a filename"""
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", name)


def select_peak_file(files, assembly=ASSEMBLY):
    candidates = [
        f for f in files
        if isinstance(f, dict)
        and f.get("assembly") == assembly
        and f.get("file_format") == "bed"
    ]
    for output_type in PEAK_PRIORITY:
        matches = [f for f in candidates if f.get("output_type") == output_type]
        if not matches:
            continue
        preferred = [f for f in matches if f.get("preferred_default")]
        if preferred:
            return preferred[0], output_type
        return matches[0], output_type
    return None, None


def fetch_experiments(cell_line, assay_title, extra_params=None):
    """query ENCODE for released experiments of a given assay in a given cell line"""
    params = {
        "type": "Experiment",
        "assay_title": assay_title,
        "biosample_ontology.term_name": cell_line,
        "assembly": ASSEMBLY,
        "status": "released",
        "format": "json",
        "limit": "all",
    }
    if extra_params:
        params.update(extra_params)
    r = requests.get(ENCODE_SEARCH_URL, params=params, headers=HEADERS)
    r.raise_for_status()
    return r.json()["@graph"]


def fetch_experiment_detail(experiment_id, max_retries=5, base_delay=2):
    url = f"https://www.encodeproject.org/experiments/{experiment_id}/?format=json"
    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            return r.json()
        except (requests.exceptions.Timeout,
                requests.exceptions.ConnectionError,
                requests.exceptions.HTTPError) as e:
            is_retryable_http_error = (
                isinstance(e, requests.exceptions.HTTPError)
                and e.response is not None
                and e.response.status_code in (502, 503, 504)
            )
            is_timeout_or_conn = isinstance(e, (requests.exceptions.Timeout,
                                                 requests.exceptions.ConnectionError))
            if (is_retryable_http_error or is_timeout_or_conn) and attempt < max_retries - 1:
                wait = base_delay * (2 ** attempt)
                print(f"  {type(e).__name__} on {experiment_id}, retrying in {wait}s "
                      f"(attempt {attempt+1}/{max_retries})")
                time.sleep(wait)
                continue
            raise

def get_target_name(experiment_data):
    return experiment_data.get("target", {}).get("label", "unknown")


def download_file(url, dest_path):
    if os.path.exists(dest_path):
        return dest_path
    resp = requests.get(url, headers=HEADERS)
    resp.raise_for_status()
    with open(dest_path, "wb") as out:
        out.write(resp.content)
    return dest_path


def download_blacklist():
    if not os.path.exists(BLACKLIST_PATH):
        download_file(BLACKLIST_URL, BLACKLIST_PATH)
    return BLACKLIST_PATH


def load_and_tag(fname, tag):
    """decompress if needed, tag every interval with `tag` in column 4, return a saved BedTool.
    the decompressed intermediate is deleted immediately - only the pybedtools-managed
    tagged copy (cleaned up later via pybedtools.cleanup) is needed downstream"""
    working_path = fname
    unzipped_path = None
    if fname.endswith(".gz"):
        unzipped_path = fname[:-3]
        with gzip.open(fname, "rb") as fin, open(unzipped_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        working_path = unzipped_path

    bt = pybedtools.BedTool(working_path)
    tagged = bt.each(
        lambda f: pybedtools.create_interval_from_list([f.chrom, f.start, f.end, tag])
    ).saveas()

    if unzipped_path and os.path.exists(unzipped_path):
        os.remove(unzipped_path)

    return tagged


# ---------------------------------------------------------------------------
# Step 1: build HOT loci per cell line from TF ChIP-seq
# ---------------------------------------------------------------------------
def build_hot_loci(cell_line):
    print(f"[{cell_line}] fetching TF ChIP-seq experiment list")
    experiments = fetch_experiments(cell_line, "TF ChIP-seq")
    print(f"[{cell_line}] {len(experiments)} experiments")

    cell_peaks_dir = os.path.join(PEAKS_DIR, safe_name(cell_line))
    os.makedirs(cell_peaks_dir, exist_ok=True)

    manifest = []
    skipped = []

    for experiment in experiments:
        experiment_id = experiment["accession"]
        experiment_data = fetch_experiment_detail(experiment_id)
        tf_name = safe_name(get_target_name(experiment_data))

        selected_file, output_type_used = select_peak_file(experiment_data.get("files", []))
        if selected_file is None:
            skipped.append((tf_name, experiment_id))
            continue

        file_url = "https://www.encodeproject.org" + selected_file["href"]
        fname = os.path.join(cell_peaks_dir, f"{tf_name}_{experiment_id}.bed.gz")
        manifest.append((tf_name, experiment_id, file_url, fname, output_type_used))

    print(f"[{cell_line}] {len(manifest)} usable peak files, {len(skipped)} skipped (no matching peak type)")

    all_beds = []
    for tf_name, experiment_id, file_url, fname, output_type_used in manifest:
        try:
            download_file(file_url, fname)
            all_beds.append(load_and_tag(fname, tf_name))
        except Exception as e:
            print(f"[{cell_line}] skipping {tf_name} ({experiment_id}): {e}")
        finally:
            if not KEEP_RAW_PEAK_FILES and os.path.exists(fname):
                os.remove(fname)

    if not all_beds:
        print(f"[{cell_line}] no usable TF peak files - skipping HOT loci generation")
        pybedtools.cleanup(remove_all=True)
        return None

    combined = pybedtools.BedTool.cat(*all_beds, postmerge=False).sort()
    merged = combined.merge(c=4, o="distinct,count_distinct")

    merged_path = os.path.join(HOT_DIR, f"merged_tf_counts_{safe_name(cell_line)}.bed")
    merged.saveas(merged_path)

    df = pd.read_csv(merged_path, sep="\t", names=["chrom", "start", "end", "tfs", "n_tfs"])
    print(f"[{cell_line}] n_tfs summary:\n{df['n_tfs'].describe()}")
    print(f"[{cell_line}] quantiles:\n{df['n_tfs'].quantile([0.90, 0.95, 0.99])}")

    threshold = df["n_tfs"].quantile(HOT_QUANTILE)
    hot_regions = df[df["n_tfs"] >= threshold]

    hot_path = os.path.join(HOT_DIR, f"hot_regions_{safe_name(cell_line)}.bed")
    hot_regions.to_csv(hot_path, sep="\t", index=False, header=False)

    blacklist = pybedtools.BedTool(download_blacklist())
    hot_bt = pybedtools.BedTool(hot_path)
    clean_hot = hot_bt.intersect(blacklist, v=True)

    clean_path = os.path.join(HOT_DIR, f"hot_regions_{safe_name(cell_line)}.clean.bed")
    clean_hot.saveas(clean_path)

    # purge every pybedtools-managed temp file created while processing this
    # cell line (merge/sort/intersect/each all spawn their own temp files)
    pybedtools.cleanup(remove_all=True)

    print(f"[{cell_line}] wrote {merged_path}, {hot_path}, {clean_path}")
    return clean_path


# ---------------------------------------------------------------------------
# Step 2: download a specified histone-mark panel per cell line
# ---------------------------------------------------------------------------
def download_histone_peaks(cell_line, marks):
    """
    for each requested mark, query ENCODE Histone ChIP-seq experiments in this
    cell line filtered server-side by target.label, download the best peak
    file per matching experiment, and write a manifest CSV recording what was
    found/downloaded/missing
    """
    cell_dir = os.path.join(HISTONE_PEAKS_DIR, safe_name(cell_line))
    os.makedirs(cell_dir, exist_ok=True)

    rows = []
    for mark in marks:
        print(f"[{cell_line}] fetching Histone ChIP-seq experiments for {mark}")
        experiments = fetch_experiments(
            cell_line, "Histone ChIP-seq", extra_params={"target.label": mark}
        )

        if not experiments:
            rows.append({
                "cell_line": cell_line, "mark": mark, "experiment_id": None,
                "output_type": None, "fname": None, "status": "no experiments found",
            })
            continue

        for experiment in experiments:
            experiment_id = experiment["accession"]
            experiment_data = fetch_experiment_detail(experiment_id)
            selected_file, output_type_used = select_peak_file(experiment_data.get("files", []))

            if selected_file is None:
                rows.append({
                    "cell_line": cell_line, "mark": mark, "experiment_id": experiment_id,
                    "output_type": None, "fname": None, "status": "no matching peak file",
                })
                continue

            file_url = "https://www.encodeproject.org" + selected_file["href"]
            fname = os.path.join(cell_dir, f"{safe_name(mark)}_{experiment_id}.bed.gz")

            try:
                download_file(file_url, fname)
                status = "downloaded"
            except Exception as e:
                fname = None
                status = f"download failed: {e}"

            rows.append({
                "cell_line": cell_line, "mark": mark, "experiment_id": experiment_id,
                "output_type": output_type_used, "fname": fname, "status": status,
            })

    manifest_df = pd.DataFrame(rows)
    manifest_path = os.path.join(HISTONE_MANIFEST_DIR, f"histone_experiments_{safe_name(cell_line)}.csv")
    manifest_df.to_csv(manifest_path, index=False)

    n_downloaded = (manifest_df["status"] == "downloaded").sum() if not manifest_df.empty else 0
    print(f"[{cell_line}] histone manifest -> {manifest_path} ({n_downloaded} files downloaded)")
    return manifest_df


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    hot_loci_paths = {}
    for cell_line in CELL_LINES:
        hot_loci_paths[cell_line] = build_hot_loci(cell_line)

    histone_manifests = {}
    for cell_line in CELL_LINES:
        histone_manifests[cell_line] = download_histone_peaks(cell_line, HISTONE_MARKS)

    # final safety sweep in case any step raised before its own cleanup ran
    pybedtools.cleanup(remove_all=True)

    print("\ndone.")
    for cell_line, path in hot_loci_paths.items():
        print(f"  {cell_line}: HOT loci -> {path}")
