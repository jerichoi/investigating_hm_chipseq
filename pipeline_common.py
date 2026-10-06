"""
Shared configuration and helper functions for the two-stage ENCODE
HOT-loci pipeline.

Stage 1 (01_plan_downloads.py):
    Queries ENCODE for TF ChIP-seq and Histone ChIP-seq experiments,
    selects the best peak file per experiment, and writes every file
    that would need to be downloaded (plus every skip and why) to a
    single manifest file: data/download_manifest.tsv
    Makes no downloads, touches no bedtools.

Stage 2 (02_download_and_process.py):
    Reads that manifest, does all the actual downloading, then runs
    the same tagging / merging / HOT-loci / blacklist-filtering logic
    the original single-script version ran inline.

Splitting the plan from the execution means the slow, API-heavy
querying step can be run once and inspected/edited as plain text, and
the slow, bandwidth-heavy downloading step can be run separately -
resumed, rate-limited, or scheduled for later - without re-querying
ENCODE.
"""

import os
import re
import time
import requests

CELL_LINES = ["A549"]
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

DOWNLOAD_MANIFEST_PATH = os.path.join(BASE_DIR, "download_manifest.tsv")
MANIFEST_FIELDS = [
    "cell_line", "kind", "name", "experiment_id",
    "output_type", "url", "dest_path", "status",
]

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


def ensure_dirs():
    for d in (PEAKS_DIR, HOT_DIR, HISTONE_PEAKS_DIR, HISTONE_MANIFEST_DIR, PYBEDTOOLS_TMPDIR):
        os.makedirs(d, exist_ok=True)


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
    """fetch experiment JSON, retrying on timeouts/connection errors/502/503/504"""
    url = ENCODE_EXPERIMENT_URL.format(experiment_id)
    for attempt in range(max_retries):
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
            r.raise_for_status()
            return r.json()
        except (requests.exceptions.Timeout,
                requests.exceptions.ConnectionError,
                requests.exceptions.HTTPError) as e:
            is_retryable_http = (
                isinstance(e, requests.exceptions.HTTPError)
                and e.response is not None
                and e.response.status_code in (502, 503, 504)
            )
            is_timeout_or_conn = isinstance(
                e, (requests.exceptions.Timeout, requests.exceptions.ConnectionError)
            )
            if (is_retryable_http or is_timeout_or_conn) and attempt < max_retries - 1:
                wait = base_delay * (2 ** attempt)
                print(f"  {type(e).__name__} on {experiment_id}, retrying in {wait}s "
                      f"(attempt {attempt + 1}/{max_retries})")
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
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with open(dest_path, "wb") as out:
        out.write(resp.content)
    return dest_path
