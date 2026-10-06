"""
Stage 2: read data/download_manifest.tsv, download every "ready" file,
then run the same tag / merge / HOT-loci / blacklist-filter logic the
original single-script pipeline ran inline.

Can be run at any later time, on a different machine, or repeatedly -
download_file() skips anything already on disk, so a partial/failed
run can just be re-invoked. It never talks to the ENCODE search API,
only the file URLs 01_plan_downloads.py already resolved.
"""

import csv
import gzip
import os
import shutil

import pandas as pd
import pybedtools

import pipeline_common as pc


def load_manifest():
    with open(pc.DOWNLOAD_MANIFEST_PATH, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def download_ready_rows(rows):
    for row in rows:
        if row["status"] != "ready":
            row["download_status"] = row["status"]
            continue
        try:
            pc.download_file(row["url"], row["dest_path"])
            row["download_status"] = "downloaded"
        except Exception as e:
            row["download_status"] = f"download_failed: {e}"
            print(f"  failed to download {row['dest_path']}: {e}")
    return rows


def load_and_tag(fname, tag):
    """decompress if needed, tag every interval with `tag` in column 4, return a saved BedTool"""
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


def build_hot_loci(cell_line, tf_rows):
    all_beds = []
    for row in tf_rows:
        if row.get("download_status") != "downloaded":
            continue
        try:
            all_beds.append(load_and_tag(row["dest_path"], row["name"]))
        except Exception as e:
            print(f"[{cell_line}] skipping {row['name']} ({row['experiment_id']}): {e}")
        finally:
            if not pc.KEEP_RAW_PEAK_FILES and os.path.exists(row["dest_path"]):
                os.remove(row["dest_path"])

    if not all_beds:
        print(f"[{cell_line}] no usable TF peak files - skipping HOT loci generation")
        pybedtools.cleanup(remove_all=True)
        return None

    combined = pybedtools.BedTool.cat(*all_beds, postmerge=False).sort()
    merged = combined.merge(c=4, o="distinct,count_distinct")

    merged_path = os.path.join(pc.HOT_DIR, f"merged_tf_counts_{pc.safe_name(cell_line)}.bed")
    merged.saveas(merged_path)

    df = pd.read_csv(merged_path, sep="\t", names=["chrom", "start", "end", "tfs", "n_tfs"])
    print(f"[{cell_line}] n_tfs summary:\n{df['n_tfs'].describe()}")
    print(f"[{cell_line}] quantiles:\n{df['n_tfs'].quantile([0.90, 0.95, 0.99])}")

    threshold = df["n_tfs"].quantile(pc.HOT_QUANTILE)
    hot_regions = df[df["n_tfs"] >= threshold]

    hot_path = os.path.join(pc.HOT_DIR, f"hot_regions_{pc.safe_name(cell_line)}.bed")
    hot_regions.to_csv(hot_path, sep="\t", index=False, header=False)

    blacklist = pybedtools.BedTool(pc.BLACKLIST_PATH)
    hot_bt = pybedtools.BedTool(hot_path)
    clean_hot = hot_bt.intersect(blacklist, v=True)

    clean_path = os.path.join(pc.HOT_DIR, f"hot_regions_{pc.safe_name(cell_line)}.clean.bed")
    clean_hot.saveas(clean_path)

    # purge every pybedtools-managed temp file created while processing this
    # cell line (merge/sort/intersect/each all spawn their own temp files)
    pybedtools.cleanup(remove_all=True)

    print(f"[{cell_line}] wrote {merged_path}, {hot_path}, {clean_path}")
    return clean_path


def write_histone_manifest(cell_line, histone_rows):
    out_rows = []
    for row in histone_rows:
        out_rows.append({
            "cell_line": row["cell_line"],
            "mark": row["name"],
            "experiment_id": row["experiment_id"] or None,
            "output_type": row["output_type"] or None,
            "fname": row["dest_path"] if row.get("download_status") == "downloaded" else None,
            "status": row.get("download_status", row["status"]),
        })

    manifest_df = pd.DataFrame(out_rows)
    manifest_path = os.path.join(pc.HISTONE_MANIFEST_DIR, f"histone_experiments_{pc.safe_name(cell_line)}.csv")
    manifest_df.to_csv(manifest_path, index=False)

    n_downloaded = (manifest_df["status"] == "downloaded").sum() if not manifest_df.empty else 0
    print(f"[{cell_line}] histone manifest -> {manifest_path} ({n_downloaded} files downloaded)")


if __name__ == "__main__":
    pc.ensure_dirs()
    pybedtools.set_tempdir(pc.PYBEDTOOLS_TMPDIR)

    rows = load_manifest()
    rows = download_ready_rows(rows)

    hot_loci_paths = {}
    for cell_line in pc.CELL_LINES:
        tf_rows = [r for r in rows if r["kind"] == "tf" and r["cell_line"] == cell_line]
        hot_loci_paths[cell_line] = build_hot_loci(cell_line, tf_rows)

    for cell_line in pc.CELL_LINES:
        histone_rows = [r for r in rows if r["kind"] == "histone" and r["cell_line"] == cell_line]
        write_histone_manifest(cell_line, histone_rows)

    pybedtools.cleanup(remove_all=True)

    print("\ndone.")
    for cell_line, path in hot_loci_paths.items():
        print(f"  {cell_line}: HOT loci -> {path}")
