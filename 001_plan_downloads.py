"""
Stage 1: query ENCODE and write data/download_manifest.tsv

Does NOT download any peak/blacklist files and does NOT import
pybedtools - it only hits the ENCODE search/experiment JSON endpoints,
decides which peak file (if any) each experiment would contribute, and
records the outcome as one row per experiment (or per mark, if no
experiments were found at all, or per blacklist file).

Run this whenever you want to (re)plan a pull - e.g. on a machine with
API access but no bandwidth/time to download right now. The manifest
it writes is a plain TSV you can open, diff, or hand-edit before
02_download_and_process.py runs, and can be handed off to run at a
different time or on a different machine.
"""

import csv
import os

import pipeline_common as pc


def plan_tf_rows(cell_line):
    print(f"[{cell_line}] fetching TF ChIP-seq experiment list")
    experiments = pc.fetch_experiments(cell_line, "TF ChIP-seq")
    print(f"[{cell_line}] {len(experiments)} experiments")

    cell_peaks_dir = os.path.join(pc.PEAKS_DIR, pc.safe_name(cell_line))
    rows = []
    for experiment in experiments:
        experiment_id = experiment["accession"]
        experiment_data = pc.fetch_experiment_detail(experiment_id)
        tf_name = pc.safe_name(pc.get_target_name(experiment_data))

        selected_file, output_type = pc.select_peak_file(experiment_data.get("files", []))
        if selected_file is None:
            rows.append({
                "cell_line": cell_line, "kind": "tf", "name": tf_name,
                "experiment_id": experiment_id, "output_type": "",
                "url": "", "dest_path": "", "status": "skipped_no_matching_peak_type",
            })
            continue

        url = "https://www.encodeproject.org" + selected_file["href"]
        dest_path = os.path.join(cell_peaks_dir, f"{tf_name}_{experiment_id}.bed.gz")
        rows.append({
            "cell_line": cell_line, "kind": "tf", "name": tf_name,
            "experiment_id": experiment_id, "output_type": output_type,
            "url": url, "dest_path": dest_path, "status": "ready",
        })

    n_ready = sum(1 for r in rows if r["status"] == "ready")
    print(f"[{cell_line}] {n_ready} usable peak files, {len(rows) - n_ready} skipped (no matching peak type)")
    return rows


def plan_histone_rows(cell_line, marks):
    cell_dir = os.path.join(pc.HISTONE_PEAKS_DIR, pc.safe_name(cell_line))
    rows = []
    for mark in marks:
        print(f"[{cell_line}] fetching Histone ChIP-seq experiments for {mark}")
        experiments = pc.fetch_experiments(
            cell_line, "Histone ChIP-seq", extra_params={"target.label": mark}
        )

        if not experiments:
            rows.append({
                "cell_line": cell_line, "kind": "histone", "name": mark,
                "experiment_id": "", "output_type": "", "url": "", "dest_path": "",
                "status": "no_experiments_found",
            })
            continue

        for experiment in experiments:
            experiment_id = experiment["accession"]
            experiment_data = pc.fetch_experiment_detail(experiment_id)
            selected_file, output_type = pc.select_peak_file(experiment_data.get("files", []))

            if selected_file is None:
                rows.append({
                    "cell_line": cell_line, "kind": "histone", "name": mark,
                    "experiment_id": experiment_id, "output_type": "",
                    "url": "", "dest_path": "", "status": "no_matching_peak_file",
                })
                continue

            url = "https://www.encodeproject.org" + selected_file["href"]
            dest_path = os.path.join(cell_dir, f"{pc.safe_name(mark)}_{experiment_id}.bed.gz")
            rows.append({
                "cell_line": cell_line, "kind": "histone", "name": mark,
                "experiment_id": experiment_id, "output_type": output_type,
                "url": url, "dest_path": dest_path, "status": "ready",
            })

    return rows


def plan_blacklist_row():
    return {
        "cell_line": "", "kind": "blacklist", "name": pc.ASSEMBLY,
        "experiment_id": "", "output_type": "",
        "url": pc.BLACKLIST_URL, "dest_path": pc.BLACKLIST_PATH, "status": "ready",
    }


if __name__ == "__main__":
    pc.ensure_dirs()

    all_rows = [plan_blacklist_row()]
    for cell_line in pc.CELL_LINES:
        all_rows.extend(plan_tf_rows(cell_line))
    for cell_line in pc.CELL_LINES:
        all_rows.extend(plan_histone_rows(cell_line, pc.HISTONE_MARKS))

    with open(pc.DOWNLOAD_MANIFEST_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=pc.MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    n_ready = sum(1 for r in all_rows if r["status"] == "ready")
    print(f"\nwrote {pc.DOWNLOAD_MANIFEST_PATH}: {len(all_rows)} rows, {n_ready} ready to download")
