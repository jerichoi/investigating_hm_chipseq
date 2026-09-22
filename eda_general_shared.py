"""
General / shared EDA plots for the HOT-loci false-positive ChIP-seq project.

This version is wired directly to the file structure produced by your two
existing notebooks (hot_loci.ipynb and histone_intersections.ipynb), rather
than to placeholder column names. Specifically it expects:

    genome                  'GRCh38_EBV.chrom.sizes.tsv'
    merged_tf_counts.bed    chrom, start, end, tfs (comma-list), n_tfs
                            -- output of hot_loci.ipynb's `merged` BedTool
    hot_regions.bed         pre-blacklist-filter top-percentile loci
    hot_regions.clean.bed   post-blacklist-filter (what histone_intersections
                            .ipynb calls `hot`)
    hg38-blacklist.v2.bed.gz  Boyle-Lab blacklist, already downloaded by
                            hot_loci.ipynb
    per-mark peak files     e.g. K562_H3K4me3.bed, K562_H3k27me3.bed, and
                            (not yet in your notebooks, but needed for the
                            full 4-mark comparison) K562_H3K27ac.bed and
                            K562_H3K9me3.bed -- assumed to be ENCODE
                            narrowPeak format (10 columns: chrom start end
                            name score strand signalValue pValue qValue
                            summit). If any of these are broadPeak (9 cols,
                            no summit) the parsing below degrades gracefully.
    manifest.csv            tf_name, experiment_id, output_type_used
                            -- NOT currently saved by hot_loci.ipynb. Add one
                            line at the end of that notebook's manifest-
                            building cell:
                                pd.DataFrame(
                                    manifest,
                                    columns=['tf_name','experiment_id',
                                             'file_url','fname',
                                             'output_type_used']
                                ).to_csv('manifest.csv', index=False)
                            This file is what makes plot 4 (replicate-
                            concordance proxy) possible -- see docstring
                            on that function for why.

Four plots, in the same order as before:
    1. Occupancy (n_tfs) histogram, log y-axis, with the actual 99th-
       percentile threshold your pipeline uses marked (plus 90th/95th for
       reference).
    2. Blacklist enrichment by occupancy bin, before vs after filtering.
    3. Peak width / signal (signalValue) by occupancy bin, faceted by mark.
    4. A replicate-concordance PROXY by occupancy bin: what fraction of the
       TFs contributing to a locus came from experiments that fell back to
       a lower-confidence ENCODE peak tier (no IDR-optimal set available)?
       This is a real, non-fabricated stand-in for "replicate concordance"
       given what's actually in your metadata -- see caveat in the function
       docstring on why this is a proxy and not literal IDR concordance.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
from pybedtools import BedTool

sns.set_theme(style="whitegrid", context="talk")

MARK_COLORS = {
    "H3K4me3": "#1b9e77",
    "H3K27ac": "#66c2a5",
    "H3K27me3": "#d95f02",
    "H3K9me3": "#7570b3",
}
MARK_ORDER = ["H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3"]

# Standard ENCODE narrowPeak column names (10 columns)
NARROWPEAK_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "signalValue", "pValue", "qValue", "summit",
]
# broadPeak has 9 columns (no summit)
BROADPEAK_COLS = NARROWPEAK_COLS[:-1]


# --------------------------------------------------------------------
# Occupancy bin assignment (shared helper)
# --------------------------------------------------------------------
def assign_occupancy_bins(df, count_col="n_tfs",
                           edges=(0, 50, 80, 90, 95, 99, 100)):
    """
    Bins loci by percentile of `count_col`, using value-based quantile
    cutoffs (matching the .quantile(0.99) logic already in hot_loci.ipynb)
    rather than rank-based cuts. Ties at bin edges are common given how
    concentrated low occupancy counts are, so duplicates='drop' merges any
    collapsed edges -- check `df['occupancy_bin'].value_counts()` after
    running this to confirm you got the resolution you expected, since a
    highly tied distribution can silently collapse e.g. the 90-95 bin.
    """
    quantiles = [e / 100 for e in edges]
    bin_edges = df[count_col].quantile(quantiles).values
    labels = [f"{edges[i]}-{edges[i+1]}" for i in range(len(edges) - 1)]
    df = df.copy()
    df["occupancy_bin"] = pd.cut(
        df[count_col], bins=bin_edges, labels=labels,
        include_lowest=True, duplicates="drop"
    )
    return df


# --------------------------------------------------------------------
# 1. Occupancy histogram
# --------------------------------------------------------------------
def plot_occupancy_histogram(merged_bed_path="merged_tf_counts.bed",
                              save_path="plot1_occupancy_histogram.png"):
    df = pd.read_csv(merged_bed_path, sep="\t",
                      names=["chrom", "start", "end", "tfs", "n_tfs"])

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.hist(df["n_tfs"], bins=100, color="#4c72b0", edgecolor="none")
    ax.set_yscale("log")
    ax.set_xlabel("Number of distinct TFs bound (n_tfs)")
    ax.set_ylabel("Number of loci (log scale)")
    ax.set_title("Genome-wide distribution of TF occupancy counts")

    # 99th percentile marked to match the actual threshold used in
    # hot_loci.ipynb to define hot_regions.bed; 90th/95th shown for context.
    for p, style in [(90, ":"), (95, "--"), (99, "-")]:
        cutoff = df["n_tfs"].quantile(p / 100)
        ax.axvline(cutoff, color="firebrick", linestyle=style, linewidth=1.3)
        ax.text(cutoff, ax.get_ylim()[1] * 0.7, f"{p}th pct",
                rotation=90, va="top", ha="right", fontsize=10, color="firebrick")

    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path}")
    return df


# --------------------------------------------------------------------
# 2. Blacklist enrichment by occupancy bin
# --------------------------------------------------------------------
def plot_blacklist_enrichment(merged_bed_path="merged_tf_counts.bed",
                               blacklist_path="hg38-blacklist.v2.bed.gz",
                               genome="GRCh38_EBV.chrom.sizes.tsv",
                               save_path="plot2_blacklist_enrichment.png"):
    """
    Bins ALL loci in merged_tf_counts.bed (not just the top-percentile hot
    regions) by occupancy, then checks blacklist overlap per bin. This is
    what actually justifies filtering specifically at the top bin(s) --
    the pipeline currently only intersects hot_regions.bed against the
    blacklist, so this plot answers a slightly broader question ("is
    blacklist contamination specific to high occupancy, or spread evenly")
    that your current notebooks don't directly show.
    """
    df = pd.read_csv(merged_bed_path, sep="\t",
                      names=["chrom", "start", "end", "tfs", "n_tfs"])
    df = assign_occupancy_bins(df)
    df = df.reset_index().rename(columns={"index": "locus_id"})

    loci_bt = BedTool.from_dataframe(
        df[["chrom", "start", "end", "locus_id"]]
    ).sort(g=genome)
    blacklist_bt = BedTool(blacklist_path).sort(g=genome)

    # -c appends a count of overlapping blacklist intervals per locus
    tagged = loci_bt.intersect(blacklist_bt, c=True)
    tagged_df = tagged.to_dataframe(
        names=["chrom", "start", "end", "locus_id", "blacklist_count"]
    )
    df = df.merge(tagged_df[["locus_id", "blacklist_count"]], on="locus_id")
    df["in_blacklist"] = df["blacklist_count"] > 0

    before = (
        df.groupby("occupancy_bin", observed=True)["in_blacklist"]
        .mean().mul(100).rename("pct_blacklisted").reset_index()
    )
    after_df = df[~df["in_blacklist"]]
    after = (
        after_df.groupby("occupancy_bin", observed=True)["in_blacklist"]
        .mean().mul(100).rename("pct_blacklisted").reset_index()
    )

    bin_order = list(df["occupancy_bin"].cat.categories)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    sns.barplot(data=before, x="occupancy_bin", y="pct_blacklisted",
                order=bin_order, ax=axes[0], color="#c44e52")
    axes[0].set_title("Before blacklist filtering")
    axes[0].set_ylabel("% of loci overlapping blacklist")
    axes[0].set_xlabel("Occupancy percentile bin")
    axes[0].tick_params(axis="x", rotation=45)

    sns.barplot(data=after, x="occupancy_bin", y="pct_blacklisted",
                order=bin_order, ax=axes[1], color="#55a868")
    axes[1].set_title("After blacklist filtering (sanity check, should be ~0)")
    axes[1].set_xlabel("Occupancy percentile bin")
    axes[1].tick_params(axis="x", rotation=45)

    fig.suptitle("Blacklist contamination across occupancy bins")
    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path}")
    return df


# --------------------------------------------------------------------
# 3. Peak width / signal by occupancy bin, faceted by mark
# --------------------------------------------------------------------
def _read_peak_bed(path):
    """Reads a mark's peak file, tolerating narrowPeak (10 col) or
    broadPeak (9 col, no summit) formats."""
    n_cols = pd.read_csv(path, sep="\t", nrows=1, header=None).shape[1]
    if n_cols >= 10:
        names = NARROWPEAK_COLS
    elif n_cols == 9:
        names = BROADPEAK_COLS
    else:
        # Minimal fallback: at least chrom/start/end, no signal available
        names = ["chrom", "start", "end"] + [f"col{i}" for i in range(3, n_cols)]
    return names


def plot_peak_quality_by_bin(binned_loci_df, mark_bed_paths,
                              genome="GRCh38_EBV.chrom.sizes.tsv",
                              save_path_prefix="plot3_peak_quality"):
    """
    Parameters
    ----------
    binned_loci_df : DataFrame with chrom, start, end, occupancy_bin
        (e.g. the return value of plot_blacklist_enrichment(), or run
        assign_occupancy_bins() on merged_tf_counts.bed yourself)
    mark_bed_paths : dict, e.g.
        {'H3K4me3': 'K562_H3K4me3.bed', 'H3K27ac': 'K562_H3K27ac.bed',
         'H3K27me3': 'K562_H3k27me3.bed', 'H3K9me3': 'K562_H3K9me3.bed'}
        Only marks present in this dict are plotted -- if you only have
        two mark files right now, pass just those two and the figure will
        show fewer panels rather than erroring.
    """
    # Preserve bin order explicitly -- casting the categorical to str for the
    # BedTool round-trip loses ordering, which otherwise leaves the x-axis
    # scrambled (e.g. "80-90" before "50-80") since seaborn then falls back
    # to order-of-first-appearance in the intersected data.
    bin_order = [b for b in binned_loci_df["occupancy_bin"].cat.categories
                 if b in binned_loci_df["occupancy_bin"].unique()]

    bins_bt = BedTool.from_dataframe(
        binned_loci_df[["chrom", "start", "end", "occupancy_bin"]].astype(
            {"occupancy_bin": str}
        )
    ).sort(g=genome)

    all_hits = []
    for mark, path in mark_bed_paths.items():
        col_names = _read_peak_bed(path)
        mark_bt = BedTool(path).sort(g=genome)
        hits = bins_bt.intersect(mark_bt, wa=True, wb=True)
        hit_names = ["chrom", "start", "end", "occupancy_bin"] + \
                    [f"m_{c}" for c in col_names]
        hits_df = hits.to_dataframe(names=hit_names)
        hits_df["peak_width"] = hits_df["m_end"] - hits_df["m_start"]
        hits_df["mark"] = mark
        # signalValue if available, else fall back to score column
        if "m_signalValue" in hits_df.columns:
            hits_df["signal"] = hits_df["m_signalValue"]
        elif "m_score" in hits_df.columns:
            hits_df["signal"] = hits_df["m_score"]
        else:
            hits_df["signal"] = np.nan
        all_hits.append(hits_df[["occupancy_bin", "mark", "peak_width", "signal"]])

    peaks_df = pd.concat(all_hits, ignore_index=True)

    present_marks = [m for m in MARK_ORDER if m in mark_bed_paths]
    for metric, ylabel, suffix in [
        ("peak_width", "Peak width (bp)", "width"),
        ("signal", "Signal (signalValue or score)", "signal"),
    ]:
        fig, axes = plt.subplots(1, len(present_marks),
                                  figsize=(5 * len(present_marks), 5.5),
                                  sharey=False, squeeze=False)
        for ax, mark in zip(axes[0], present_marks):
            sub = peaks_df[peaks_df["mark"] == mark]
            sns.boxplot(data=sub, x="occupancy_bin", y=metric, ax=ax,
                        order=bin_order, color=MARK_COLORS[mark], showfliers=False)
            ax.set_title(mark)
            ax.set_xlabel("Occupancy bin")
            ax.set_ylabel(ylabel if mark == present_marks[0] else "")
            ax.tick_params(axis="x", rotation=45)

        fig.suptitle(f"{ylabel} across occupancy bins, by mark")
        fig.tight_layout()
        out_path = f"{save_path_prefix}_{suffix}.png"
        fig.savefig(out_path, dpi=300)
        plt.close(fig)
        print(f"Saved: {out_path}")

    return peaks_df


# --------------------------------------------------------------------
# 4. Replicate-concordance PROXY by occupancy bin
# --------------------------------------------------------------------
def plot_concordance_proxy_by_bin(binned_loci_df, manifest_csv="manifest.csv",
                                   save_path="plot4_concordance_proxy.png"):
    """
    CAVEAT -- read before using: hot_loci.ipynb's `select_peak_file()`
    picks ONE peak set per experiment from a priority list (optimal IDR >
    IDR > conservative IDR > replicated peaks). It doesn't retain a
    numeric IDR score, so a literal "% peaks passing an IDR threshold"
    plot (as originally sketched) isn't available from what your pipeline
    currently saves.

    What IS available: which tier of file ENCODE provided for each TF.
    Falling back to 'replicated peaks' (no IDR set at all) or a
    'conservative' IDR set is itself a signal of weaker cross-replicate
    agreement for that TF's ChIP-seq. This function uses that as a proxy:
    for each locus, what fraction of its contributing TFs come from
    experiments that did NOT have an 'optimal idr thresholded peaks' file
    available? Then compares that fraction across occupancy bins.

    This is a coarser, TF-level proxy, not a peak-level IDR reproducibility
    metric -- worth flagging as a limitation if used in the actual proposal
    text, and worth revisiting once/if per-replicate IDR scores are pulled
    from ENCODE directly (they're available per-experiment as separate
    files, just not currently downloaded by hot_loci.ipynb).

    Requires binned_loci_df to have a 'tfs' column (comma-separated TF
    names), which merged_tf_counts.bed provides but hot_regions.bed (after
    the .quantile() filter) does not retain unless you keep that column
    through the filtering step.
    """
    manifest_df = pd.read_csv(manifest_csv)
    tier_by_tf = manifest_df.groupby("tf_name")["output_type_used"].first()
    low_conc_tfs = set(
        tier_by_tf[tier_by_tf != "optimal idr thresholded peaks"].index
    )

    def frac_low_concordance(tfs_str):
        if pd.isna(tfs_str) or tfs_str == "":
            return np.nan
        tfs = tfs_str.split(",")
        return sum(t in low_conc_tfs for t in tfs) / len(tfs)

    df = binned_loci_df.copy()
    df["frac_low_concordance_tfs"] = df["tfs"].apply(frac_low_concordance)

    summary = (
        df.groupby("occupancy_bin", observed=True)["frac_low_concordance_tfs"]
        .mean().reset_index()
    )

    bin_order = list(df["occupancy_bin"].cat.categories)

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.barplot(data=summary, x="occupancy_bin", y="frac_low_concordance_tfs",
                order=bin_order, ax=ax, color="#8172b2")
    ax.set_xlabel("Occupancy percentile bin")
    ax.set_ylabel("Mean fraction of contributing TFs\nfrom non-optimal-IDR experiments")
    ax.set_title("Replicate-concordance proxy across occupancy bins")
    ax.tick_params(axis="x", rotation=45)

    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path}")
    return df


# --------------------------------------------------------------------
# Main -- adjust paths to match your working directory
# --------------------------------------------------------------------
if __name__ == "__main__":
    GENOME = "GRCh38_EBV.chrom.sizes.tsv"
    MERGED_BED = "merged_tf_counts.bed"
    BLACKLIST = "hg38-blacklist.v2.bed.gz"
    MARK_BEDS = {
        "H3K4me3": "K562_H3K4me3.bed",
        "H3K27me3": "K562_H3k27me3.bed",
        # Add these once you've pulled them (not yet in your notebooks):
        # "H3K27ac": "K562_H3K27ac.bed",
        # "H3K9me3": "K562_H3K9me3.bed",
    }
    MANIFEST_CSV = "manifest.csv"

    plot_occupancy_histogram(MERGED_BED)
    binned_df = plot_blacklist_enrichment(MERGED_BED, BLACKLIST, GENOME)
    plot_peak_quality_by_bin(binned_df, MARK_BEDS, GENOME)

    # Only runs if you've added the manifest.csv-saving line to
    # hot_loci.ipynb as described in the module docstring above.
    try:
        plot_concordance_proxy_by_bin(binned_df, MANIFEST_CSV)
    except FileNotFoundError:
        print(f"Skipping plot 4: {MANIFEST_CSV} not found. "
              "See module docstring for the one-line addition needed "
              "in hot_loci.ipynb to produce it.")
