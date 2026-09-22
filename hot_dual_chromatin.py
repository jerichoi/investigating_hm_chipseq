#!/usr/bin/env python3
"""
hot_dual_chromatin.py

Given a HOT loci BED file for a cell line (e.g. K562) and a set of
histone-mark ChIP-seq peak files, flag HOT loci that carry BOTH active
(euchromatic) and repressive (heterochromatic) marks -- and score them
by how likely that combination is to represent a false-positive TF
binding call, rather than genuine biology.

Design rationale:
  - active + CONSTITUTIVE heterochromatin (H3K9me3, H4K20me3)  -> mechanistically
    implausible for a real TF to bind; strong false-positive evidence.
  - active + FACULTATIVE heterochromatin (H3K27me3, Polycomb)  -> can reflect
    genuine bivalent/poised chromatin; weaker, more ambiguous evidence.
  - optionally, also overlapping the ChIP input/IgG control track is
    independent, direct evidence of nonspecific pulldown (strongest signal,
    only used if you supply a "control" mark).

Uses pybedtools exclusively for interval work. No pandas.

-------------------------------------------------------------------
1. DOWNLOAD (from ENCODE, K562, GRCh38, "optimal IDR thresholded peaks")
-------------------------------------------------------------------
  eu:              H3K4me3, H3K27ac  (H3K4me1 optional)
  het_facultative: H3K27me3
  het_constitutive: H3K9me3  (H4K20me3 optional)
  control (optional): matched ChIP input / IgG peaks, if available

  Blacklist (fixed file, not per cell line):
  https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz

-------------------------------------------------------------------
2. CONFIG FILE (JSON) -- tell the script which file is which mark/class
-------------------------------------------------------------------
{
  "H3K4me3":  {"path": "downloads/H3K4me3.bed.gz",  "class": "eu"},
  "H3K27ac":  {"path": "downloads/H3K27ac.bed.gz",  "class": "eu"},
  "H3K27me3": {"path": "downloads/H3K27me3.bed.gz", "class": "het_facultative"},
  "H3K9me3":  {"path": "downloads/H3K9me3.bed.gz",  "class": "het_constitutive"}
}
(class must be one of: eu, het_facultative, het_constitutive, control)
Peak files may be gzipped (.gz) or plain BED/narrowPeak/broadPeak -- pybedtools
handles gzip transparently.

-------------------------------------------------------------------
3. RUN
-------------------------------------------------------------------
python hot_dual_chromatin.py \
    --hot hot_regions.clean.bed \
    --config config.json \
    --outdir results/ \
    --blacklist hg38-blacklist.v2.bed.gz

Add --shuffle-control --genome hg38.chrom.sizes to test whether the
observed rate of dual-chromatin loci exceeds a random-placement baseline
(chrom.sizes: e.g. from UCSC, two columns "chrom<TAB>size").

-------------------------------------------------------------------
OUTPUT (in --outdir)
-------------------------------------------------------------------
hot_marks_matrix.tsv       per-locus overlap bp per mark + flags + fp_score
fp_candidates_ranked.bed   all loci with fp_score > 0, sorted high -> low
dual_constitutive.bed      strongest false-positive candidates
dual_facultative_only.bed  ambiguous -- possible genuine bivalency
eu_only.bed / het_only.bed / neither.bed
shuffle_control_fp_fractions.txt  (only with --shuffle-control)
"""

import argparse
import json
import os
import sys
from collections import defaultdict

import pybedtools
from pybedtools import BedTool

VALID_CLASSES = ("eu", "het_facultative", "het_constitutive", "control")

# Scoring weights -- transparent, easy to justify/tune in a methods section.
SCORE_EU_PLUS_CONSTITUTIVE = 2.0   # eu + constitutive het: strong FP evidence
SCORE_EU_PLUS_FACULTATIVE = 1.0    # eu + facultative het only: weaker/ambiguous
SCORE_CONTROL_ENRICHED = 2.0       # also enriched in input/IgG: independent, strong evidence


def locus_key(interval):
    return f"{interval.chrom}:{interval.start}-{interval.end}"


def load_config(config_path):
    with open(config_path) as fh:
        cfg = json.load(fh)
    for mark, info in cfg.items():
        if "path" not in info or "class" not in info:
            sys.exit(f"[config error] mark '{mark}' needs 'path' and 'class'")
        if info["class"] not in VALID_CLASSES:
            sys.exit(f"[config error] mark '{mark}' class must be one of {VALID_CLASSES}")
        if not os.path.exists(info["path"]):
            sys.exit(f"[config error] file not found for mark '{mark}': {info['path']}")
    return cfg


def clean_hot_loci(hot_bed, blacklist_bed=None, mappability_exclude_bed=None):
    """Load HOT loci; optionally drop ENCODE blacklist and/or low-mappability overlaps."""
    hot = BedTool(hot_bed).sort()

    if blacklist_bed:
        bl = BedTool(blacklist_bed).sort()
        n_before = hot.count()
        hot = hot.intersect(bl, v=True)
        print(f"[blacklist] removed {n_before - hot.count()} / {n_before} HOT loci "
              f"overlapping blacklist regions", file=sys.stderr)

    if mappability_exclude_bed:
        mb = BedTool(mappability_exclude_bed).sort()
        n_before = hot.count()
        hot = hot.sort().intersect(mb, v=True)
        print(f"[mappability] removed {n_before - hot.count()} / {n_before} HOT loci "
              f"overlapping low-mappability regions", file=sys.stderr)

    return hot.sort()


def overlap_counts(hot_bt, mark_bt, min_frac=0.0):
    """dict: locus_key -> total overlap bp between HOT locus and mark peaks (bedtools -wo)."""
    counts = defaultdict(int)
    kwargs = dict(wo=True)
    if min_frac > 0:
        kwargs["f"] = min_frac
    for row in hot_bt.intersect(mark_bt, **kwargs):
        chrom, start, end = row.fields[0], row.fields[1], row.fields[2]
        key = f"{chrom}:{start}-{end}"
        counts[key] += int(row.fields[-1])
    return counts


def classify_locus(has_eu, has_fac, has_con, has_control):
    """
    Returns (category_label, fp_score).
      dual_constitutive      eu + constitutive het present (strongest FP candidates)
      dual_facultative_only  eu + facultative het present, NO constitutive het (ambiguous)
      eu_only                eu present, no het at all
      het_only                het present (either kind), no eu
      neither                 no eu, no het overlap
    """
    has_het_any = has_fac or has_con
    score = 0.0

    if has_eu and has_con:
        score += SCORE_EU_PLUS_CONSTITUTIVE
        category = "dual_constitutive"
    elif has_eu and has_fac:
        score += SCORE_EU_PLUS_FACULTATIVE
        category = "dual_facultative_only"
    elif has_eu:
        category = "eu_only"
    elif has_het_any:
        category = "het_only"
    else:
        category = "neither"

    if has_control:
        score += SCORE_CONTROL_ENRICHED

    return category, score


def write_matrix(hot_bt, mark_names, mark_of_class, per_mark_counts, has_control_class, outpath):
    """
    Write TSV with per-mark overlap bp + eu/fac/con/control flags + category + fp_score.
    Returns: dict locus_key -> (category, fp_score), dict locus_key -> interval
    """
    results = {}
    interval_lookup = {}

    header = ["chrom", "start", "end"] + mark_names + [
        "has_eu", "has_het_facultative", "has_het_constitutive", "has_control",
        "category", "fp_score",
    ]

    with open(outpath, "w") as out:
        out.write("\t".join(header) + "\n")

        for iv in hot_bt:
            key = locus_key(iv)
            interval_lookup[key] = iv

            row = [iv.chrom, str(iv.start), str(iv.end)]
            has_eu = has_fac = has_con = has_ctrl = False

            for mark in mark_names:
                bp = per_mark_counts[mark].get(key, 0)
                row.append(str(bp))
                if bp > 0:
                    cls = mark_of_class[mark]
                    if cls == "eu":
                        has_eu = True
                    elif cls == "het_facultative":
                        has_fac = True
                    elif cls == "het_constitutive":
                        has_con = True
                    elif cls == "control":
                        has_ctrl = True

            category, fp_score = classify_locus(has_eu, has_fac, has_con, has_ctrl)
            row += [str(int(has_eu)), str(int(has_fac)), str(int(has_con)),
                    str(int(has_ctrl)) if has_control_class else "NA",
                    category, f"{fp_score:.1f}"]
            out.write("\t".join(row) + "\n")

            results[key] = (category, fp_score)

    return results, interval_lookup


def write_bed_subset(results, interval_lookup, target_category, outpath):
    n = 0
    with open(outpath, "w") as out:
        for key, (category, score) in results.items():
            if category == target_category:
                iv = interval_lookup[key]
                out.write(f"{iv.chrom}\t{iv.start}\t{iv.end}\t.\t{score:.1f}\n")
                n += 1
    return n


def write_ranked_candidates(results, interval_lookup, outpath):
    """All loci with fp_score > 0, sorted descending -- prioritized false-positive
    candidate list for manual review / orthogonal follow-up."""
    scored = [(key, cat, score) for key, (cat, score) in results.items() if score > 0]
    scored.sort(key=lambda x: x[2], reverse=True)

    with open(outpath, "w") as out:
        for key, category, score in scored:
            iv = interval_lookup[key]
            out.write(f"{iv.chrom}\t{iv.start}\t{iv.end}\t{category}\t{score:.1f}\n")

    return len(scored)


def shuffle_control(hot_bt, genome_file, mark_bts, mark_of_class, mark_names,
                     n_iter, seed, outdir, blacklist_bed=None):
    """Randomize HOT loci genome-wide and recompute fraction with fp_score > 0,
    to test whether observed dual-chromatin enrichment exceeds chance."""
    excl_kwargs = {}
    if blacklist_bed:
        excl_kwargs["excl"] = blacklist_bed

    eu_marks = [m for m in mark_names if mark_of_class[m] == "eu"]
    fac_marks = [m for m in mark_names if mark_of_class[m] == "het_facultative"]
    con_marks = [m for m in mark_names if mark_of_class[m] == "het_constitutive"]

    eu_union = BedTool.cat(*[mark_bts[m] for m in eu_marks], postmerge=True) if eu_marks else None
    fac_union = BedTool.cat(*[mark_bts[m] for m in fac_marks], postmerge=True) if fac_marks else None
    con_union = BedTool.cat(*[mark_bts[m] for m in con_marks], postmerge=True) if con_marks else None

    fp_fracs = []
    for i in range(n_iter):
        shuffled = hot_bt.shuffle(g=genome_file, chrom=True, seed=seed + i, **excl_kwargs).sort()
        n_total = shuffled.count()
        if n_total == 0:
            fp_fracs.append(0.0)
            continue

        eu_keys = {locus_key(iv) for iv in shuffled.intersect(eu_union, u=True)} if eu_union else set()
        fac_keys = {locus_key(iv) for iv in shuffled.intersect(fac_union, u=True)} if fac_union else set()
        con_keys = {locus_key(iv) for iv in shuffled.intersect(con_union, u=True)} if con_union else set()

        n_fp = len(eu_keys & (fac_keys | con_keys))
        fp_fracs.append(n_fp / n_total)

    obs_path = os.path.join(outdir, "shuffle_control_fp_fractions.txt")
    with open(obs_path, "w") as out:
        out.write("iteration\tfp_fraction\n")
        for i, frac in enumerate(fp_fracs):
            out.write(f"{i}\t{frac:.6f}\n")

    mean_frac = sum(fp_fracs) / len(fp_fracs) if fp_fracs else 0.0
    print(f"[shuffle control] {n_iter} iterations, mean fp-candidate fraction "
          f"under random placement = {mean_frac:.4f} (see {obs_path})", file=sys.stderr)
    return mean_frac


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hot", required=True, help="HOT loci BED file, e.g. hot_regions.clean.bed")
    ap.add_argument("--config", required=True, help="JSON config: mark -> {path, class}")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--blacklist", default=None, help="ENCODE hg38 blacklist BED(.gz)")
    ap.add_argument("--mappability-exclude", default=None,
                     help="BED of low-mappability regions to exclude (optional)")
    ap.add_argument("--min-frac", type=float, default=0.0,
                     help="Min fraction of HOT locus a peak must cover to count (bedtools -f). Default: any overlap.")
    ap.add_argument("--shuffle-control", action="store_true")
    ap.add_argument("--genome", default=None, help="chrom sizes file, required for --shuffle-control")
    ap.add_argument("--n-shuffle", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if args.shuffle_control and not args.genome:
        sys.exit("--shuffle-control requires --genome (chrom sizes file)")

    os.makedirs(args.outdir, exist_ok=True)
    pybedtools.set_tempdir(args.outdir)

    cfg = load_config(args.config)
    mark_names = list(cfg.keys())
    mark_of_class = {m: cfg[m]["class"] for m in mark_names}
    has_control_class = any(c == "control" for c in mark_of_class.values())

    n_eu = sum(1 for c in mark_of_class.values() if c == "eu")
    n_fac = sum(1 for c in mark_of_class.values() if c == "het_facultative")
    n_con = sum(1 for c in mark_of_class.values() if c == "het_constitutive")
    n_ctrl = sum(1 for c in mark_of_class.values() if c == "control")
    print(f"[config] {len(mark_names)} marks: {n_eu} eu / {n_fac} het_facultative / "
          f"{n_con} het_constitutive / {n_ctrl} control", file=sys.stderr)
    if n_eu == 0 or (n_fac == 0 and n_con == 0):
        sys.exit("[config error] need >=1 'eu' mark and >=1 'het_facultative' or "
                  "'het_constitutive' mark to evaluate dual-chromatin false positives")

    print("[step 1] loading + cleaning HOT loci", file=sys.stderr)
    hot_bt = clean_hot_loci(args.hot, blacklist_bed=args.blacklist,
                             mappability_exclude_bed=args.mappability_exclude)
    n_hot = hot_bt.count()
    print(f"[step 1] {n_hot} HOT loci retained", file=sys.stderr)

    print("[step 2] loading mark peak files + intersecting", file=sys.stderr)
    mark_bts, per_mark_counts = {}, {}
    for mark in mark_names:
        mbt = BedTool(cfg[mark]["path"]).sort()
        mark_bts[mark] = mbt
        per_mark_counts[mark] = overlap_counts(hot_bt, mbt, min_frac=args.min_frac)
        n_loci_hit = len(per_mark_counts[mark])
        print(f"    {mark:12s} ({mark_of_class[mark]:16s}): {n_loci_hit} / {n_hot} HOT loci overlap",
              file=sys.stderr)

    print("[step 3] building matrix + classifying loci", file=sys.stderr)
    matrix_path = os.path.join(args.outdir, "hot_marks_matrix.tsv")
    results, interval_lookup = write_matrix(hot_bt, mark_names, mark_of_class,
                                             per_mark_counts, has_control_class, matrix_path)

    counts_by_category = defaultdict(int)
    for category, _ in results.values():
        counts_by_category[category] += 1

    print("[step 4] writing BED subsets + ranked FP candidate list", file=sys.stderr)
    for category, fname in [("dual_constitutive", "dual_constitutive.bed"),
                             ("dual_facultative_only", "dual_facultative_only.bed"),
                             ("eu_only", "eu_only.bed"),
                             ("het_only", "het_only.bed"),
                             ("neither", "neither.bed")]:
        outpath = os.path.join(args.outdir, fname)
        n = write_bed_subset(results, interval_lookup, category, outpath)
        print(f"    {fname}: {n} loci", file=sys.stderr)

    n_ranked = write_ranked_candidates(results, interval_lookup,
                                        os.path.join(args.outdir, "fp_candidates_ranked.bed"))
    print(f"    fp_candidates_ranked.bed: {n_ranked} loci (fp_score > 0, sorted high->low)",
          file=sys.stderr)

    print("\n[summary]", file=sys.stderr)
    print(f"  total HOT loci analyzed : {n_hot}", file=sys.stderr)
    for category in ("dual_constitutive", "dual_facultative_only", "eu_only", "het_only", "neither"):
        n = counts_by_category.get(category, 0)
        pct = 100 * n / n_hot if n_hot else 0
        print(f"  {category:22s}: {n:6d} ({pct:5.1f}%)", file=sys.stderr)
    print("\n  Interpretation guide:", file=sys.stderr)
    print("  - dual_constitutive      : strongest false-positive candidates", file=sys.stderr)
    print("  - dual_facultative_only  : ambiguous -- may be genuine bivalent binding", file=sys.stderr)

    if args.shuffle_control:
        print("\n[step 5] running shuffle control", file=sys.stderr)
        n_fp_obs = counts_by_category.get("dual_constitutive", 0) + \
            counts_by_category.get("dual_facultative_only", 0)
        obs_fp_frac = n_fp_obs / n_hot if n_hot else 0
        mean_shuffled_frac = shuffle_control(
            hot_bt, args.genome, mark_bts, mark_of_class, mark_names,
            n_iter=args.n_shuffle, seed=args.seed, outdir=args.outdir,
            blacklist_bed=args.blacklist,
        )
        fold = obs_fp_frac / mean_shuffled_frac if mean_shuffled_frac > 0 else float("inf")
        print(f"\n[shuffle control result]", file=sys.stderr)
        print(f"  observed fp-candidate fraction   : {obs_fp_frac:.4f}", file=sys.stderr)
        print(f"  mean random-shuffle fraction     : {mean_shuffled_frac:.4f}", file=sys.stderr)
        print(f"  fold enrichment                  : {fold:.2f}x", file=sys.stderr)

    pybedtools.cleanup(remove_all=True)
    print(f"\nDone. Outputs in {args.outdir}/", file=sys.stderr)


if __name__ == "__main__":
    main()
