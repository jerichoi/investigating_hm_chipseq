#!/usr/bin/env python3
"""
Find HOT loci that overlap both active and repressive histone marks,
classify them, and print/save the results.

Edit the file paths below to point at your ENCODE downloads, then run:
    python core_classify.py
"""

import pybedtools
from pybedtools import BedTool

# ---------------------------------------------------------------------
# 1. INPUT FILES -- edit these paths
# ---------------------------------------------------------------------
HOT_LOCI = "hot_regions.clean.bed"

ACTIVE_MARKS = {
    "H3K4me3": "H3K4me3.bed.gz",
    "H3K27ac": "H3K27ac.bed.gz",
}

FACULTATIVE_HET_MARKS = {
    "H3K27me3": "H3K27me3.bed.gz",
}

CONSTITUTIVE_HET_MARKS = {
    "H3K9me3": "H3K9me3.bed.gz",
}

OUTDIR = "results"

# ---------------------------------------------------------------------
# 2. LOAD HOT LOCI
# ---------------------------------------------------------------------
hot = BedTool(HOT_LOCI).sort()
n_hot = hot.count()
print(f"Loaded {n_hot} HOT loci from {HOT_LOCI}")

# ---------------------------------------------------------------------
# 3. INTERSECT EACH MARK AGAINST THE HOT LOCI
#    -> for each HOT locus, record whether it overlaps at least one
#       peak from each mark (boolean presence)
# ---------------------------------------------------------------------
def loci_with_overlap(hot_bt, mark_path):
    """Return the set of HOT-locus keys (chrom:start-end) that overlap
    at least one peak in the given mark's BED file."""
    mark_bt = BedTool(mark_path).sort()
    hits = set()
    for iv in hot_bt.intersect(mark_bt, u=True):   # -u: report each HOT locus once if it has >=1 overlap
        hits.add(f"{iv.chrom}:{iv.start}-{iv.end}")
    return hits

print("\nIntersecting marks...")
active_hits = set()
for name, path in ACTIVE_MARKS.items():
    hits = loci_with_overlap(hot, path)
    print(f"  {name:10s} (active)              : {len(hits)} / {n_hot} loci overlap")
    active_hits |= hits

facultative_hits = set()
for name, path in FACULTATIVE_HET_MARKS.items():
    hits = loci_with_overlap(hot, path)
    print(f"  {name:10s} (facultative het)      : {len(hits)} / {n_hot} loci overlap")
    facultative_hits |= hits

constitutive_hits = set()
for name, path in CONSTITUTIVE_HET_MARKS.items():
    hits = loci_with_overlap(hot, path)
    print(f"  {name:10s} (constitutive het)     : {len(hits)} / {n_hot} loci overlap")
    constitutive_hits |= hits

# ---------------------------------------------------------------------
# 4. CLASSIFY EACH LOCUS
# ---------------------------------------------------------------------
classification = {}   # locus_key -> category
interval_lookup = {}  # locus_key -> pybedtools Interval

for iv in hot:
    key = f"{iv.chrom}:{iv.start}-{iv.end}"
    interval_lookup[key] = iv

    is_active = key in active_hits
    is_facultative = key in facultative_hits
    is_constitutive = key in constitutive_hits

    if is_active and is_constitutive:
        category = "dual_constitutive"        # strongest false-positive evidence
    elif is_active and is_facultative:
        category = "dual_facultative_only"     # ambiguous -- possible real bivalency
    elif is_active:
        category = "eu_only"
    elif is_facultative or is_constitutive:
        category = "het_only"
    else:
        category = "neither"

    classification[key] = category

# ---------------------------------------------------------------------
# 5. PRESENT RESULTS
# ---------------------------------------------------------------------
import os
from collections import Counter

os.makedirs(OUTDIR, exist_ok=True)

counts = Counter(classification.values())
print("\n--- Summary ---")
print(f"Total HOT loci: {n_hot}")
for category in ["dual_constitutive", "dual_facultative_only", "eu_only", "het_only", "neither"]:
    n = counts.get(category, 0)
    pct = 100 * n / n_hot if n_hot else 0
    print(f"  {category:22s}: {n:5d} ({pct:5.1f}%)")

# write full annotated table
table_path = os.path.join(OUTDIR, "classified_loci.tsv")
with open(table_path, "w") as out:
    out.write("chrom\tstart\tend\tcategory\n")
    for key, iv in interval_lookup.items():
        out.write(f"{iv.chrom}\t{iv.start}\t{iv.end}\t{classification[key]}\n")
print(f"\nFull table written to {table_path}")

# write one BED file per category
for category in ["dual_constitutive", "dual_facultative_only", "eu_only", "het_only", "neither"]:
    path = os.path.join(OUTDIR, f"{category}.bed")
    with open(path, "w") as out:
        for key, cat in classification.items():
            if cat == category:
                iv = interval_lookup[key]
                out.write(f"{iv.chrom}\t{iv.start}\t{iv.end}\n")
    print(f"  {category}.bed written ({counts.get(category, 0)} loci)")

pybedtools.cleanup(remove_all=True)
