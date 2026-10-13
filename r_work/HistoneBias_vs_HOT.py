import numpy as np
import pandas as pd
from pathlib import Path
from pybedtools import BedTool

# settings
out_dir      = Path("r_work/output")                   # written by the R script
tf_loci_fn   = "data/merged_tf_counts.bed"
blacklist_fn = "data/hg38-blacklist.v2.bed"
hot_fns      = {"K562": "data/hot_regions.clean.bed"}  # add HOT loci for other cell lines here
default_hot  = "K562"                                  # used for cell lines without their own HOT loci
n_perm       = 200
require_clean_null = True                              # skip pairs whose shuffled control picked K > 1
rng = np.random.default_rng(1)

# TF-bound loci (blacklist removed once), and hot / not-hot sets per HOT cell line
tf_loci = BedTool(tf_loci_fn).intersect(BedTool(blacklist_fn), v=True).saveas()
universe = {}
def get_universe(cell):
    if cell not in universe:
        hot = BedTool(hot_fns[cell])
        universe[cell] = (hot, tf_loci.intersect(hot, v=True).saveas())
    return universe[cell]

def odds_ratio(regions, hot, not_hot):
    a = hot.intersect(regions, u=True).count()
    c = not_hot.intersect(regions, u=True).count()
    b, d = hot.count() - a, not_hot.count() - c
    if 0 in (a, b, c, d):
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    return a * d / (b * c)

def block_null(blocks_by_chrom, target_bp):
    # relabel random blocks as "biased": same chromosomes, same total bp per chromosome
    picked = []
    for chrom, target in target_bp.items():
        blocks, bp, chosen = blocks_by_chrom[chrom], 0, []
        for i in rng.permutation(len(blocks)):
            chosen.append(blocks[i]); bp += blocks[i][2] - blocks[i][1]
            if bp >= target:
                if bp - target > (blocks[i][2] - blocks[i][1]) / 2 and len(chosen) > 1:
                    chosen.pop()
                break
        picked += chosen
    return BedTool("\n".join(f"{c}\t{s}\t{e}" for c, s, e in picked), from_string=True)

def p_values(obs, null):
    pe = (1 + sum(n >= obs for n in null)) / (1 + len(null))
    pd_ = (1 + sum(n <= obs for n in null)) / (1 + len(null))
    return pe, pd_, min(1, 2 * min(pe, pd_))

def bh(p):                                              # Benjamini-Hochberg FDR
    p = np.asarray(p); order = np.argsort(p); n = len(p)
    q = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(q, 1)
    return out

# one test per biased pair
summary = pd.read_csv(out_dir / "bias_summary.csv")
todo = summary.dropna(subset=["biased_bed"])
if require_clean_null:
    todo = todo[todo["K_shuffled"] == 1]

rows = []
for r in todo.itertuples():
    hot_cell = r.cell_line if r.cell_line in hot_fns else default_hot
    hot, not_hot = get_universe(hot_cell)
    biased = BedTool(str(out_dir / r.biased_bed))
    blocks = [(iv.chrom, iv.start, iv.end) for iv in BedTool(str(out_dir / r.all_bed))]
    blocks_by_chrom = {}
    for b in blocks: blocks_by_chrom.setdefault(b[0], []).append(b)
    target_bp = {}
    for iv in biased: target_bp[iv.chrom] = target_bp.get(iv.chrom, 0) + iv.length

    obs = odds_ratio(biased, hot, not_hot)
    null_b = [odds_ratio(block_null(blocks_by_chrom, target_bp), hot, not_hot) for _ in range(n_perm)]
    null_g = [odds_ratio(biased.shuffle(genome="hg38", chrom=True, excl=blacklist_fn, seed=i).saveas(), hot, not_hot)
              for i in range(n_perm)]
    pe_b, pd_b, p2_b = p_values(obs, null_b)
    pe_g, pd_g, p2_g = p_values(obs, null_g)
    rows.append(dict(cell_line=r.cell_line, mark=r.mark, exp_A=r.exp_A, exp_B=r.exp_B,
                     hot_cell=hot_cell, matched=(hot_cell == r.cell_line),
                     share_one_sided=r.share_one_sided, biased_Mb=r.biased_Mb, OR=obs,
                     null_median_blocks=np.median(null_b), p_enrich_blocks=pe_b, p_deplete_blocks=pd_b, p_blocks=p2_b,
                     null_median_genome=np.median(null_g), p_enrich_genome=pe_g, p_deplete_genome=pd_g, p_genome=p2_g))

res = pd.DataFrame(rows)
res["q_blocks"] = bh(res["p_blocks"])
res["q_genome"] = bh(res["p_genome"])
res.to_csv(out_dir / "hot_overlap_results.csv", index=False)

print(res.groupby(["cell_line", "mark", "matched"]).agg(
    pairs_tested=("OR", "size"), median_OR=("OR", "median"),
    sig_blocks=("q_blocks", lambda q: int((q < 0.05).sum())),
    sig_genome=("q_genome", lambda q: int((q < 0.05).sum()))))
