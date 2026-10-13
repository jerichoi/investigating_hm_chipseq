library(rtracklayer)
library(ChIPRegions)

script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
                       error = function(e) tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
                                                    error = function(e) getwd()))
root <- normalizePath(file.path(script_dir, ".."), winslash = "/")   # parent folder, which holds data/

mark <- "H3K27me3"
meta <- read.csv(file.path(root, "data/histone_manifests/histone_experiments_A549.csv"))
meta <- meta[meta$mark == mark, ][1:2, ]
Bedlist <- lapply(file.path(root, meta$fname), import, format = "narrowPeak")

UnifiedSequence <- Unify(Bedlist)
Chrom <- ChromGet(UnifiedSequence)

BICreal <- sapply(1:6, function(k) hmm.discnp::hmm(Chrom[["ChromSplit"]], K = k, itmax = 1000)$BIC)
BestK <- which.min(BICreal)
if (BestK == 1) stop("K = 1 is best: no regional bias between these experiments")

BiasHMM <- hmm.discnp::hmm(Chrom[["ChromSplit"]], K = BestK, itmax = 1000)
BiasViterbi <- hmm.discnp::viterbi(Chrom[["ChromSplit"]], BiasHMM)
ContiguousRegions <- GetRegions(Chrom[["ChromSplit"]], BiasViterbi, Chrom[["Unified Sequence"]])

Rho <- BiasHMM$Rho.matrix
SperUn <- apply(Rho[rownames(Rho) %in% c("TRUE FALSE", "FALSE TRUE"), , drop = FALSE], 2, max)
WhichState <- which.max(SperUn)
BiasedRegions <- ContiguousRegions[ContiguousRegions$State == WhichState]

export(BiasedRegions, file.path(script_dir, paste0(mark, "_biased_regions.bed")))
export(ContiguousRegions, file.path(script_dir, paste0(mark, "_all_regions.bed")))

cat("Best K:", BestK, "\n")
cat("Share of each state's peaks found by only one experiment:\n"); print(round(SperUn, 3))
cat("Biased state:", WhichState, "\n")
cat(length(BiasedRegions), "biased of", length(ContiguousRegions), "regions;",
    round(sum(width(BiasedRegions)) / 1e6), "Mb biased\n")