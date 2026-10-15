library(rtracklayer)
library(ChIPRegions)

targets <- list(K562 = c("H3K27me3", "H3K9me3"))   # cell line = its marks; add more cell lines as extra entries

script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) getwd()))
root <- normalizePath(file.path(script_dir, ".."), winslash = "/")

for (cell in names(targets)) {
  manifest <- read.csv(file.path(root, "data/histone_manifests", paste0("histone_experiments_", cell, ".csv")))

  for (mark in targets[[cell]]) {
    cat("\n=====", cell, mark, "=====\n")
    meta <- manifest[manifest$mark == mark, ]
    if (nrow(meta) < 2) { cat("fewer than 2 experiments, skipped\n"); next }
    Bedlist <- lapply(file.path(root, meta$fname), import, format = "narrowPeak")
    names(Bedlist) <- meta$experiment_id
    Bedlist <- lapply(Bedlist, function(x) { seqlevels(x, pruning.mode = "coarse") <- paste0("chr", c(1:22, "X", "Y")); x })   # same chromosomes in every file

    UnifiedSequence <- Unify(Bedlist)
    Chrom <- ChromGet(UnifiedSequence)

    BICreal <- sapply(1:6, function(k) hmm.discnp::hmm(Chrom[["ChromSplit"]], K = k, itmax = 1000)$BIC)
    BestK <- which.min(BICreal)

    cat("Experiments:", length(Bedlist), " Peak unions:", length(UnifiedSequence), "\n")
    cat("BIC, K = 1..6:", round(BICreal), "\n")
    cat("Best K:", BestK, "\n")
    if (BestK == 1) { cat("No regional bias\n"); next }

    BiasHMM <- hmm.discnp::hmm(Chrom[["ChromSplit"]], K = BestK, itmax = 1000)
    BiasViterbi <- hmm.discnp::viterbi(Chrom[["ChromSplit"]], BiasHMM)
    ContiguousRegions <- GetRegions(Chrom[["ChromSplit"]], BiasViterbi, Chrom[["Unified Sequence"]])

    Rho <- BiasHMM$Rho.matrix
    Calls <- do.call(rbind, strsplit(rownames(Rho), " ")) == "TRUE"
    Only <- matrix(sapply(seq_along(Bedlist), function(k) colSums(Rho[Calls[, k] & rowSums(Calls) == 1, , drop = FALSE])), nrow = BestK)
    SperUn <- apply(Only, 1, max)
    WhichState <- which.max(SperUn)
    BiasedRegions <- ContiguousRegions[ContiguousRegions$State == WhichState]

    export(BiasedRegions, file.path(script_dir, paste0(cell, "_", mark, "_biased.bed")))
    export(ContiguousRegions, file.path(script_dir, paste0(cell, "_", mark, "_all.bed")))

    cat("Share of each state's peaks found by only one experiment:\n"); print(round(SperUn, 3))
    cat("Biased state:", WhichState, "\n")
    cat(length(BiasedRegions), "biased of", length(ContiguousRegions), "regions;", round(sum(width(BiasedRegions)) / 1e6), "Mb biased\n")

    cat("Emission matrix, Rho (rows = which experiments called a region, columns = states):\n"); print(round(Rho, 3))
    Called <- t(Calls) %*% Rho; rownames(Called) <- names(Bedlist)
    cat("Share of each state's regions called by each experiment:\n"); print(round(Called, 3))
    cat("Share of all regions in each state:\n"); print(round(BiasHMM$ispd, 3))
  }
}
