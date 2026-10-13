library(rtracklayer)
library(ChIPRegions)

## settings -------------------------------------------------------------------
cell_lines <- c("A549")   # needs data/histone_manifests/histone_experiments_<cell_line>.csv
marks      <- NULL        # NULL = every mark in the manifest, or e.g. c("H3K27me3", "H3K4me3")
max_pairs  <- 20          # random pairs per mark; NULL = all pairs (H3K27ac has 190)
Kmax       <- 6
min_peaks  <- 1000        # drop experiments with fewer peaks than this
match_peak_numbers <- FALSE   # TRUE = keep each pair's top-N peaks by signal, N = the smaller peak count
set.seed(1)

script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) getwd()))
root    <- normalizePath(file.path(script_dir, ".."), winslash = "/")   # parent folder, which holds data/
out_dir <- file.path(script_dir, "output")

## functions ------------------------------------------------------------------
Quiet <- function(x) { invisible(capture.output(r <- suppressWarnings(suppressMessages(x)))); r }

FitBIC <- function(seqs) {
  fits <- lapply(1:Kmax, function(k)
    Quiet(tryCatch(hmm.discnp::hmm(seqs, K = k, itmax = 1000), error = function(e) NULL)))
  list(fits = fits, bic = sapply(fits, function(m) if (is.null(m)) Inf else m$BIC))
}

BiasPair <- function(peaks, ids, reldir) {
  if (match_peak_numbers) {
    n <- min(lengths(peaks))
    peaks <- lapply(peaks, function(p) p[order(-p$signalValue)][seq_len(n)])
  }
  U <- Unify(peaks)
  Chrom <- Quiet(ChromGet(U))
  real <- FitBIC(Chrom[["ChromSplit"]])
  K <- which.min(real$bic)

  U$OverlapPattern <- sample(U$OverlapPattern)                  # shuffled control
  shuf <- FitBIC(Quiet(ChromGet(U))[["ChromSplit"]])

  row <- data.frame(exp_A = ids[1], exp_B = ids[2],
                    peaks_A = length(peaks[[1]]), peaks_B = length(peaks[[2]]), unions = length(U),
                    K_real = K, K_shuffled = which.min(shuf$bic), bic_gain = real$bic[1] - min(real$bic),
                    share_one_sided = NA, n_biased_regions = NA, biased_Mb = NA,
                    biased_bed = NA, all_bed = NA)
  if (K == 1) return(row)                                       # no regional bias in this pair

  m <- real$fits[[K]]
  Viterbi <- hmm.discnp::viterbi(Chrom[["ChromSplit"]], m)
  Regions <- GetRegions(Chrom[["ChromSplit"]], Viterbi, Chrom[["Unified Sequence"]])

  Rho <- m$Rho.matrix
  SperUn <- apply(Rho[rownames(Rho) %in% c("TRUE FALSE", "FALSE TRUE"), , drop = FALSE], 2, max)
  Biased <- Regions[Regions$State == which.max(SperUn)]

  dir.create(file.path(out_dir, reldir), recursive = TRUE, showWarnings = FALSE)
  tag <- paste0(ids[1], "_", ids[2])
  row$biased_bed <- file.path(reldir, paste0(tag, "_biased.bed"))
  row$all_bed    <- file.path(reldir, paste0(tag, "_all.bed"))
  export(Biased,  file.path(out_dir, row$biased_bed))
  export(Regions, file.path(out_dir, row$all_bed))
  row$share_one_sided <- max(SperUn)
  row$n_biased_regions <- length(Biased)
  row$biased_Mb <- sum(width(Biased)) / 1e6
  row
}

## run ------------------------------------------------------------------------
results <- list()
for (cell in cell_lines) {
  meta <- read.csv(file.path(root, "data/histone_manifests", paste0("histone_experiments_", cell, ".csv")))
  meta <- meta[meta$status == "downloaded", ]
  for (mk in if (is.null(marks)) unique(meta$mark) else marks) {
    sub <- meta[meta$mark == mk, ]
    peaks <- setNames(lapply(file.path(root, sub$fname), import, format = "narrowPeak"), sub$experiment_id)
    peaks <- peaks[lengths(peaks) >= min_peaks]
    if (length(peaks) < 2) next
    pairs <- combn(names(peaks), 2, simplify = FALSE)
    if (!is.null(max_pairs) && length(pairs) > max_pairs) pairs <- sample(pairs, max_pairs)

    for (p in pairs) {
      message(cell, " ", mk, ": ", p[1], " vs ", p[2])
      row <- tryCatch(BiasPair(peaks[p], p, file.path(cell, mk)),
                      error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
      if (!is.null(row)) results[[length(results) + 1]] <- cbind(cell_line = cell, mark = mk, row)
    }
    write.csv(do.call(rbind, results), file.path(out_dir, "bias_summary.csv"), row.names = FALSE)
  }
}

res <- do.call(rbind, results)
print(aggregate(cbind(pairs = 1, biased = K_real > 1, unreliable = K_shuffled > 1) ~ cell_line + mark,
                data = res, FUN = sum))
cat("Files written to:", out_dir, "\n")
