library(rtracklayer)
library(ChIPRegions)

cell_lines <- c("A549")   # needs data/histone_manifests/histone_experiments_<cell_line>.csv
max_pairs  <- 20          # random pairs per mark; NULL = all pairs
Kmax       <- 6
min_peaks  <- 1000
set.seed(1)

script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) getwd()))
root <- normalizePath(file.path(script_dir, ".."), winslash = "/")

Quiet <- function(x) { invisible(capture.output(r <- suppressWarnings(suppressMessages(x)))); r }

BICs <- function(seqs) sapply(1:Kmax, function(k) {
  m <- Quiet(tryCatch(hmm.discnp::hmm(seqs, K = k, itmax = 1000), error = function(e) NULL))
  if (is.null(m)) Inf else m$BIC
})

OptimalK <- function(Bedlist) {
  U <- Unify(Bedlist)
  real <- BICs(Quiet(ChromGet(U))[["ChromSplit"]])
  U$OverlapPattern <- sample(U$OverlapPattern)                 # shuffled control
  shuf <- BICs(Quiet(ChromGet(U))[["ChromSplit"]])
  data.frame(K_real = which.min(real), K_shuffled = which.min(shuf), bic_gain = real[1] - min(real))
}

results <- list()
for (cell in cell_lines) {
  meta <- read.csv(file.path(root, "data/histone_manifests", paste0("histone_experiments_", cell, ".csv")))
  meta <- meta[meta$status == "downloaded", ]
  for (mk in unique(meta$mark)) {
    sub <- meta[meta$mark == mk, ]
    peaks <- setNames(lapply(file.path(root, sub$fname), import, format = "narrowPeak"), sub$experiment_id)
    peaks <- peaks[lengths(peaks) >= min_peaks]
    if (length(peaks) < 2) next
    pairs <- combn(names(peaks), 2, simplify = FALSE)
    if (!is.null(max_pairs) && length(pairs) > max_pairs) pairs <- sample(pairs, max_pairs)
    for (p in pairs) {
      message(cell, " ", mk, ": ", p[1], " vs ", p[2])
      r <- tryCatch(OptimalK(peaks[p]), error = function(e) NULL)
      if (!is.null(r)) results[[length(results) + 1]] <- cbind(cell_line = cell, mark = mk, exp_A = p[1], exp_B = p[2], r)
    }
  }
}

res <- do.call(rbind, results)
res$biased <- res$K_real > 1 & res$K_shuffled == 1
write.csv(res, file.path(script_dir, "optimal_K_pairs.csv"), row.names = FALSE)

summary_tab <- do.call(rbind, lapply(split(res, list(res$cell_line, res$mark), drop = TRUE), function(d)
  data.frame(cell_line = d$cell_line[1], mark = d$mark[1], pairs = nrow(d), biased = sum(d$biased),
             pct_biased = round(100 * mean(d$biased)), median_K = median(d$K_real), max_K = max(d$K_real),
             unreliable = sum(d$K_shuffled > 1))))
rownames(summary_tab) <- NULL
write.csv(summary_tab, file.path(script_dir, "optimal_K_summary.csv"), row.names = FALSE)
print(summary_tab)
