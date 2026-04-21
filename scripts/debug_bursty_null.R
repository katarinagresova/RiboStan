#!/usr/bin/env Rscript
# A harder null test: bursty/clustered coverage. Real uORFs often have
# reads concentrated at a few positions (start codon peak, stall sites)
# without true 3-nt periodicity. Would ftestvect false-positive here?

suppressPackageStartupMessages({
  library(Ribostan)
  library(dplyr)
})

set.seed(2)

# Several burstiness models
bursty_uniform_k <- function(n, k) {
  # place k random nonzero positions, each with Poisson(3)
  v <- rep(0L, n)
  pos <- sample.int(n, size = min(k, n))
  v[pos] <- rpois(length(pos), 3)
  v
}
bursty_contiguous <- function(n, k) {
  # k consecutive positions, Poisson(3)
  v <- rep(0L, n)
  if (k >= n) return(rpois(n, 3))
  start <- sample.int(n - k + 1, 1)
  v[start:(start + k - 1)] <- rpois(k, 3)
  v
}
bursty_two_peaks <- function(n, total = 20) {
  # all counts concentrated at two random positions
  v <- rep(0L, n)
  pos <- sample.int(n, 2)
  v[pos] <- rmultinom(1, total, prob = c(0.5, 0.5))
  v
}

do_sim <- function(generator, label, n_sim = 500, lens = c(30, 60, 120)) {
  out <- lapply(lens, function(L) {
    ps <- vapply(seq_len(n_sim), function(i) {
      v <- generator(L)
      tryCatch(Ribostan:::ftestvect(v)["p.value"],
               error = function(e) NA_real_)
    }, numeric(1))
    ps <- ps[!is.na(ps)]
    data.frame(model = label, len = L,
               n = length(ps),
               med = median(ps),
               frac05 = mean(ps < 0.05),
               frac01 = mean(ps < 0.01))
  })
  bind_rows(out)
}

results <- bind_rows(
  do_sim(function(L) bursty_uniform_k(L, k = 5),     "5 random peaks"),
  do_sim(function(L) bursty_uniform_k(L, k = 10),    "10 random peaks"),
  do_sim(function(L) bursty_contiguous(L, k = 10),   "10 contiguous peaks"),
  do_sim(function(L) bursty_two_peaks(L, total = 20),"2 random peaks (concentrated)"),
  do_sim(function(L) rpois(L, 2),                    "Poisson (smooth, not bursty)")
)

cat("==== p-value calibration under BURSTY non-periodic models ====\n")
cat("Expected under correct null: median ~0.5, frac p<0.05 ~5%\n\n")
print(results %>% as.data.frame())

cat("\n")
cat("If the test is robust, all these should show frac p<0.05 ~= 5%.\n")
cat("False positives on these would mean the test mistakes clustering\n")
cat("for periodicity in short/sparse data.\n")
