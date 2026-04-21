#!/usr/bin/env Rscript
# Test ftestvect calibration on synthetic short non-periodic vectors.
#
# Concern: with k = 24 tapers and n < 50 data points (padded to 50 with
# zeros), the F-test's nominal (2, 2k - 2) = (2, 46) null distribution
# may not hold — the effective degrees of freedom are bounded by the
# real data length, and short + padded data creates spurious harmonic
# structure at low frequencies that can leak into the 1/3 Hz band.
#
# Simulate vectors of various lengths from a Poisson distribution
# (no periodicity), run ftestvect, and check if p-values are
# approximately uniform on [0, 1].

suppressPackageStartupMessages({
  library(Ribostan)
  library(dplyr)
  library(ggplot2)
})

set.seed(1)

n_sim <- 500
rates <- c(0.5, 1, 2, 5)   # Poisson rate parameters
lens  <- c(15, 30, 60, 120, 300, 600)  # ORF lengths in nt

run_one <- function(n, rate) {
  v <- rpois(n, rate)
  tryCatch(Ribostan:::ftestvect(v),
           error   = function(e) c(Fstat = NA, spec_coef = NA, p.value = NA),
           warning = function(w) suppressWarnings(Ribostan:::ftestvect(v)))
}

grid <- expand.grid(len = lens, rate = rates, rep = seq_len(n_sim))
results <- vapply(seq_len(nrow(grid)),
                  function(i) run_one(grid$len[i], grid$rate[i]),
                  numeric(3))
results <- as.data.frame(t(results))
grid$p <- results$p.value

cat("==== F-test p-value calibration on Poisson (non-periodic) vectors ====\n")
cat("Each cell: median p  /  frac p<0.05  (expected 0.5 / 5% under null)\n\n")
tbl <- grid %>%
  filter(!is.na(p)) %>%
  group_by(len, rate) %>%
  summarise(
    n = n(),
    med = median(p),
    frac05 = mean(p < 0.05),
    frac01 = mean(p < 0.01),
    .groups = "drop"
  )
print(as.data.frame(tbl))

cat("\n==== Aggregated by length ====\n")
agg_len <- grid %>%
  filter(!is.na(p)) %>%
  group_by(len) %>%
  summarise(
    n = n(),
    med = median(p),
    frac05 = mean(p < 0.05),
    frac01 = mean(p < 0.01)
  )
print(as.data.frame(agg_len))

# Histograms by length
p_plot <- ggplot(grid %>% filter(!is.na(p)),
                 aes(x = p)) +
  geom_histogram(boundary = 0, binwidth = 0.05,
                 fill = "steelblue", colour = "white") +
  facet_wrap(~ paste0("len = ", len, " nt"), scales = "free_y") +
  geom_hline(yintercept = n_sim * length(rates) / 20,
             linetype = "dashed", colour = "red") +
  labs(title = "ftestvect p-value distribution on synthetic non-periodic Poisson data",
       subtitle = "Under proper null calibration histograms should be flat (red dashed line)",
       x = "p.value", y = "count") +
  theme_bw(base_size = 11)
ggsave("scripts/short_orf_pval_histograms.png", p_plot,
       width = 10, height = 6, dpi = 150)
cat("\nPlot saved: scripts/short_orf_pval_histograms.png\n")

# ALSO: check what happens with periodic data to confirm the test
# detects real signal
cat("\n==== Sanity: periodic data (spike every 3 positions) ====\n")
make_periodic <- function(n, rate = 2) {
  v <- rep(0L, n)
  v[seq(1, n, by = 3)] <- rpois(length(v[seq(1, n, by = 3)]), rate * 3)
  v
}
sim_periodic <- function(len, rate, nsim = 200) {
  out <- vapply(seq_len(nsim), function(i) {
    v <- make_periodic(len, rate)
    tryCatch(Ribostan:::ftestvect(v)["p.value"],
             error = function(e) NA_real_)
  }, numeric(1))
  out[!is.na(out)]
}
for (ln in c(30, 60, 120, 300)) {
  ps <- sim_periodic(ln, rate = 1, nsim = 200)
  cat(sprintf("  len=%3d  median p=%.3g  frac p<0.05=%.1f%%\n",
              ln, median(ps), 100 * mean(ps < 0.05)))
}
