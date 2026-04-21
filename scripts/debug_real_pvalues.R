#!/usr/bin/env Rscript
# Re-run the F-test but compute ACTUAL p-values from the F statistic,
# not the (mislabeled) spectral power.
#
# Then check: are null p-values uniform on [0,1] when we use
# phase-independent psites? If yes, we've identified 3 distinct bugs:
#   1. shift(-phase) flattens real signal (tautology)
#   2. merge() reorders phaseshift, scrambling row alignment
#   3. "p.value" column is actually spectral power, not a p-value

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(dplyr)
  library(multitaper)
})

data(chr22_anno)
data(rpfs)
data(offsets_df)

set.seed(1)

# Correct F-test per vector
ftest_with_real_pval <- function(psit, k = 24, bw = 12) {
  psit <- as.numeric(psit)
  if (length(psit) < 25) {
    remain  <- 50 - length(psit)
    halfrmn <- as.integer(remain / 2)
    psit <- c(rep(0, halfrmn), psit, rep(0, remain %% 2 + halfrmn))
  }
  padding <- if (length(psit) < 512) 1024 else "default"
  sv <- multitaper::dpss(n = length(psit), k = k, nw = bw)
  rs <- multitaper::spec.mtm(as.ts(psit), k = k, nw = bw, nFFT = padding,
         centreWithSlepians = TRUE, Ftest = TRUE,
         maxAdaptiveIterations = 100, returnZeroFreq = FALSE,
         plot = FALSE, dpssIN = sv)
  ci <- which.min(abs(rs$freq - 1/3))
  Fstat <- rs$mtm$Ftest[ci]
  # spec.mtm's Ftest is a standard F with (2, 2k-2) df
  pval <- pf(Fstat, df1 = 2, df2 = 2*k - 2, lower.tail = FALSE)
  c(Fstat = Fstat, p = pval)
}

# Build psites two ways -------------------------------------------------
# (A) current pipeline
ps_pipeline_real <- suppressWarnings(get_psite_gr(rpfs, offsets_df, chr22_anno))

# randomize positions
tx_lens  <- seqlengths(rpfs)
rpf_lens <- tx_lens[as.character(seqnames(rpfs))]
max_start <- pmax(rpf_lens - rpfs$readlen + 1L, 1L)
rs_start <- pmin(pmax(1L, as.integer(round(runif(length(rpfs), 1, max_start)))), max_start)
rpfs_rand <- GRanges(
  seqnames(rpfs),
  IRanges(start = rs_start, width = rpfs$readlen),
  strand(rpfs), seqinfo = seqinfo(rpfs)
)
mcols(rpfs_rand) <- mcols(rpfs)

ps_pipeline_rand <- suppressWarnings(get_psite_gr(rpfs_rand, offsets_df, chr22_anno))

# (B) phase-independent
rl_off <- offsets_df %>% filter(phase == 0) %>% distinct(readlen, p_offset)
build_phase_indep <- function(rr, anno, rl) {
  orfs <- anno$trspacecds
  ov <- findOverlaps(rr, orfs, select="first", ignore.strand=TRUE)
  rr$orf <- S4Vectors::Rle(names(orfs)[ov])
  rr <- rr[!is.na(rr$orf)]
  off <- as.data.frame(mcols(rr)[,"readlen",drop=FALSE]) %>%
    left_join(rl, by="readlen")
  rr <- rr[!is.na(off$p_offset)]
  off <- off[!is.na(off$p_offset),]
  ps <- rr %>% resize(1,"start") %>% shift(off$p_offset)
  sl <- seqlengths(ps); lens <- sl[as.character(seqnames(ps))]
  ps[(start(ps) >= 1) & (is.na(lens) | start(ps) <= lens)]
}
ps_indep_real <- build_phase_indep(rpfs, chr22_anno, rl_off)
ps_indep_rand <- build_phase_indep(rpfs_rand, chr22_anno, rl_off)

# Shared tester: build per-orf psite coverage then run F-test with real p
run_real_ftest <- function(ps, anno) {
  orfs <- anno$trspacecds
  orfs_with_ps <- intersect(unique(as.character(mcols(ps)$orf)), names(orfs))
  orfs_sub <- orfs[orfs_with_ps]
  mapped <- GenomicFeatures::mapToTranscripts(ps, orfs_sub, ignore.strand = TRUE)
  mapped <- mapped[as.character(mcols(ps)$orf)[mapped$xHits] ==
                   names(orfs_sub)[mapped$transcriptsHits]]
  cov <- coverage(mapped)
  cov <- cov[sum(cov > 0) > 1]
  stats <- lapply(cov, function(v) {
    tryCatch(ftest_with_real_pval(v),
             error = function(e) c(Fstat = NA_real_, p = NA_real_))
  })
  res <- do.call(rbind, stats) %>% as.data.frame()
  res$orf_id <- names(cov)
  res
}

message("Running real-pvalue F-test: pipeline real ...")
rp_pr <- run_real_ftest(ps_pipeline_real, chr22_anno)
message("Running real-pvalue F-test: pipeline randomized ...")
rp_prd <- run_real_ftest(ps_pipeline_rand, chr22_anno)
message("Running real-pvalue F-test: phase-indep real ...")
rp_ir <- run_real_ftest(ps_indep_real, chr22_anno)
message("Running real-pvalue F-test: phase-indep randomized ...")
rp_ird <- run_real_ftest(ps_indep_rand, chr22_anno)

summarize_p <- function(label, df) {
  p <- df$p[!is.na(df$p)]
  cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  tested        : %d\n", length(p)))
  cat(sprintf("  median p      : %.4g\n", median(p)))
  cat(sprintf("  frac p < 0.05 : %.1f%%\n", 100*mean(p < 0.05)))
  cat(sprintf("  frac p < 0.01 : %.1f%%\n", 100*mean(p < 0.01)))
  cat("\n")
}

cat("\n==== Real F-test p-values ====\n\n")
summarize_p("Pipeline psites, REAL reads",              rp_pr)
summarize_p("Pipeline psites, RANDOMIZED reads (null)", rp_prd)
summarize_p("Phase-indep psites, REAL reads",              rp_ir)
summarize_p("Phase-indep psites, RANDOMIZED reads (null)", rp_ird)

# Save histograms for visual inspection
library(ggplot2)
plot_df <- bind_rows(
  rp_pr  %>% mutate(cond = "Pipeline",    data = "Real"),
  rp_prd %>% mutate(cond = "Pipeline",    data = "Randomized"),
  rp_ir  %>% mutate(cond = "Phase-indep", data = "Real"),
  rp_ird %>% mutate(cond = "Phase-indep", data = "Randomized")
) %>% filter(!is.na(p))

p_plot <- ggplot(plot_df, aes(p)) +
  geom_histogram(boundary = 0, binwidth = 0.02,
                 fill = "steelblue", colour = "white") +
  facet_grid(cond ~ data) +
  labs(title = "ACTUAL F-test p-values (not the mislabeled spectral power)",
       subtitle = "Calibrated null should be UNIFORM on [0,1]",
       x = "p", y = "count") +
  theme_bw(base_size = 12)
ggsave("scripts/real_pvalue_histograms.png", p_plot, width = 10, height = 6, dpi = 150)
message("Plot saved: scripts/real_pvalue_histograms.png")
