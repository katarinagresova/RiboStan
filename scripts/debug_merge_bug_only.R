#!/usr/bin/env Rscript
# Verify: if we fix ONLY the merge-ordering bug (keep shift(-phase) +
# shift(shft), but apply shft with correct per-row alignment), does
# the null test become calibrated?
#
# If yes: the merge bug is the real culprit, shift(-phase) is fine.
# If no: shift(-phase) is also problematic.

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

# -----------------------------------------------------------------
# Fixed get_psite_gr: identical to current except the phaseshift uses
# a left_join (order-preserving) instead of merge().
# -----------------------------------------------------------------
get_psite_gr_fixed <- function(rpfs, offsets_df, anno) {
  orfs <- c(anno$trspacecds)
  orfov_ind <- findOverlaps(rpfs, orfs, select = "first", ignore.strand = TRUE)
  rpfs$orf   <- S4Vectors::Rle(names(orfs)[orfov_ind])
  rpfs$phase <- (start(rpfs) - start(orfs)[orfov_ind]) %% 3

  rpfs$p_offset <- as.data.frame(mcols(rpfs)[, c("readlen", "phase")]) %>%
    left_join(offsets_df %>% select(readlen, phase, p_offset),
              by = c("readlen", "phase")) %>%
    .$p_offset
  rpfs <- subset(rpfs, !is.na(p_offset))

  psites <- rpfs %>% resize(1, "start") %>% shift(.$p_offset)
  sl <- seqlengths(psites); lens <- sl[as.character(seqnames(psites))]
  psites <- psites[(start(psites) >= 1) &
                   (is.na(lens) | start(psites) <= lens)]

  shiftdf <- as.data.frame(mcols(psites)[, c("phase", "readlen")]) %>%
    count(phase, readlen) %>%
    group_by(readlen) %>% mutate(shft = rank(-n) - 1) %>% ungroup() %>%
    select(phase, readlen, shft)

  # Order-preserving lookup (the bugfix)
  shft_per_row <- as.data.frame(mcols(psites)[, c("phase", "readlen")]) %>%
    left_join(shiftdf, by = c("phase", "readlen")) %>%
    .$shft

  psites %>% shift(-.$phase) %>% shift(shft_per_row)
}

# -----------------------------------------------------------------
# Randomize RPFs to produce null data
# -----------------------------------------------------------------
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

message("Computing psites: fixed (no merge bug), real ...")
ps_fixed_real <- suppressWarnings(get_psite_gr_fixed(rpfs, offsets_df, chr22_anno))
message("Computing psites: fixed (no merge bug), randomized ...")
ps_fixed_rand <- suppressWarnings(get_psite_gr_fixed(rpfs_rand, offsets_df, chr22_anno))

# Real F-test p-value
ftest_real_p <- function(psit, k = 24, bw = 12) {
  psit <- as.numeric(psit)
  if (length(psit) < 25) {
    remain <- 50 - length(psit); halfrmn <- as.integer(remain/2)
    psit <- c(rep(0, halfrmn), psit, rep(0, remain %% 2 + halfrmn))
  }
  padding <- if (length(psit) < 512) 1024 else "default"
  sv <- multitaper::dpss(n = length(psit), k = k, nw = bw)
  rs <- multitaper::spec.mtm(as.ts(psit), k = k, nw = bw, nFFT = padding,
         centreWithSlepians = TRUE, Ftest = TRUE,
         maxAdaptiveIterations = 100, returnZeroFreq = FALSE,
         plot = FALSE, dpssIN = sv)
  ci <- which.min(abs(rs$freq - 1/3))
  pf(rs$mtm$Ftest[ci], df1 = 2, df2 = 2*k - 2, lower.tail = FALSE)
}

run_ftest <- function(ps, anno) {
  orfs <- anno$trspacecds
  sub  <- intersect(unique(as.character(mcols(ps)$orf)), names(orfs))
  orfs_sub <- orfs[sub]
  m <- GenomicFeatures::mapToTranscripts(ps, orfs_sub, ignore.strand = TRUE)
  m <- m[as.character(mcols(ps)$orf)[m$xHits] == names(orfs_sub)[m$transcriptsHits]]
  cov <- coverage(m); cov <- cov[sum(cov > 0) > 1]
  p <- vapply(cov, function(v) tryCatch(ftest_real_p(v),
                                         error = function(e) NA_real_),
              numeric(1))
  data.frame(orf_id = names(cov), p = p, row.names = NULL)
}

message("F-test: fixed-pipeline real ...")
rf_real <- run_ftest(ps_fixed_real, chr22_anno)
message("F-test: fixed-pipeline randomized ...")
rf_rand <- run_ftest(ps_fixed_rand, chr22_anno)

summ <- function(label, df) {
  p <- df$p[!is.na(df$p)]
  cat(sprintf("%-45s  n=%4d  median=%.4g  p<0.05=%.1f%%  p<0.01=%.1f%%\n",
              label, length(p), median(p),
              100*mean(p < 0.05), 100*mean(p < 0.01)))
}

cat("\n==== Fixed pipeline (shift(-phase) + shift(shft) with correct alignment) ====\n")
summ("FIXED pipeline, REAL reads",       rf_real)
summ("FIXED pipeline, RANDOMIZED reads", rf_rand)

cat("\n==== Comparison summary ====\n")
cat("If only the merge bug matters, FIXED randomized should give\n")
cat("uniform p-values (median ~0.5, p<0.05 ~5%).\n")
cat("If shift(-phase) is also problematic, we'll still see over-calling.\n")

# Frame distribution check
orf_st <- function(ids) start(chr22_anno$trspacecds[ids])
fr_r <- (start(ps_fixed_real) -
         orf_st(as.character(mcols(ps_fixed_real)$orf))) %% 3
fr_n <- (start(ps_fixed_rand) -
         orf_st(as.character(mcols(ps_fixed_rand)$orf))) %% 3
cat(sprintf("\nFrame distribution, FIXED pipeline:\n"))
cat(sprintf("REAL       : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100*mean(fr_r==0), 100*mean(fr_r==1), 100*mean(fr_r==2)))
cat(sprintf("RANDOMIZED : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100*mean(fr_n==0), 100*mean(fr_n==1), 100*mean(fr_n==2)))
