#!/usr/bin/env Rscript
# Re-run the null test using SHIFTS COMPUTED ON REAL DATA,
# applied to both real and randomized psites.
#
# This is the realistic workflow: calibrate shift parameters on
# data with strong periodicity, then apply those fixed parameters
# to any test set (including novel ORFs or null data).

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
# Step 1: Compute shiftdf ONCE on real data
# -----------------------------------------------------------------
compute_shiftdf <- function(rpfs, offsets_df, anno) {
  orfs <- anno$trspacecds
  ov <- findOverlaps(rpfs, orfs, select = "first", ignore.strand = TRUE)
  rpfs$orf   <- S4Vectors::Rle(names(orfs)[ov])
  rpfs$phase <- (start(rpfs) - start(orfs)[ov]) %% 3
  rpfs$p_offset <- as.data.frame(mcols(rpfs)[, c("readlen", "phase")]) %>%
    left_join(offsets_df %>% select(readlen, phase, p_offset),
              by = c("readlen", "phase")) %>% .$p_offset
  rpfs <- subset(rpfs, !is.na(p_offset))
  psites <- rpfs %>% resize(1, "start") %>% shift(.$p_offset)
  sl <- seqlengths(psites); lens <- sl[as.character(seqnames(psites))]
  psites <- psites[(start(psites) >= 1) &
                   (is.na(lens) | start(psites) <= lens)]
  as.data.frame(mcols(psites)[, c("phase", "readlen")]) %>%
    count(phase, readlen) %>%
    group_by(readlen) %>% mutate(shft = rank(-n) - 1) %>% ungroup() %>%
    select(phase, readlen, shft)
}

message("Computing shiftdf from real data ...")
shiftdf_fixed <- compute_shiftdf(rpfs, offsets_df, chr22_anno)
cat("Fixed shiftdf (from real data):\n")
print(as.data.frame(shiftdf_fixed %>% arrange(readlen, phase)))

# -----------------------------------------------------------------
# Step 2: psite builder that applies FIXED shift parameters
# -----------------------------------------------------------------
get_psites_apply_fixed <- function(rpfs, offsets_df, anno, shiftdf) {
  orfs <- anno$trspacecds
  ov <- findOverlaps(rpfs, orfs, select = "first", ignore.strand = TRUE)
  rpfs$orf   <- S4Vectors::Rle(names(orfs)[ov])
  rpfs$phase <- (start(rpfs) - start(orfs)[ov]) %% 3
  rpfs$p_offset <- as.data.frame(mcols(rpfs)[, c("readlen", "phase")]) %>%
    left_join(offsets_df %>% select(readlen, phase, p_offset),
              by = c("readlen", "phase")) %>% .$p_offset
  rpfs <- subset(rpfs, !is.na(p_offset))
  psites <- rpfs %>% resize(1, "start") %>% shift(.$p_offset)
  sl <- seqlengths(psites); lens <- sl[as.character(seqnames(psites))]
  psites <- psites[(start(psites) >= 1) &
                   (is.na(lens) | start(psites) <= lens)]
  # apply fixed shiftdf, order-preserving
  shft_per_row <- as.data.frame(mcols(psites)[, c("phase", "readlen")]) %>%
    left_join(shiftdf, by = c("phase", "readlen")) %>% .$shft
  psites %>% shift(-.$phase) %>% shift(shft_per_row)
}

# -----------------------------------------------------------------
# Step 3: Randomize and apply
# -----------------------------------------------------------------
tx_lens <- seqlengths(rpfs)
rpf_lens <- tx_lens[as.character(seqnames(rpfs))]
max_start <- pmax(rpf_lens - rpfs$readlen + 1L, 1L)
rs_start <- pmin(pmax(1L, as.integer(round(runif(length(rpfs), 1, max_start)))), max_start)
rpfs_rand <- GRanges(
  seqnames(rpfs),
  IRanges(start = rs_start, width = rpfs$readlen),
  strand(rpfs), seqinfo = seqinfo(rpfs)
)
mcols(rpfs_rand) <- mcols(rpfs)

message("Building psites with FIXED shifts (real) ...")
ps_real <- suppressWarnings(
  get_psites_apply_fixed(rpfs, offsets_df, chr22_anno, shiftdf_fixed)
)
message("Building psites with FIXED shifts (randomized) ...")
ps_rand <- suppressWarnings(
  get_psites_apply_fixed(rpfs_rand, offsets_df, chr22_anno, shiftdf_fixed)
)

# -----------------------------------------------------------------
# Step 4: F-test with real p-values
# -----------------------------------------------------------------
ftest_real_p <- function(psit, k = 24, bw = 12) {
  psit <- as.numeric(psit)
  if (length(psit) < 25) {
    remain <- 50 - length(psit); halfrmn <- as.integer(remain / 2)
    psit <- c(rep(0, halfrmn), psit, rep(0, remain %% 2 + halfrmn))
  }
  padding <- if (length(psit) < 512) 1024 else "default"
  sv <- multitaper::dpss(n = length(psit), k = k, nw = bw)
  rs <- multitaper::spec.mtm(as.ts(psit), k = k, nw = bw, nFFT = padding,
         centreWithSlepians = TRUE, Ftest = TRUE,
         maxAdaptiveIterations = 100, returnZeroFreq = FALSE,
         plot = FALSE, dpssIN = sv)
  ci <- which.min(abs(rs$freq - 1/3))
  pf(rs$mtm$Ftest[ci], df1 = 2, df2 = 2 * k - 2, lower.tail = FALSE)
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

message("F-test with real p-values: FIXED-SHIFT pipeline, real ...")
pr <- run_ftest(ps_real, chr22_anno)
message("F-test with real p-values: FIXED-SHIFT pipeline, randomized ...")
pn <- run_ftest(ps_rand, chr22_anno)

summ <- function(label, df) {
  p <- df$p[!is.na(df$p)]
  cat(sprintf("%-55s  n=%4d  median=%.4g  p<0.05=%.1f%%  p<0.01=%.1f%%\n",
              label, length(p), median(p),
              100 * mean(p < 0.05), 100 * mean(p < 0.01)))
}

cat("\n==== Pipeline with SHIFTS FROM REAL DATA, properly aligned ====\n")
summ("FIXED-SHIFT pipeline, REAL reads",       pr)
summ("FIXED-SHIFT pipeline, RANDOMIZED reads", pn)

# Frame distributions
orf_st <- function(ids) start(chr22_anno$trspacecds[ids])
fr_r <- (start(ps_real) - orf_st(as.character(mcols(ps_real)$orf))) %% 3
fr_n <- (start(ps_rand) - orf_st(as.character(mcols(ps_rand)$orf))) %% 3
cat("\nFrame distribution (FIXED-SHIFT pipeline):\n")
cat(sprintf("REAL       : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100 * mean(fr_r == 0), 100 * mean(fr_r == 1), 100 * mean(fr_r == 2)))
cat(sprintf("RANDOMIZED : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100 * mean(fr_n == 0), 100 * mean(fr_n == 1), 100 * mean(fr_n == 2)))

cat("\n\n==== Comparison summary across all variants ====\n")
cat(sprintf("%-55s  %-10s  %-10s\n", "method", "real p<.05", "null p<.05"))
cat(sprintf("%-55s  %-10s  %-10s\n", "-----", "----------", "----------"))
cat(sprintf("%-55s  %-10s  %-10s\n",
            "Original pipeline (merge bug, spec_coef as pval)",
            "n/a (label wrong)", "n/a"))
cat(sprintf("%-55s  %-10s  %-10s\n",
            "Original pipeline + real p-values",
            "72%", "77%"))
cat(sprintf("%-55s  %-10s  %-10s\n",
            "Merge bug fixed, shifts recomputed on null",
            "21%", "14%"))
cat(sprintf("%-55s  %-10s  %-10s\n",
            "Merge bug fixed, FIXED SHIFTS from real data (this run)",
            sprintf("%.1f%%", 100 * mean(pr$p < 0.05, na.rm = TRUE)),
            sprintf("%.1f%%", 100 * mean(pn$p < 0.05, na.rm = TRUE))))
cat(sprintf("%-55s  %-10s  %-10s\n",
            "Phase-independent (no shift(-phase))",
            "21%", "4.1%"))
