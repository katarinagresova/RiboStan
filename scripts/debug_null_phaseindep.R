#!/usr/bin/env Rscript
# Extension of the null test: also run it on PHASE-INDEPENDENT psites
# (single offset per readlen, no shift(-phase), no shift(shft)).
#
# If this version gives ~uniform p-values on random reads, we've
# identified exactly which step of get_psite_gr is generating the
# false positives.

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(dplyr)
})

data(chr22_anno)
data(rpfs)
data(offsets_df)

set.seed(1)

# -----------------------------------------------------------------
# Randomize positions
# -----------------------------------------------------------------
tx_lens <- seqlengths(rpfs)
rpf_lens <- tx_lens[as.character(seqnames(rpfs))]
max_start <- pmax(rpf_lens - rpfs$readlen + 1L, 1L)
rand_start <- pmin(
  pmax(1L, as.integer(round(runif(length(rpfs), 1, max_start)))),
  max_start
)
rpfs_rand <- GRanges(
  seqnames = seqnames(rpfs),
  ranges   = IRanges(start = rand_start, width = rpfs$readlen),
  strand   = strand(rpfs),
  seqinfo  = seqinfo(rpfs)
)
mcols(rpfs_rand) <- mcols(rpfs)

# -----------------------------------------------------------------
# Phase-independent psite builder
# -----------------------------------------------------------------
rl_offset <- offsets_df %>%
  filter(phase == 0) %>%
  select(readlen, p_offset) %>%
  distinct()

build_phase_indep_psites <- function(rpfs_in, anno, rl_offset) {
  orfs <- anno$trspacecds
  ov <- findOverlaps(rpfs_in, orfs, select = "first", ignore.strand = TRUE)
  rpfs2 <- rpfs_in
  rpfs2$orf <- S4Vectors::Rle(names(orfs)[ov])
  rpfs2 <- rpfs2[!is.na(rpfs2$orf)]
  off <- as.data.frame(mcols(rpfs2)[, "readlen", drop = FALSE]) %>%
    left_join(rl_offset, by = "readlen")
  keep <- !is.na(off$p_offset)
  rpfs2 <- rpfs2[keep]
  off   <- off[keep, ]
  ps <- rpfs2 %>% resize(1, "start") %>% shift(off$p_offset)
  # Bounds check
  sl <- seqlengths(ps)
  lens <- sl[as.character(seqnames(ps))]
  ok <- (start(ps) >= 1) & (is.na(lens) | start(ps) <= lens)
  ps[ok]
}

message("Building phase-indep psites (real) ...")
psA_real <- suppressWarnings(
  build_phase_indep_psites(rpfs, chr22_anno, rl_offset)
)
message("Building phase-indep psites (randomized) ...")
psA_rand <- suppressWarnings(
  build_phase_indep_psites(rpfs_rand, chr22_anno, rl_offset)
)

# -----------------------------------------------------------------
# Run ftest_orfs
# -----------------------------------------------------------------
message("ftest_orfs on phase-indep real ...")
ft_real <- suppressMessages(ftest_orfs(psA_real, chr22_anno))
message("ftest_orfs on phase-indep randomized ...")
ft_rand <- suppressMessages(ftest_orfs(psA_rand, chr22_anno))

summarize_ft <- function(label, ft) {
  tested <- ft[!is.na(ft$p.value), ]
  cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  ORFs tested      : %d\n", nrow(tested)))
  cat(sprintf("  median p         : %.4g\n", median(tested$p.value)))
  cat(sprintf("  frac p < 0.05    : %.1f%%\n", 100*mean(tested$p.value < 0.05)))
  cat(sprintf("  frac p < 0.01    : %.1f%%\n", 100*mean(tested$p.value < 0.01)))
  cat("\n")
  tested
}

cat("\n==== Phase-independent psite construction ====\n\n")
summarize_ft("REAL reads (phase-indep)", ft_real)
summarize_ft("RANDOMIZED reads (phase-indep, expected NULL)", ft_rand)

# Frame distributions
orf_start_of <- function(ids) start(chr22_anno$trspacecds[ids])
fr_real <- (start(psA_real) -
            orf_start_of(as.character(mcols(psA_real)$orf))) %% 3
fr_rand <- (start(psA_rand) -
            orf_start_of(as.character(mcols(psA_rand)$orf))) %% 3

cat("\nFrame distribution of phase-indep psites:\n")
cat(sprintf("REAL       : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100*mean(fr_real == 0), 100*mean(fr_real == 1), 100*mean(fr_real == 2)))
cat(sprintf("RANDOMIZED : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100*mean(fr_rand == 0), 100*mean(fr_rand == 1), 100*mean(fr_rand == 2)))
cat("\nExpect REAL: strong f0 bias. Expect RANDOMIZED: ~33/33/33.\n")
