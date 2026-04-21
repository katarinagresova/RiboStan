#!/usr/bin/env Rscript
# Sanity check: after applying all three fixes, does the full
# real pipeline give a calibrated null?

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(dplyr)
})

data(chr22_anno)
data(rpfs)
data(offsets_df)

set.seed(1)

# Randomize read positions (null)
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

message("Computing psites (real)  ...")
ps_real <- suppressWarnings(get_psite_gr(rpfs, offsets_df, chr22_anno))
message("Computing psites (null)  ...")
ps_null <- suppressWarnings(get_psite_gr(rpfs_rand, offsets_df, chr22_anno))

message("Running ftest_orfs (real) ...")
ft_real <- suppressMessages(ftest_orfs(ps_real, chr22_anno))
message("Running ftest_orfs (null) ...")
ft_null <- suppressMessages(ftest_orfs(ps_null, chr22_anno))

summ <- function(lbl, df) {
  p <- df$p.value[!is.na(df$p.value)]
  cat(sprintf("%-35s  n=%4d  median p=%.4g  p<0.05=%.1f%%  p<0.01=%.1f%%\n",
              lbl, length(p), median(p),
              100 * mean(p < 0.05), 100 * mean(p < 0.01)))
}

cat("\n==== After fixes: get_psite_gr + ftest_orfs ====\n")
summ("Real reads",       ft_real)
summ("Randomized reads", ft_null)
cat("\nExpected for calibrated null: median p ~0.5, p<0.05 ~5%, p<0.01 ~1%\n")

# Frame distribution sanity check
orf_st <- function(ids) start(chr22_anno$trspacecds[ids])
fr_r <- (start(ps_real) - orf_st(as.character(mcols(ps_real)$orf))) %% 3
fr_n <- (start(ps_null) - orf_st(as.character(mcols(ps_null)$orf))) %% 3
cat(sprintf("\nFrame distribution of psites:\n"))
cat(sprintf("Real       : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100 * mean(fr_r == 0), 100 * mean(fr_r == 1), 100 * mean(fr_r == 2)))
cat(sprintf("Randomized : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100 * mean(fr_n == 0), 100 * mean(fr_n == 1), 100 * mean(fr_n == 2)))
cat("Expected for null: 33/33/33.\n")

# ftest_orfs columns sanity
cat("\nftest_orfs columns: ", paste(colnames(ft_real), collapse = ", "), "\n")
cat("Head of real ftest_orfs output:\n")
print(head(ft_real %>% filter(!is.na(p.value))))
