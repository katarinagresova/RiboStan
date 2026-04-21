#!/usr/bin/env Rscript
# Null test: if we randomize read positions (destroying any real
# 3-nt periodicity), does the current pipeline still report
# significant periodicity via ftest_orfs?
#
# If the pipeline is tautological, it will call many/most ORFs
# periodic even on random data. If it is honest, the F-test p-values
# should be roughly uniformly distributed on [0, 1] under the null.

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
})

data(chr22_anno)
data(rpfs)
data(offsets_df)

set.seed(1)

# -----------------------------------------------------------------
# 1. Scramble rpfs: randomize read 5' positions uniformly within
#    the host transcript, keeping readlen unchanged.
# -----------------------------------------------------------------
tx_lens <- seqlengths(rpfs)
rpf_seqs <- as.character(seqnames(rpfs))
rpf_lens <- tx_lens[rpf_seqs]

# Random uniform 5' end within [1, tx_len - readlen + 1]
max_start <- rpf_lens - rpfs$readlen + 1L
max_start <- pmax(max_start, 1L)
rand_start <- pmin(
  pmax(1L, as.integer(round(runif(length(rpfs), 1, max_start)))),
  max_start
)

# Rebuild ranges in the correct order so width>=1 is never violated
rpfs_rand <- GRanges(
  seqnames = seqnames(rpfs),
  ranges   = IRanges(start = rand_start, width = rpfs$readlen),
  strand   = strand(rpfs),
  seqinfo  = seqinfo(rpfs)
)
mcols(rpfs_rand) <- mcols(rpfs)

# -----------------------------------------------------------------
# 2. Run pipeline on real vs randomized
# -----------------------------------------------------------------
message("Computing psites (real) ...")
psites_real <- suppressWarnings(get_psite_gr(rpfs, offsets_df, chr22_anno))
message("Computing psites (randomized positions) ...")
psites_rand <- suppressWarnings(get_psite_gr(rpfs_rand, offsets_df, chr22_anno))

cat(sprintf("Real     psites: %d\n", length(psites_real)))
cat(sprintf("Random   psites: %d\n\n", length(psites_rand)))

# -----------------------------------------------------------------
# 3. Run ftest_orfs on each, compare p-value distributions
# -----------------------------------------------------------------
message("Running ftest_orfs (real) ...")
ft_real <- suppressMessages(ftest_orfs(psites_real, chr22_anno, n_cores = 1))
message("Running ftest_orfs (randomized) ...")
ft_rand <- suppressMessages(ftest_orfs(psites_rand, chr22_anno, n_cores = 1))

# Focus only on ORFs that actually got tested (non-NA p-values)
summarize_ft <- function(label, ft) {
  tested <- ft[!is.na(ft$p.value), ]
  cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  ORFs tested              : %d\n", nrow(tested)))
  cat(sprintf("  median p                 : %.4g\n", median(tested$p.value)))
  cat(sprintf("  frac p < 0.05            : %.1f%%\n", 100*mean(tested$p.value < 0.05)))
  cat(sprintf("  frac p < 0.01            : %.1f%%\n", 100*mean(tested$p.value < 0.01)))
  cat(sprintf("  frac p < 0.001           : %.1f%%\n", 100*mean(tested$p.value < 0.001)))
  cat("\n")
  tested
}

cat("\n==== ftest_orfs results ====\n\n")
t_real <- summarize_ft("REAL data (pipeline)", ft_real)
t_rand <- summarize_ft("RANDOMIZED reads (pipeline)", ft_rand)

# Split by uORF vs CDS
is_uorf <- chr22_anno$uORF
t_real$is_uorf <- is_uorf[t_real$orf_id]
t_rand$is_uorf <- is_uorf[t_rand$orf_id]

cat("\n--- Breakdown by uORF vs CDS (randomized data) ---\n")
rand_cds  <- t_rand %>% filter(!is_uorf)
rand_uorf <- t_rand %>% filter(is_uorf)
cat(sprintf("CDS  (null): %d tested, frac p<0.05 = %.1f%%\n",
            nrow(rand_cds),  100*mean(rand_cds$p.value  < 0.05, na.rm=TRUE)))
cat(sprintf("uORF (null): %d tested, frac p<0.05 = %.1f%%\n",
            nrow(rand_uorf), 100*mean(rand_uorf$p.value < 0.05, na.rm=TRUE)))

# -----------------------------------------------------------------
# 4. Distribution plot of p-values
# -----------------------------------------------------------------
plot_df <- bind_rows(
  t_real %>% mutate(cond = "Real reads"),
  t_rand %>% mutate(cond = "Randomized reads (null)")
)
p <- ggplot(plot_df, aes(p.value)) +
  geom_histogram(boundary = 0, binwidth = 0.02, fill = "steelblue", colour = "white") +
  facet_grid(cond ~ ifelse(is_uorf, "uORFs", "CDS")) +
  labs(
    title = "F-test p-value distribution: real reads vs positionally randomized reads",
    subtitle = "If pipeline is calibrated, randomized data should give UNIFORM p-values (dashed)",
    x = "p.value", y = "count"
  ) +
  theme_bw(base_size = 12)
ggsave("scripts/null_periodicity_pvals.png", p, width = 10, height = 6, dpi = 150)
message("\nPlot saved: scripts/null_periodicity_pvals.png")

# -----------------------------------------------------------------
# 5. Also check the *raw* psite frame distribution under randomization
# -----------------------------------------------------------------
orf_start_of <- function(ids) start(chr22_anno$trspacecds[ids])
fr_real <- (start(psites_real) -
            orf_start_of(as.character(mcols(psites_real)$orf))) %% 3
fr_rand <- (start(psites_rand) -
            orf_start_of(as.character(mcols(psites_rand)$orf))) %% 3

cat("\n==== Frame distribution of pipeline psites ====\n")
cat(sprintf("REAL       : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100*mean(fr_real == 0), 100*mean(fr_real == 1), 100*mean(fr_real == 2)))
cat(sprintf("RANDOMIZED : f0=%.1f%% f1=%.1f%% f2=%.1f%%\n",
            100*mean(fr_rand == 0), 100*mean(fr_rand == 1), 100*mean(fr_rand == 2)))
cat("If pipeline is tautological, RANDOMIZED should still show frame bias.\n")
