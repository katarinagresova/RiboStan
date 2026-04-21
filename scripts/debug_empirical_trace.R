#!/usr/bin/env Rscript
# Empirical trace of the phase/offset/phaseshift logic.
# Goal: verify by hand how a known read of known length/phase lands.

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(dplyr)
})

data(chr22_anno)
data(rpfs)
data(offsets_df)

message("==== 1. Are offsets strictly multiples of 3? ====")
print(offsets_df %>% mutate(mod3 = p_offset %% 3) %>% as.data.frame())
cat(sprintf("All offsets %% 3 == 0: %s\n\n", all((offsets_df$p_offset %% 3) == 0)))

message("==== 2. Does merge() reorder phaseshifts? ====")
# Reproduce the merge inside get_psite_gr exactly to see
set.seed(1)
n <- 20
demo <- data.frame(
  phase   = sample(0:2, n, replace = TRUE),
  readlen = sample(26:29, n, replace = TRUE)
)
demo$idx <- seq_len(n)
shiftdf <- expand.grid(phase = 0:2, readlen = 26:29) %>%
  mutate(shft = (phase * 1 + readlen) %% 3)  # dummy shft values
m <- merge(demo, shiftdf, all.x = TRUE)
cat("Input row order (phase, readlen, idx):\n")
print(demo)
cat("After merge() — note the row reshuffle:\n")
print(m)
cat(sprintf("Row order preserved? %s\n\n",
            identical(m$idx, demo$idx)))

message("==== 3. Trace one phase=0 and one phase=1 read from real data ====")
# Take one read of each phase that overlaps the top uORF
top_uorf_id <- "TxID:214938_1"
top_orf <- chr22_anno$trspacecds[top_uorf_id]
orf_start <- start(top_orf)
orf_end   <- end(top_orf)
cat(sprintf("uORF %s on transcript %s, tx-range [%d, %d], len=%d\n",
            top_uorf_id, as.character(seqnames(top_orf)),
            orf_start, orf_end, width(top_orf)))

# Find RPFs on the host transcript overlapping the uORF
host_tx <- as.character(seqnames(top_orf))
rpfs_on_host <- rpfs[as.character(seqnames(rpfs)) == host_tx]
rpfs_on_uorf <- rpfs_on_host[
  start(rpfs_on_host) >= orf_start & start(rpfs_on_host) <= orf_end
]
cat(sprintf("RPFs with 5' end inside uORF: %d\n", length(rpfs_on_uorf)))

# compute phase relative to uORF start
rpf_phase <- (start(rpfs_on_uorf) - orf_start) %% 3
cat("Phase distribution of RPFs on uORF:\n")
print(table(phase = rpf_phase, readlen = rpfs_on_uorf$readlen))

message("\n==== 4. Run get_psite_gr and locate the same reads after processing ====")
psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
# psites on this uORF
ps_on_uorf <- psites[as.character(mcols(psites)$orf) == top_uorf_id]
cat(sprintf("P-sites assigned to top uORF: %d\n", length(ps_on_uorf)))
cat("Frame of those P-sites:\n")
print(table((start(ps_on_uorf) - orf_start) %% 3))

message("\n==== 5. Check if phaseshift-merge-bug actually corrupts shifts ====")
# Reconstruct what get_psite_gr does step by step
rpfs2 <- rpfs
orfs <- chr22_anno$trspacecds
ov <- findOverlaps(rpfs2, orfs, select="first", ignore.strand=TRUE)
rpfs2$orf   <- S4Vectors::Rle(names(orfs)[ov])
rpfs2$phase <- (start(rpfs2) - start(orfs)[ov]) %% 3
rpfs2$p_offset <- as.data.frame(mcols(rpfs2)[,c("readlen","phase")]) %>%
  left_join(offsets_df %>% select(readlen, phase, p_offset),
            by=c("readlen","phase")) %>% .$p_offset
rpfs2 <- subset(rpfs2, !is.na(p_offset))
# post initial shift
ps_pre <- rpfs2 %>% resize(1,"start") %>% shift(.$p_offset)
# bounds check against orf end (treat out-of-ORF as OK for this demo)
seq_lens <- seqlengths(ps_pre)
has_len  <- !is.na(seq_lens[as.character(seqnames(ps_pre))])
keep     <- (start(ps_pre) >= 1) &
            (is.na(seq_lens[as.character(seqnames(ps_pre))]) |
             start(ps_pre) <= seq_lens[as.character(seqnames(ps_pre))])
ps_pre <- ps_pre[keep]

# Without the merge bug, what would frames be?
# Do the "correct" phaseshift: per-row lookup
rl_phs <- as.data.frame(mcols(ps_pre)[,c("phase","readlen")]) %>%
  count(phase, readlen) %>%
  group_by(readlen) %>% mutate(shft = rank(-n) - 1) %>% ungroup()
cat("Per-(phase,readlen) shift table:\n")
print(as.data.frame(rl_phs %>% arrange(readlen, phase)))

# Aligned lookup (correct)
lookup <- as.data.frame(mcols(ps_pre)[,c("phase","readlen")]) %>%
  left_join(rl_phs %>% select(phase, readlen, shft), by=c("phase","readlen"))
stopifnot(nrow(lookup) == length(ps_pre))
ps_correct <- ps_pre %>% shift(-.$phase) %>% shift(lookup$shft)

# Buggy version (as currently in get_psite_gr — merge reorders)
buggy_merge <- merge(
  as.data.frame(mcols(ps_pre)[,c("phase","readlen")]),
  as.data.frame(rl_phs %>% select(phase, readlen, shft)),
  all.x = TRUE
)
ps_buggy <- ps_pre %>% shift(-.$phase) %>% shift(buggy_merge$shft)

# Now measure frame distribution on ANNOTATED CDS (where we have good signal)
# against both methods
cds_ids <- names(chr22_anno$uORF)[!chr22_anno$uORF]
in_cds <- as.character(mcols(ps_pre)$orf) %in% cds_ids
orf_starts_per_ps <- start(orfs)[match(as.character(mcols(ps_pre)$orf), names(orfs))]
frame_correct <- (start(ps_correct) - orf_starts_per_ps) %% 3
frame_buggy   <- (start(ps_buggy)   - orf_starts_per_ps) %% 3

cat("\nCDS frame distribution (correct per-row lookup):\n")
print(table(frame_correct[in_cds]))
cat(sprintf("Fraction frame 0: %.1f%%\n", 100 * mean(frame_correct[in_cds] == 0)))

cat("\nCDS frame distribution (buggy merge):\n")
print(table(frame_buggy[in_cds]))
cat(sprintf("Fraction frame 0: %.1f%%\n", 100 * mean(frame_buggy[in_cds] == 0)))

message("\n==== 6. What fraction of shft values differ between correct vs buggy? ====")
cat(sprintf("Psites where shft differs: %d / %d (%.1f%%)\n",
            sum(lookup$shft != buggy_merge$shft),
            nrow(lookup),
            100 * mean(lookup$shft != buggy_merge$shft)))
