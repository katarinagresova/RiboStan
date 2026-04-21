#!/usr/bin/env Rscript
# Demonstrate that the current psite pipeline is tautological for the
# periodicity test, and show what a phase-independent alternative
# looks like on a specific uORF.
#
# Two psite constructions to compare:
#   (A) "Pipeline": current get_psite_gr output
#         - applies phase-dependent offset (all %% 3 = 0)
#         - then shift(-phase) forces P-site to frame 0 of assigned ORF
#         - then shift(shft) redistributes by global rank
#   (B) "Phase-independent": just use per-readlen offset (ignore phase),
#         preserving the raw 3-nt phase signal of the 5' ends.

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(dplyr)
})

data(chr22_anno)
data(rpfs)
data(offsets_df)

# -----------------------------------------------------------------
# (A) Pipeline P-sites
# -----------------------------------------------------------------
psites_A <- get_psite_gr(rpfs, offsets_df, chr22_anno)

# -----------------------------------------------------------------
# (B) Phase-independent P-sites: one offset per read length
#     (use the phase=0 offset, i.e. the dominant one)
# -----------------------------------------------------------------
rl_offset <- offsets_df %>%
  group_by(readlen) %>%
  slice_max(order_by = -phase, n = 1, with_ties = FALSE) %>%  # phase=0 row
  ungroup() %>%
  select(readlen, fixed_offset = p_offset)

# Assign reads to ORFs (same logic as pipeline), but shift by a
# phase-independent per-readlen offset.
orfs <- chr22_anno$trspacecds
rpfs_B <- rpfs
ov <- findOverlaps(rpfs_B, orfs, select="first", ignore.strand=TRUE)
rpfs_B$orf <- S4Vectors::Rle(names(orfs)[ov])
rpfs_B <- rpfs_B[!is.na(rpfs_B$orf)]
# merge fixed offset
mb <- as.data.frame(mcols(rpfs_B)[, "readlen", drop=FALSE]) %>%
  left_join(rl_offset, by="readlen")
rpfs_B <- rpfs_B[!is.na(mb$fixed_offset)]
mb <- mb[!is.na(mb$fixed_offset), ]
psites_B <- rpfs_B %>% resize(1, "start") %>% shift(mb$fixed_offset)

# -----------------------------------------------------------------
# Compare on the top uORF
# -----------------------------------------------------------------
top_id    <- "TxID:214938_1"
top_orf   <- chr22_anno$trspacecds[top_id]
orf_start <- start(top_orf)
orf_end   <- end(top_orf)

get_frames <- function(ps, orf_id) {
  on_orf <- ps[as.character(mcols(ps)$orf) == orf_id]
  (start(on_orf) - orf_start) %% 3
}

fA <- get_frames(psites_A, top_id)
fB <- get_frames(psites_B, top_id)

cat("==== Top uORF =", top_id, "====\n")
cat(sprintf("uORF range [%d, %d], length %d\n\n", orf_start, orf_end, width(top_orf)))

cat("--- (A) Pipeline P-sites (post shift(-phase) + shft) ---\n")
print(table(frame = fA))
cat(sprintf("Fraction frame 0: %.1f%%  (n=%d)\n\n", 100 * mean(fA == 0), length(fA)))

cat("--- (B) Phase-independent P-sites (fixed offset per readlen) ---\n")
print(table(frame = fB))
cat(sprintf("Fraction frame 0: %.1f%%  (n=%d)\n\n", 100 * mean(fB == 0), length(fB)))

# -----------------------------------------------------------------
# Same comparison aggregated across ALL uORFs vs all annotated CDS
# -----------------------------------------------------------------
is_uorf   <- chr22_anno$uORF
uorf_ids  <- names(is_uorf)[is_uorf]
cds_ids   <- names(is_uorf)[!is_uorf]

frame_dist_for_orfs <- function(ps, orf_ids, anno) {
  on <- ps[as.character(mcols(ps)$orf) %in% orf_ids]
  starts <- start(anno$trspacecds[as.character(mcols(on)$orf)])
  (start(on) - starts) %% 3
}

fA_uorf <- frame_dist_for_orfs(psites_A, uorf_ids, chr22_anno)
fB_uorf <- frame_dist_for_orfs(psites_B, uorf_ids, chr22_anno)
fA_cds  <- frame_dist_for_orfs(psites_A, cds_ids,  chr22_anno)
fB_cds  <- frame_dist_for_orfs(psites_B, cds_ids,  chr22_anno)

cat("==== Aggregated frame distributions ====\n")
show_dist <- function(label, f) {
  t <- table(factor(f, levels = 0:2))
  pct <- 100 * prop.table(t)
  cat(sprintf("  %s  n=%d  frame0=%.1f%%  frame1=%.1f%%  frame2=%.1f%%\n",
              label, length(f), pct[1], pct[2], pct[3]))
}
cat("(A) Pipeline:\n")
show_dist("all CDS ", fA_cds)
show_dist("all uORFs", fA_uorf)
cat("(B) Phase-independent:\n")
show_dist("all CDS ", fB_cds)
show_dist("all uORFs", fB_uorf)

# -----------------------------------------------------------------
# Null test: random transcript windows that are NOT ORFs — would
# the pipeline call these "periodic" too?
# -----------------------------------------------------------------
cat("\n==== Shuffled-ORF null: assign each psite to a wrong ORF ====\n")
# For each pipeline psite, compute (start - random_other_orf_start) %% 3
set.seed(42)
ps <- psites_A
orfs_vec <- chr22_anno$trspacecds
# sample a random orf_start for each psite
rand_ids  <- sample(names(orfs_vec), length(ps), replace = TRUE)
rand_sts  <- start(orfs_vec[rand_ids])
fA_null   <- (start(ps) - rand_sts) %% 3
show_dist("pipeline, shuffled ORF", fA_null)
# Same with B
ps2 <- psites_B
rand_ids2 <- sample(names(orfs_vec), length(ps2), replace = TRUE)
rand_sts2 <- start(orfs_vec[rand_ids2])
fB_null   <- (start(ps2) - rand_sts2) %% 3
show_dist("phase-indep, shuffled ORF", fB_null)

cat("\nIf pipeline is tautological, pipeline frame0% should be ~33% on shuffled ORFs\n")
cat("(it acts on global read phases, not per-ORF evidence).\n")
cat("If phase-indep preserves raw signal, shuffled-ORF should also be ~33%.\n")
