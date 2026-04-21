#!/usr/bin/env Rscript
# Diagnostic: pick the single uORF with the most P-sites and verify
# that the coverage profile shows genuine 3-nt periodicity in the
# correct frame relative to the uORF AUG.

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
})

message("Loading data ...")
data(chr22_anno)
data(rpfs)
data(offsets_df)

message("Computing P-sites ...")
psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
message(sprintf("Total P-sites: %d", length(psites)))

# --------------------------------------------------------------------------
# 1. Identify uORF IDs from the anno
# --------------------------------------------------------------------------
uorf_ids <- names(chr22_anno$uORF)[chr22_anno$uORF]
message(sprintf("uORF IDs in anno: %d", length(uorf_ids)))

# --------------------------------------------------------------------------
# 2. Count psites per uORF
#    psites$orf is the ORF ID; just tally
# --------------------------------------------------------------------------
psite_orf <- as.character(GenomicRanges::mcols(psites)$orf)
ps_on_uorfs <- psites[psite_orf %in% uorf_ids]
message(sprintf("P-sites assigned to a uORF: %d", length(ps_on_uorfs)))

counts <- sort(table(as.character(GenomicRanges::mcols(ps_on_uorfs)$orf)),
               decreasing = TRUE)
message("Top 10 uORFs by P-site count:")
print(head(counts, 10))

if (length(counts) == 0) stop("No P-sites on any uORF.")

top_id   <- names(counts)[1]
n_top    <- as.integer(counts[1])
top_orf  <- chr22_anno$trspacecds[top_id]   # GRanges in tx-space
message(sprintf("\nTop uORF: %s  (%d P-sites, length %d nt)",
                top_id, n_top, width(top_orf)))

# --------------------------------------------------------------------------
# 3. Get psites for this uORF; compute position & frame relative to ORF AUG
# --------------------------------------------------------------------------
ps_top <- ps_on_uorfs[as.character(GenomicRanges::mcols(ps_on_uorfs)$orf) == top_id]

orf_start <- start(top_orf)
ps_start  <- start(ps_top)
pos_in_orf <- ps_start - orf_start + 1L      # 1 = first nt of AUG
frame      <- (ps_start - orf_start) %% 3L   # 0 = in-frame

message("\nPosition distribution (transcript space, relative to ORF AUG):")
print(sort(table(pos_in_orf)))
message("\nFrame distribution (0 = in-frame with AUG):")
print(table(frame))
message(sprintf("Fraction frame 0: %.1f%%", 100 * mean(frame == 0)))

# --------------------------------------------------------------------------
# 4. Detailed per-read table: readlen, phase, p_offset, shft applied,
#    raw psite frame
# --------------------------------------------------------------------------
ps_df <- as.data.frame(GenomicRanges::mcols(ps_top))
ps_df$pos_in_orf <- pos_in_orf
ps_df$frame      <- frame

# Reconstruct the phaseshift that was applied globally
rl_phs_ <- as.data.frame(GenomicRanges::mcols(psites))[, c("phase","readlen")] %>%
  dplyr::count(phase, readlen) %>%
  dplyr::group_by(readlen) %>%
  dplyr::mutate(shft = rank(-n) - 1) %>%
  dplyr::ungroup()

message("\nGlobal phaseshift table (rank of each phase per readlen):")
print(rl_phs_ %>% dplyr::arrange(readlen, phase) %>% as.data.frame())

ps_df <- ps_df %>%
  left_join(rl_phs_, by = c("readlen", "phase")) %>%
  mutate(
    net_shift        = p_offset - phase + shft,
    # what frame does the offset+phaseshift put us in, relative to phase-0 read?
    frame_from_offsets = (p_offset + shft) %% 3
  )

message("\nPer-readlen/phase summary for top uORF psites:")
ps_df %>%
  dplyr::count(readlen, phase, p_offset, shft, net_shift, frame_from_offsets) %>%
  dplyr::arrange(readlen, phase) %>%
  as.data.frame() %>%
  print()

message("\nFrame distribution by read phase:")
print(with(ps_df, table(phase = factor(phase), frame = factor(frame))))

# Direct position check: show raw transcript positions for a few reads of each phase
message("\nDirect position check (orf_start = ", orf_start, "):")
for (ph in 0:2) {
  idx <- which(ps_df$phase == ph)
  if (length(idx) == 0) next
  take <- head(idx, 5)
  ps_tx_pos <- start(ps_top[take])
  cat(sprintf("  phase=%d | transcript positions: %s\n",
              ph, paste(ps_tx_pos, collapse=", ")))
  cat(sprintf("           pos-orf_start: %s\n",
              paste(ps_tx_pos - orf_start, collapse=", ")))
  cat(sprintf("           frame (tx-orf_start)%%3: %s\n",
              paste((ps_tx_pos - orf_start) %% 3, collapse=", ")))
}

# Also: for phase=1 reads, what offset+phaseshift was applied?
# Verify by back-computing the original read start
message("\nBack-compute check for phase=1 reads (first 5):")
ph1 <- ps_df[ps_df$phase == 1, ][1:min(5, sum(ps_df$phase == 1)), ]
for (i in seq_len(nrow(ph1))) {
  psite_pos <- start(ps_top[rownames(ph1)[i] %in% names(ps_top)])
  cat(sprintf("  net_shift=%d  → psite_pos - orf_start = %d  (frame %d)\n",
              ph1$net_shift[i],
              ph1$pos_in_orf[i] - 1L,
              (ph1$pos_in_orf[i] - 1L) %% 3))
}

# --------------------------------------------------------------------------
# 5. Position × frame plot (coverage barplot)
# --------------------------------------------------------------------------
orf_len <- width(top_orf)
cov_df <- ps_df %>%
  dplyr::filter(pos_in_orf >= 1, pos_in_orf <= orf_len) %>%
  dplyr::count(pos_in_orf, frame = factor(frame, levels = 0:2)) %>%
  tidyr::complete(pos_in_orf = 1:orf_len, frame, fill = list(n = 0L))

p_cov <- ggplot(cov_df, aes(pos_in_orf, n, fill = frame)) +
  geom_col() +
  scale_fill_manual(values = c("0" = "#2166ac", "1" = "#d6604d", "2" = "#74c476")) +
  labs(
    title    = sprintf("P-site coverage: uORF %s", top_id),
    subtitle = sprintf("%d nt, %d P-sites assigned", orf_len, n_top),
    x = "Position relative to uORF AUG (nt, 1-based)",
    y = "P-site count",
    fill = "Frame"
  ) +
  theme_bw(base_size = 12)

out <- "scripts/debug_uorf_coverage.png"
ggsave(out, p_cov, width = 12, height = 4, dpi = 150)
message("Coverage plot: ", out)

# --------------------------------------------------------------------------
# 6. Compare annotated CDS frame distribution
# --------------------------------------------------------------------------
cds_ids  <- names(chr22_anno$uORF)[!chr22_anno$uORF]
ps_cds   <- psites[as.character(GenomicRanges::mcols(psites)$orf) %in% cds_ids]
cds_orf  <- chr22_anno$trspacecds[as.character(GenomicRanges::mcols(ps_cds)$orf)]
cds_frame <- (start(ps_cds) - start(cds_orf)) %% 3L
message("\nAnnotated CDS frame distribution:")
print(table(cds_frame))
message(sprintf("Fraction frame 0 (CDS): %.1f%%", 100 * mean(cds_frame == 0)))
message(sprintf("Fraction frame 0 (uORF top): %.1f%%", 100 * mean(frame == 0)))

message("\nDone.")
