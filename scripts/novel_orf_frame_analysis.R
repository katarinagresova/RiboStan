#!/usr/bin/env Rscript
# =============================================================================
# Novel ORF Frame Analysis
# =============================================================================
# Demonstrates translation evidence for ORFs outside annotated CDS by showing
# 3-nt P-site periodicity in three ORF categories:
#   1. Main ORFs          (annotated CDS — positive control)
#   2. Periodic uORFs     (subset passing multitaper F-test, p < 0.05)
#   3. Non-periodic uORFs (rest of the uORFs)
#
# Runs entirely on shipped package data — no downloads required.
#
# Background on the x-axis:
#   After P-site assignment, Ribostan applies a per-read-length phaseshift so
#   that reads of different lengths land on a common codon-relative frame.
#   Each (readlen, phase) bin is assigned a shift in {0, 1, 2} equal to the
#   rank of its count among the three phases for that read length. Reads of
#   the dominant phase (rank 0) stay at frame 0 of the ORF; other phases are
#   shifted to frames 1 and 2. This is why the x-axis here is called
#   "read-length adjusted frame": reads of different lengths but the same
#   underlying codon position appear in the same frame bin.
# =============================================================================

suppressPackageStartupMessages({
  library(Ribostan)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

# -----------------------------------------------------------------------------
# 1. Load shipped example data
# -----------------------------------------------------------------------------
message("Loading package data...")
data(chr22_anno)   # annotated CDS + ORFik-called uORFs
data(rpfs)         # ribosome footprint GRanges (transcript coords, chr22)
data(offsets_df)   # P-site offsets calibrated from CDS-inclusion

message(sprintf(
  "Annotation: %d annotated ORFs + %d uORFs on %d transcripts",
  sum(!chr22_anno$uORF), sum(chr22_anno$uORF), length(chr22_anno$exonsgrl)
))

# -----------------------------------------------------------------------------
# 2. Assign P-sites
# -----------------------------------------------------------------------------
message("Assigning P-sites...")
psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
message(sprintf("  %d P-sites across %d ORFs",
  length(psites), length(unique(as.character(psites$orf)))))

uorf_ids <- names(chr22_anno$uORF[chr22_anno$uORF])
cds_ids  <- names(chr22_anno$uORF[!chr22_anno$uORF])

# -----------------------------------------------------------------------------
# 3. Periodicity test on uORFs
# -----------------------------------------------------------------------------
message("Running multitaper F-tests on uORFs...")
filt_anno <- periodicity_filter_uORFs(
  psites, chr22_anno, remove = FALSE, n_cores = 1
)
ptest_df <- mcols(filt_anno$trspacecds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("orf_id") %>%
  filter(orf_id %in% uorf_ids)

periodic_ids     <- ptest_df %>%
  filter(!is.na(p.value), p.value <  0.05) %>% pull(orf_id)
non_periodic_ids <- ptest_df %>%
  filter(!is.na(p.value), p.value >= 0.05) %>% pull(orf_id)

message(sprintf(
  "  %d periodic uORFs (p < 0.05) / %d non-periodic uORFs / %d uORFs total tested",
  length(periodic_ids), length(non_periodic_ids),
  length(periodic_ids) + length(non_periodic_ids)
))

# -----------------------------------------------------------------------------
# 4. Frame computation: ORF-relative position of each (post-phaseshift) P-site
# -----------------------------------------------------------------------------
# The same logic used inside ftest_orfs: map psites back to transcript-space
# ORF coordinates, then (pos - 1) %% 3 is the 0-based frame within the ORF.
compute_frame <- function(ps, anno, orf_ids) {
  covered <- intersect(as.character(unique(ps$orf)), orf_ids)
  if (length(covered) == 0) return(integer(0))
  orfs_gr <- anno$trspacecds[covered]
  mapped  <- GenomicFeatures::mapToTranscripts(ps, orfs_gr, ignore.strand = TRUE)
  keep    <- as.character(ps$orf[mapped$xHits]) ==
             names(orfs_gr)[mapped$transcriptsHits]
  mapped  <- mapped[keep]
  (GenomicRanges::start(mapped) - 1L) %% 3L
}

make_frame_df <- function(ps, anno, orf_ids, label) {
  frame <- compute_frame(ps, anno, orf_ids)
  if (length(frame) == 0) {
    return(data.frame(frame = factor(integer(0), levels = 0:2),
                      label = character(0), stringsAsFactors = FALSE))
  }
  data.frame(
    frame = factor(frame, levels = 0:2),
    label = label,
    stringsAsFactors = FALSE
  )
}

ps_cds      <- psites[as.character(psites$orf) %in% cds_ids]
ps_periodic <- psites[as.character(psites$orf) %in% periodic_ids]
ps_nonperio <- psites[as.character(psites$orf) %in% non_periodic_ids]

frame_df <- bind_rows(
  make_frame_df(ps_cds,      chr22_anno, cds_ids,          "Main ORFs"),
  make_frame_df(ps_periodic, chr22_anno, periodic_ids,     "Periodic uORFs"),
  make_frame_df(ps_nonperio, chr22_anno, non_periodic_ids, "Non-periodic uORFs")
) %>%
  mutate(label = factor(
    label, levels = c("Main ORFs", "Periodic uORFs", "Non-periodic uORFs")
  )) %>%
  count(label, frame, name = "count", .drop = FALSE) %>%
  group_by(label) %>%
  mutate(pct = if (sum(count) > 0) count / sum(count) * 100 else count * 0) %>%
  ungroup()

total_counts <- frame_df %>%
  group_by(label) %>%
  summarise(total = sum(count), n_orfs = NA_integer_, .groups = "drop")
total_counts$n_orfs <- c(
  length(intersect(cds_ids, as.character(unique(psites$orf)))),
  length(intersect(periodic_ids, as.character(unique(psites$orf)))),
  length(intersect(non_periodic_ids, as.character(unique(psites$orf))))
)

message("\nP-site counts per category:")
print(as.data.frame(total_counts))

# -----------------------------------------------------------------------------
# 5. Plot — three panels side-by-side
# -----------------------------------------------------------------------------
fills <- c(
  "Main ORFs"          = "#2166AC",
  "Periodic uORFs"     = "#2CA25F",
  "Non-periodic uORFs" = "#D6604D"
)

facet_labels <- total_counts %>%
  mutate(facet_label = sprintf("%s\n(n = %s P-sites, %d ORFs)",
    label, format(total, big.mark = ","), n_orfs)) %>%
  { setNames(.$facet_label, as.character(.$label)) }

frame_df <- frame_df %>%
  mutate(label_n = factor(
    facet_labels[as.character(label)],
    levels = facet_labels[levels(label)]
  ))

p <- ggplot(frame_df, aes(x = frame, y = pct, fill = label)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = 33.3, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  annotate("text", x = 0.55, y = 36, label = "random",
           size = 3, colour = "grey50", hjust = 0) +
  facet_wrap(~label_n, nrow = 1) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 33.3, 50, 75, 100),
    labels = c("0%", "33%", "50%", "75%", "100%")
  ) +
  labs(
    title    = "3-nt P-site periodicity: main ORFs vs uORFs (chr22, hg38)",
    subtitle = paste0(
      "x-axis: read-length adjusted frame within the ORF.\n",
      "get_psite_gr applies a per-(readlen, phase) shift in {0, 1, 2} equal to ",
      "the rank of that bin's\n",
      "count for that read length, so reads of different lengths sharing the ",
      "same codon land in the\n",
      "same frame bin. Frame 0 = in-frame with the ORF start codon."
    ),
    x = "Read-length adjusted frame (0 = in-frame with ORF AUG)",
    y = "% of P-sites in frame"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background   = element_rect(fill = "#EEF4FB"),
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.subtitle      = element_text(size = 9, colour = "grey40",
                                      lineheight = 1.3)
  )

outfile <- "novel_orf_frame_periodicity.png"
ggsave(outfile, p, width = 10, height = 5, dpi = 150)
message(sprintf("\nPlot saved to: %s",
                normalizePath(outfile, mustWork = FALSE)))

cat("\n=== Frame distribution summary (% of P-sites per frame) ===\n")
frame_df %>%
  mutate(pct = round(pct, 1)) %>%
  select(label, frame, pct) %>%
  tidyr::pivot_wider(names_from = frame, values_from = pct,
                     names_prefix = "frame_") %>%
  as.data.frame() %>%
  print()
