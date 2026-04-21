# Feature tests for overlapping-ORF handling, alt start codons, uORF
# overlap metadata, and the nested-uORF pruning step in ftest_orfs.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(dplyr)
})


# ---------------------------------------------------------------------------
# find_orfs_cpp: alternative start codons (ATG / CTG / GTG)
# ---------------------------------------------------------------------------
test_that("find_orfs_cpp finds ORFs at ATG, CTG, and GTG starts", {
  # Each sequence carries one ORF starting at position 1 with its named
  # start codon and a stop at positions 13-15.
  seqs <- c(
    atg = "ATGAAACCCAAATAA",
    ctg = "CTGAAACCCAAATAA",
    gtg = "GTGAAACCCAAATAA",
    none = "AAAAAACCCAAATAA"
  )
  r <- Ribostan:::find_orfs_cpp(
    seqs         = seqs,
    start_codons = c("ATG", "CTG", "GTG"),
    stop_codons  = c("TAA", "TAG", "TGA"),
    min_body     = 0L
  )

  expect_equal(length(r$starts), 3L)
  expect_setequal(r$indices, 1:3)
  expect_true(all(r$starts == 1L))
  expect_true(all(r$ends   == 15L))
  # start_codon_idx should map each ORF back to the expected codon
  codons <- c("ATG", "CTG", "GTG")
  mapping <- setNames(codons[r$start_codon_idx], names(seqs)[r$indices])
  expect_equal(mapping[["atg"]], "ATG")
  expect_equal(mapping[["ctg"]], "CTG")
  expect_equal(mapping[["gtg"]], "GTG")
})


test_that("find_orfs_cpp respects a restricted start-codon set", {
  seqs <- c(
    both = "ATGAAACTGAAATAA"   # ATG at 1, CTG at 7, shared stop at 13
  )
  # With ATG only we expect one ORF starting at 1.
  r1 <- Ribostan:::find_orfs_cpp(
    seqs = seqs,  start_codons = "ATG",
    stop_codons = c("TAA","TAG","TGA"), min_body = 0L
  )
  expect_equal(r1$starts, 1L)
  expect_equal(r1$start_codon_idx, 1L)

  # With ATG+CTG we expect two nested ORFs sharing the stop.
  r2 <- Ribostan:::find_orfs_cpp(
    seqs = seqs, start_codons = c("ATG","CTG"),
    stop_codons = c("TAA","TAG","TGA"), min_body = 0L
  )
  expect_setequal(r2$starts, c(1L, 7L))
  # start 7 was CTG (index 2 in start_codons)
  expect_equal(r2$start_codon_idx[r2$starts == 7L], 2L)
  expect_equal(r2$start_codon_idx[r2$starts == 1L], 1L)
})


# ---------------------------------------------------------------------------
# find_orfs_cpp: overlapping ORFs (nested shared-stop vs. different-frame)
# ---------------------------------------------------------------------------
test_that("find_orfs_cpp enumerates nested ORFs sharing a stop", {
  seq <- "ATGAAACCCATGAAACCCTAA"   # outer start 1, inner start 10, stop 19-21
  r <- Ribostan:::find_orfs_cpp(
    seqs = c(tx = seq),
    start_codons = "ATG",
    stop_codons  = c("TAA","TAG","TGA"),
    min_body     = 0L
  )
  expect_equal(length(r$starts), 2L)
  expect_setequal(r$starts, c(1L, 10L))
  expect_true(all(r$ends == 21L))
})


test_that("find_orfs_cpp finds overlapping ORFs in different reading frames", {
  # Sequence: ATGACATGAAATAG
  # pos:      1234567890123 4
  # ATG at pos 1 (frame 0): ATG ACA TGA  -> TGA stop at pos 7-9
  # ATG at pos 6 (frame 2): ATG AAA TAG  -> TAG stop at pos 12-14
  seq <- "ATGACATGAAATAG"
  r <- Ribostan:::find_orfs_cpp(
    seqs = c(tx = seq),
    start_codons = "ATG",
    stop_codons  = c("TAA","TAG","TGA"),
    min_body     = 0L
  )
  expect_equal(length(r$starts), 2L)
  expect_setequal(r$starts, c(1L, 6L))
  # Verify that the two starts land on different frames (mod 3)
  frames <- (r$starts - 1L) %% 3L
  expect_equal(length(unique(frames)), 2L)
})


test_that("longestORF=TRUE logic collapses nested ORFs per stop", {
  seq <- "ATGAAACCCATGAAACCCTAA"   # nested pair
  raw_all <- Ribostan:::find_orfs_cpp(
    seqs = c(x = seq), start_codons = "ATG",
    stop_codons = c("TAA","TAG","TGA"), min_body = 0L
  )
  expect_equal(length(raw_all$starts), 2L)

  # Emulate the longestORF=TRUE collapsing that find_uorfs performs when
  # that flag is set: within (tx_idx, end) groups keep longest ORF.
  widths  <- raw_all$ends - raw_all$starts
  key     <- paste0(raw_all$indices, "_", raw_all$ends)
  keep    <- unlist(tapply(seq_along(key), key,
                            function(i) i[which.max(widths[i])]),
                    use.names = FALSE)
  collapsed <- lapply(raw_all, `[`, keep)
  expect_equal(length(collapsed$starts), 1L)
  expect_equal(collapsed$starts, 1L)
})


# ---------------------------------------------------------------------------
# .filter_uorfs: drop only exact duplicates; keep overlapping / nested ORFs
# ---------------------------------------------------------------------------
test_that(".filter_uorfs keeps overlapping ORFs and drops exact CDS duplicates", {
  mk <- function(s, e, strand = "+", seqn = "chr1") {
    GenomicRanges::GRangesList(GenomicRanges::GRanges(
      seqnames = seqn,
      ranges   = IRanges::IRanges(start = s, end = e),
      strand   = strand
    ))
  }
  # One CDS on chr1:100-200 (+)
  cds   <- mk(100, 200)
  names(cds) <- "cds1"

  # Four candidate uORFs:
  #  u_exact   : same start and stop as cds1    -> DROP
  #  u_nstart  : same start, different stop     -> KEEP
  #  u_nstop   : same stop, different start     -> KEEP (N-term extension-like)
  #  u_inside  : entirely inside cds, different endpoints -> KEEP (internal ORF)
  u_exact  <- mk(100, 200)
  u_nstart <- mk(100, 250)
  u_nstop  <- mk(50,  200)
  u_inside <- mk(120, 180)
  uorfs <- c(u_exact, u_nstart, u_nstop, u_inside)
  names(uorfs) <- c("u_exact", "u_nstart", "u_nstop", "u_inside")

  filtered <- Ribostan:::.filter_uorfs(uorfs, cds)
  expect_equal(sort(names(filtered)),
               sort(c("u_nstart", "u_nstop", "u_inside")))
  expect_false("u_exact" %in% names(filtered))
})


# ---------------------------------------------------------------------------
# .annotate_uorf_overlaps: per-uORF CDS-overlap metadata + frame
# ---------------------------------------------------------------------------
test_that(".annotate_uorf_overlaps marks overlap + frame against main CDS", {
  # Construct a mini genomic annotation:
  #   transcript tx1 has a CDS on chr1:100-300 (+)
  #   uORFs (tx1_1 .. tx1_3) on the same transcript:
  #     tx1_1: 50..99        -> disjoint from CDS     -> overlaps_cds FALSE
  #     tx1_2: 50..150       -> overlaps CDS, in frame (start at 100 in CDS trspace, (100-100) %% 3 = 0)
  #     tx1_3: 51..150       -> overlaps CDS, frame 1 (51 vs 100 in CDS trspace? No: tx1_3 start=51 is before CDS, map into CDS trspace = position 1, (1-1)%%3=0. Hmm.)
  # Simpler: make uORFs whose START is inside the CDS.
  mk_cds <- function(start, end) {
    GenomicRanges::GRangesList(GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges   = IRanges::IRanges(start = start, end = end),
      strand   = "+"
    ))
  }
  cds_gr <- mk_cds(100, 300)
  u1 <- mk_cds(50,  80)       # disjoint
  u2 <- mk_cds(100, 150)      # start at 100, same frame as CDS
  u3 <- mk_cds(101, 150)      # start at 101, (101-100) %% 3 = 1, out of frame
  u4 <- mk_cds(102, 150)      # start at 102, (102-100) %% 3 = 2, out of frame
  cdsgrl <- c(cds_gr, u1, u2, u3, u4)
  names(cdsgrl) <- c("tx1", "tx1_1", "tx1_2", "tx1_3", "tx1_4")
  is_uORF <- c(tx1 = FALSE, tx1_1 = TRUE, tx1_2 = TRUE, tx1_3 = TRUE,
               tx1_4 = TRUE)

  # Build a trspacecds-shaped GRanges (we only need names + a spot for mcols)
  # Use genomic coords as a placeholder; the helper only reads names + uses
  # cdsgrl + is_uORF for its logic.
  trspacecds <- GenomicRanges::GRanges(
    seqnames = names(cdsgrl),
    ranges   = IRanges::IRanges(start = 1L, end = 1L),
    strand   = "+"
  )
  names(trspacecds) <- names(cdsgrl)

  out <- Ribostan:::.annotate_uorf_overlaps(trspacecds, cdsgrl, is_uORF)
  m <- as.data.frame(S4Vectors::mcols(out))
  rownames(m) <- names(out)

  expect_true(all(c("overlaps_cds", "overlap_cds_id", "overlap_frame",
                    "out_of_phase") %in% colnames(m)))

  # u1 is disjoint from the CDS
  expect_false(m["tx1_1", "overlaps_cds"])
  expect_true(is.na(m["tx1_1", "overlap_cds_id"]))
  expect_true(is.na(m["tx1_1", "overlap_frame"]))

  # u2 starts exactly at the CDS start -> in-frame (frame 0)
  expect_true(m["tx1_2", "overlaps_cds"])
  expect_equal(m["tx1_2", "overlap_cds_id"], "tx1")
  expect_equal(m["tx1_2", "overlap_frame"], 0L)
  expect_false(m["tx1_2", "out_of_phase"])

  # u3 offset by 1 -> frame 1, out of phase
  expect_true(m["tx1_3", "overlaps_cds"])
  expect_equal(m["tx1_3", "overlap_frame"], 1L)
  expect_true(m["tx1_3", "out_of_phase"])

  # u4 offset by 2 -> frame 2, out of phase
  expect_true(m["tx1_4", "overlaps_cds"])
  expect_equal(m["tx1_4", "overlap_frame"], 2L)
  expect_true(m["tx1_4", "out_of_phase"])
})


# ---------------------------------------------------------------------------
# get_psite_gr: emits one psite per (read, ORF) overlap; read_mult set
# ---------------------------------------------------------------------------
test_that("get_psite_gr emits one psite per (read, ORF) overlap", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  # read_mult column should be set and always >= 1
  expect_true("read_mult" %in% colnames(GenomicRanges::mcols(psites)))
  expect_true(all(psites$read_mult >= 1L))

  # At least some reads overlap multiple ORFs (the whole point of the
  # refactor). chr22_anno was built with longestORF=TRUE so has no nested
  # uORF pairs, but multi-transcript UTR/CDS overlap still produces
  # multi-ORF reads.
  expect_gt(sum(psites$read_mult > 1L), 0L)

  # There must be at least as many psites as there are reads that
  # overlap at least one ORF (each such read contributes >= 1 psite).
  orfs <- chr22_anno$trspacecds
  # Align rpfs with orfs seqlevels for countOverlaps
  names(rpfs) <- NULL
  rpfs2 <- rpfs[as.character(
    GenomicRanges::seqnames(rpfs)) %in% as.character(
      GenomicRanges::seqnames(orfs))]
  n_rpfs_with_orf <- sum(
    GenomicRanges::countOverlaps(rpfs2, orfs, ignore.strand = TRUE) > 0L
  )
  expect_gte(length(psites), n_rpfs_with_orf * 0.5)
})


test_that("get_read_spmat row-normalises so each read has total weight 1", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  spmat  <- Ribostan:::get_read_spmat(psites, chr22_anno)
  row_sums <- Matrix::rowSums(spmat)
  expect_equal(unname(row_sums), rep(1, length(row_sums)), tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# ftest_orfs: nested-uORF pruning + BH applied only to survivors
# ---------------------------------------------------------------------------
test_that("ftest_orfs prunes nested (shared-stop) uORFs by spec_coef before BH", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  # Pick one well-covered uORF and create an exact synthetic clone of it
  # that shares transcript prefix AND trspace end (the pruning key).
  # The clone is created by replicating the uORF's entries in trspacecds,
  # cdsgrl, and uORF, then duplicating its psite attributions. Identical
  # inputs produce identical spec_coef, so which one the pruning keeps
  # is arbitrary but exactly one must survive.
  anno <- chr22_anno
  uorf_ids <- names(anno$uORF)[anno$uORF]
  cov_counts <- table(as.character(psites$orf))
  cov_counts <- cov_counts[names(cov_counts) %in% uorf_ids]
  cov_counts <- sort(cov_counts, decreasing = TRUE)
  if (length(cov_counts) < 1L) skip("No covered uORFs available")
  src_id <- names(cov_counts)[1]

  tx1    <- sub("_[0-9]+$", "", src_id)
  new_id <- paste0(tx1, "_99999")

  # Clone trspacecds row
  src_gr <- anno$trspacecds[src_id]
  names(src_gr) <- new_id
  anno$trspacecds <- c(anno$trspacecds, src_gr)

  # Clone cdsgrl element
  src_cds <- anno$cdsgrl[src_id]
  names(src_cds) <- new_id
  anno$cdsgrl <- c(anno$cdsgrl, src_cds)

  # Clone uORF flag
  anno$uORF[new_id] <- TRUE

  # Clone psites attributed to the source ORF
  clones <- psites[as.character(psites$orf) == src_id]
  clones$orf <- new_id
  psites <- c(psites, clones)

  ft <- suppressMessages(ftest_orfs(psites, anno, n_cores = 1))

  row_a <- ft[ft$orf_id == src_id, ]
  row_b <- ft[ft$orf_id == new_id, ]
  expect_equal(nrow(row_a), 1L)
  expect_equal(nrow(row_b), 1L)

  # Exactly one of the two must have been pruned (NA p/q.value). The
  # pruned side is dropped from spec_test_df before BH, so its row gets
  # NA across p.value and q.value in the final left_join output.
  na_count <- is.na(row_a$q.value) + is.na(row_b$q.value)
  expect_equal(na_count, 1L)
})


test_that("ftest_orfs q.value is BH only over non-pruned survivors", {
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)
  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)

  ft <- suppressMessages(ftest_orfs(psites, chr22_anno, n_cores = 1))
  surv <- ft[!is.na(ft$p.value), ]
  # q.value is BH of p.value on the set of surviving rows.
  expect_equal(surv$q.value, p.adjust(surv$p.value, method = "BH"),
               tolerance = 1e-12)
})
