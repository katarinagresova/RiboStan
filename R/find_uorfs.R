#' Find upstream ORFs (uORFs) from 5' UTR annotation
#'
#' Self-contained replacement for \code{ORFik::findUORFs}. Finds ORFs within
#' 5' UTR regions, optionally extending the search space into the downstream
#' CDS and filtering out ORFs that coincide with the annotated CDS.
#'
#' Procedure:
#' 1. Build search space: 5' UTRs, extended into the CDS if \code{cds} is
#'    supplied (so that ORFs bridging the 5' UTR / CDS boundary are found).
#' 2. Extract transcript sequences from the FASTA file.
#' 3. Scan each sequence for in-frame start → stop codon pairs.
#' 4. Map ORF positions from transcript coordinates back to genomic coordinates.
#' 5. Filter out ORFs that are fully within, share a start with, share a stop
#'    with, or have a start inside any annotated CDS element.
#'
#' @param fiveUTRs \code{GRangesList} of 5' UTR exons, one element per
#'   transcript.
#' @param fa Character path to an indexed genome FASTA file, or a
#'   \code{FaFile} object.
#' @param startCodon Character vector; start codons to search for. Defaults
#'   to \code{c("ATG","CTG","GTG")} so that non-canonical starts (CUG/GUG
#'   in RNA) are found in addition to ATG.
#' @param stopCodons Character vector; stop codons to recognise
#'   (default \code{c("TAA", "TAG", "TGA")}).
#' @param longestORF Logical; if \code{TRUE}, keep only the longest ORF sharing
#'   each stop codon position.
#' @param minimumLength Integer; minimum body length in nucleotides between the
#'   end of the start codon and the beginning of the stop codon.  \code{0}
#'   (the default) accepts any ORF consisting of just a start and a stop.
#' @param cds \code{GRangesList} of CDS exons.  When supplied, the CDS is
#'   appended to the 5' UTR search space and detected ORFs that are exact
#'   duplicates of an annotated CDS (sharing both start and stop position)
#'   are removed. Overlapping/nested ORFs are retained.
#' @return A \code{GRangesList} of uORFs in genomic coordinates, named
#'   \code{"<transcript_id>_<N>"}. Each element carries a metadata column
#'   \code{start_codon} giving the start codon that was matched.
#' @keywords internal
find_uorfs <- function(fiveUTRs, fa,
                       startCodon    = c("ATG", "CTG", "GTG"),
                       stopCodons    = c("TAA", "TAG", "TGA"),
                       longestORF    = FALSE,
                       minimumLength = 0L,
                       cds           = NULL) {
  .log_msg <- function(msg) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(str_interp("[${ts}] [find_uorfs] ${msg}"))
  }
  .log_msg(str_interp("starting uORF search: ${length(fiveUTRs)} transcripts, start=${paste(startCodon,',')}, stop=${paste(stopCodons,',')}, minLen=${minimumLength}, longest=${longestORF}"))

  ## 1. Build search space ------------------------------------------------
  .log_msg("building search space from 5' UTRs")
  search_space <- fiveUTRs
  if (!is.null(cds) && length(cds) > 0L) {
    .log_msg("extending search space with CDS regions")
    search_space <- .append_cds_to_utrs(fiveUTRs, cds)
  }
  # Ensure exon_rank is numeric (required by extractTranscriptSeqs)
  if ("exon_rank" %in% colnames(mcols(search_space@unlistData))) {
    mcols(search_space@unlistData)$exon_rank <- as.numeric(as.character(
      mcols(search_space@unlistData)$exon_rank
    ))
  }
  .log_msg(str_interp("search space: ${length(search_space)} transcripts"))

  ## 2. Extract transcript sequences ---------------------------------------
  .log_msg("extracting transcript sequences from FASTA")
  fa_file <- if (is.character(fa)) Rsamtools::FaFile(fa) else fa
  seqs <- toupper(as.character(
    GenomicFeatures::extractTranscriptSeqs(fa_file, search_space)
  ))

  ## 3. Find ORFs in transcript coordinates (via C++) ----------------------
  .log_msg("scanning sequences for ORFs (C++)")
  start_codons <- toupper(as.character(startCodon))
  raw <- find_orfs_cpp(
    seqs         = seqs,
    start_codons = start_codons,
    stop_codons  = toupper(stopCodons),
    min_body     = as.integer(minimumLength)
  )
  .log_msg(str_interp("found ${length(raw$starts)} raw ORFs"))

  if (length(raw$starts) == 0L) {
    .log_msg("no ORFs found, returning empty GRangesList")
    return(GenomicRanges::GRangesList())
  }

  # Optionally retain only the longest ORF per (sequence, stop-site) pair
  if (isTRUE(longestORF)) {
    .log_msg("filtering to keep only longest ORF per stop site")
    widths  <- raw$ends - raw$starts
    key     <- paste0(raw$indices, "_", raw$ends)
    keep    <- unlist(
      tapply(seq_along(key), key, function(i) i[which.max(widths[i])]),
      use.names = FALSE
    )
    raw <- lapply(raw, `[`, keep)
    .log_msg(str_interp("retained ${length(raw$starts)} ORFs after longestORF filtering"))
  }

  start_codon_hit <- start_codons[raw$start_codon_idx]

  ## 4. Map to genomic coordinates -----------------------------------------
  .log_msg("mapping ORFs from transcript to genomic coordinates")
  # tx_idx: which element of search_space each ORF belongs to (1-based)
  tx_idx     <- raw$indices
  all_orf_ir <- IRanges::IRanges(start = raw$starts, end = raw$ends)

  has_orfs   <- logical(length(search_space))
  has_orfs[unique(tx_idx)] <- TRUE
  tx_lengths <- tabulate(tx_idx, nbins = length(search_space))[unique(tx_idx)]

  # pmapFromTranscripts maps orf_gr[i] through search_space[tx_idx[i]]
  orf_gr      <- GenomicRanges::GRanges(
    seqnames = names(search_space)[tx_idx],
    ranges   = all_orf_ir
  )
  genomic_grl <- GenomicFeatures::pmapFromTranscripts(orf_gr,
                                                       search_space[tx_idx])
  .log_msg(str_interp("mapped ${length(genomic_grl)} ORFs to genome"))

  ## 5. Name ORFs: <txname>_<rank> -----------------------------------------
  .log_msg("assigning ORF names and start codon metadata")
  # Sequential rank within each transcript
  tx_names_rep <- names(search_space)[tx_idx]
  rank         <- unlist(
    tapply(seq_along(tx_idx), tx_idx, seq_along),
    use.names = FALSE
  )
  names(genomic_grl)  <- paste0(tx_names_rep, "_", rank)
  seqinfo(genomic_grl) <- seqinfo(fiveUTRs)
  # Record which start codon produced each uORF. Attach at both the GRL
  # level (one value per uORF) and the inner GRanges level (one value
  # per exon, same value for all exons of a uORF) so the column survives
  # the unlist() / split() cycle in load_annotation.
  S4Vectors::mcols(genomic_grl)$start_codon <- start_codon_hit
  S4Vectors::mcols(genomic_grl@unlistData)$start_codon <-
    rep(start_codon_hit, S4Vectors::elementNROWS(genomic_grl))

  ## 6. Filter ORFs that are exact duplicates of annotated CDS -------------
  if (!is.null(cds) && length(genomic_grl) > 0L) {
    .log_msg("filtering out ORFs that duplicate annotated CDS")
    n_before <- length(genomic_grl)
    genomic_grl <- .filter_uorfs(genomic_grl, cds)
    .log_msg(str_interp("removed ${n_before - length(genomic_grl)} duplicate ORFs"))
  }

  .log_msg(str_interp("uORF search complete: ${length(genomic_grl)} final uORFs"))
  genomic_grl
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Scan a single nucleotide string for ORFs
#'
#' @param seq Character; uppercase DNA sequence.
#' @param start_codon Character; start codon (default \code{"ATG"}).
#' @param stop_codons Character vector; stop codons.
#' @param min_length Integer; minimum body length in nucleotides.
#' @param longest_orf Logical; return only longest ORF per stop site.
#' @return \code{IRanges} of ORF coordinates (1-based; both start and stop
#'   codons included).
#' @keywords internal
.find_orfs_in_seq <- function(seq,
                               start_codon = "ATG",
                               stop_codons = c("TAA", "TAG", "TGA"),
                               min_length  = 0L,
                               longest_orf = FALSE) {
  n       <- nchar(seq)
  starts  <- gregexpr(start_codon, seq, fixed = TRUE)[[1L]]
  if (starts[[1L]] == -1L) return(IRanges::IRanges())

  orf_starts <- integer(0)
  orf_ends   <- integer(0)

  for (s in as.integer(starts)) {
    pos <- s + 3L
    while (pos + 2L <= n) {
      codon <- substr(seq, pos, pos + 2L)
      if (codon %in% stop_codons) {
        # body = nt between end of start and start of stop = pos - s - 3
        if ((pos - s - 3L) >= min_length) {
          orf_starts <- c(orf_starts, s)
          orf_ends   <- c(orf_ends, pos + 2L)
        }
        break
      }
      pos <- pos + 3L
    }
  }

  if (length(orf_starts) == 0L) return(IRanges::IRanges())

  if (longest_orf) {
    widths <- orf_ends - orf_starts
    keep   <- unlist(
      tapply(seq_along(orf_ends), orf_ends,
             function(i) i[which.max(widths[i])]),
      use.names = FALSE
    )
    orf_starts <- orf_starts[keep]
    orf_ends   <- orf_ends[keep]
  }

  IRanges::IRanges(start = orf_starts, end = orf_ends)
}


#' Append CDS exons to each 5' UTR to extend the ORF search space
#'
#' Allows detection of ORFs that start in the 5' UTR and extend into the CDS.
#' Only transcripts present in both \code{fiveUTRs} and \code{cds} are
#' extended; others are returned unchanged.
#'
#' @param fiveUTRs GRangesList of 5' UTR exons.
#' @param cds GRangesList of CDS exons.
#' @return GRangesList with CDS appended (and reduced) for matching transcripts.
#' @keywords internal
.append_cds_to_utrs <- function(fiveUTRs, cds) {
  shared <- intersect(names(fiveUTRs), names(cds))
  if (length(shared) == 0L) return(fiveUTRs)

  # Combine element-wise (equivalent to ORFik::pc)
  combined <- GenomicRanges::GRangesList(
    mapply(c, fiveUTRs[shared], cds[shared], SIMPLIFY = FALSE)
  )
  names(combined) <- shared

  reduced <- GenomicRanges::reduce(combined, ignore.strand = FALSE)
  names(reduced)   <- shared
  seqinfo(reduced) <- seqinfo(fiveUTRs)

  only_utrs <- fiveUTRs[setdiff(names(fiveUTRs), shared)]
  c(reduced, only_utrs)
}


#' Remove uORFs that exactly duplicate an annotated CDS
#'
#' A uORF is dropped only if it shares BOTH its start site and its stop site
#' with the same annotated CDS (i.e. it \emph{is} the CDS). Overlapping,
#' nested, and in-frame-extension uORFs are retained so that downstream
#' analyses can reason about them (see \code{load_annotation} which adds
#' per-uORF overlap metadata).
#'
#' @param uorfs GRangesList of candidate uORFs.
#' @param cds GRangesList of annotated CDS.
#' @return Filtered GRangesList.
#' @keywords internal
.filter_uorfs <- function(uorfs, cds) {
  if (length(uorfs) == 0L) return(uorfs)

  start_hits <- GenomicRanges::findOverlaps(
    .start_sites_gr(uorfs), .start_sites_gr(cds), type = "within"
  )
  stop_hits  <- GenomicRanges::findOverlaps(
    .stop_sites_gr(uorfs),  .stop_sites_gr(cds),  type = "within"
  )
  if (length(start_hits) == 0L || length(stop_hits) == 0L) return(uorfs)

  # uORF i duplicates CDS j iff (i, j) appears in BOTH hit sets.
  start_pairs <- paste(S4Vectors::queryHits(start_hits),
                       S4Vectors::subjectHits(start_hits), sep = "_")
  stop_pairs  <- paste(S4Vectors::queryHits(stop_hits),
                       S4Vectors::subjectHits(stop_hits),  sep = "_")
  dup_rows <- unique(
    S4Vectors::queryHits(start_hits)[start_pairs %in% stop_pairs]
  )
  if (length(dup_rows) > 0L) uorfs[-dup_rows] else uorfs
}


#' Width-1 GRanges at the genomic start (5'-most base) of each GRL element
#' @keywords internal
.start_sites_gr <- function(grl) {
  ul     <- unlist(grl, use.names = FALSE)
  is_neg <- as.character(GenomicRanges::strand(ul)) == "-"
  # + strand: start of the leftmost exon; - strand: end of the rightmost exon
  coords <- ifelse(is_neg,
                   GenomicRanges::end(ul),
                   GenomicRanges::start(ul))
  grp    <- factor(rep(seq_along(grl), S4Vectors::elementNROWS(grl)),
                   levels = seq_along(grl))
  sel    <- unlist(
    tapply(seq_along(ul), grp, function(g) {
      if (is_neg[g[[1L]]]) g[which.max(coords[g])]
      else                  g[which.min(coords[g])]
    }),
    use.names = FALSE
  )
  GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(ul)[sel],
    ranges   = IRanges::IRanges(start = coords[sel], width = 1L),
    strand   = GenomicRanges::strand(ul)[sel]
  )
}


#' Width-1 GRanges at the genomic stop (3'-most base) of each GRL element
#' @keywords internal
.stop_sites_gr <- function(grl) {
  ul     <- unlist(grl, use.names = FALSE)
  is_neg <- as.character(GenomicRanges::strand(ul)) == "-"
  # + strand: end of rightmost exon; - strand: start of leftmost exon
  coords <- ifelse(is_neg,
                   GenomicRanges::start(ul),
                   GenomicRanges::end(ul))
  grp    <- factor(rep(seq_along(grl), S4Vectors::elementNROWS(grl)),
                   levels = seq_along(grl))
  sel    <- unlist(
    tapply(seq_along(ul), grp, function(g) {
      if (is_neg[g[[1L]]]) g[which.min(coords[g])]
      else                  g[which.max(coords[g])]
    }),
    use.names = FALSE
  )
  GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(ul)[sel],
    ranges   = IRanges::IRanges(start = coords[sel], width = 1L),
    strand   = GenomicRanges::strand(ul)[sel]
  )
}
