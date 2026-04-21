# Ribostan — Developer Notes

> Working notes for ongoing development. For the public-facing vignette, see `vignettes/Ribostan.rmd`.

## What This Package Does

**Ribostan** is an R/Bioconductor package for analysing ribosomal profiling (ribo-seq) data. Starting from aligned RPF (ribosome-protected fragment) BAM files and a genome annotation, it provides:

- **P/A-site alignment** — offsets are chosen by maximising CDS inclusion (Ahmed et al. 2019) or by KL-divergence peaks in metacodon profiles (O'Connor et al. 2016).
- **Isoform-aware quantification** — Stan-based optimisation distributes multimapping reads across isoforms, producing RITPM values (ribosome-footprint TPMs).
- **uORF identification** — ORFik-based upstream ORF detection followed by multitaper 3-nt periodicity tests to filter for bona fide translation.
- **Codon occupancy / elongation** — RUST or GLM-based models for per-codon dwell times, aggregated to ORF- and gene-level elongation estimates.
- **Metacodon / KL plots** — diagnostic plots that verify offset choice independent of CDS-inclusion method.

### Pipeline in brief

```
BAM  ──► get_readgr()  ──► get_offsets()  ──► get_psite_gr()
                                                      │
                               ┌──────────────────────┤
                               ▼                      ▼
                       get_ritpms()          get_metacodon_profs()
                       gene_level_expr()     get_kl_df() / plot_kl_dv()
                               │
                               ▼
                       periodicity_filter_uORFs()
                               │
                               ▼
                       get_codon_occs() ──► get_orf_elong() ──► gene_level_elong()
```

---

## Package State (as of Apr 2026)

### Data / extdata

The source tree does **not** include binary test data (`inst/extdata/`, `data/*.rda`). These are expected to be present in an installed copy built from a release tarball. The data-raw scripts (`data-raw/create_testdata.R`, `data-raw/testing.R`) describe how to recreate them from chr22 of hg38.

Datasets shipped with the package (documented in `R/data.R`):

| Object | Description |
|--------|-------------|
| `chr22_anno` | Annotation list from `load_annotation()` on chr22 GENCODE v32, with uORFs |
| `rpfs` | Example `GRanges` of RPF reads |
| `offsets_df` | CDS-occupancy offsets for `rpfs` |
| `metacodondf` | Pre-computed metacodon tibble (~23 k rows) |
| `ms_df` | Mass-spec protein abundance for chr22 genes (used in `test-quant.R`) |

### Dependencies

Install status (R 4.5, Bioc 3.21):

| Package | Status |
|---------|--------|
| `rtracklayer`, `Rsamtools`, `GenomicFeatures`, `GenomicAlignments` | Install on first run via `BiocManager` |
| `ORFik` | Bioconductor — required for uORF detection |
| `rstan`, `multitaper` | CRAN |
| All others | Generally pre-installed in a Bioc environment |

```r
BiocManager::install(c('rtracklayer','Rsamtools','GenomicFeatures',
                       'GenomicAlignments','ORFik','multitaper','rstan'))
devtools::install_local(".")   # or remotes::install_github("ohlerlab/RiboStan")
```

---

## Known Issues / TODOs

### ~~CRITICAL: Three distinct bugs make periodicity testing unreliable~~ (FIXED)

All three bugs identified and fixed in this session. The pipeline is
now properly calibrated: null false-positive rate drops from **77% to
3.9%** at p<0.05, and the real data still correctly identifies
periodic ORFs. Regression tests live in
`tests/testthat/test-regression.R` and assert:

- `get_psite_gr` applies phaseshift in correct per-row order
  (phase-0 reads must land at frame 0 of their assigned ORF ≥97% of
  the time)
- `ftest_orfs` returns actual F-test p-values (bounded in [0,1] and
  monotone in Fstat)
- The full pipeline gives roughly uniform p-values under a
  randomized-read null (median p > 0.3, fraction p<0.05 < 10%)

Details preserved below for reference.



Empirical null test: randomize 5' positions of all RPFs uniformly
within host transcripts (no true periodicity), compute ACTUAL F-test
p-values:

| Method | Real p<0.05 | Null p<0.05 | Null median p |
|---|---|---|---|
| Current pipeline (merge bug) | 72% | **77%** | 0.0009 |
| Merge bug fixed, shifts recomputed per call | 21% | 14% | 0.34 |
| **Merge bug fixed + shifts calibrated once on real data** | **20.5%** | **4.5%** | **0.46** |
| Phase-independent (no `shift(-phase)`) | 21% | 4.1% | 0.47 |

The existing `shift(-phase) + shift(shft)` design is mathematically
sound — it's a clean permutation of `(readlen, phase) → frame` that
preserves real signal and gives uniform null when `shft` is
correctly applied. The three root causes of the broken behaviour
are all mechanical:

**Bug 1 — `merge()` scrambles phaseshift row alignment
(`get_ritpms.R:287`) — PRIMARY cause of the tautology.**
```r
phaseshifts <- merge(mcols(psites)[, c("phase","readlen")], shiftdf, all.x = TRUE)
psites <- psites %>% shift(-.$phase) %>% shift(phaseshifts$shft)
```
`merge()` sorts by the join keys, so `phaseshifts$shft` is no longer
aligned with `psites` — every P-site gets a shift value meant for
some other P-site. Under this bug the final frame of a P-site is
essentially a random draw from the global `shft` distribution,
independent of its own phase. Fixing only this bug (replacing
`merge()` with an order-preserving `left_join`) drops the null
false-positive rate from **77% → 14% at p<0.05**, which demonstrates
the `shift(-phase) + shift(shft)` design was conceptually sound but
destroyed by the wrong join. A proper left-join preserving input
order is the minimum required fix.

**Bug 2 — `shiftdf` is recomputed on every call to `get_psite_gr`
(`get_ritpms.R:278-285`).**
The `shiftdf` (per-`(readlen, phase)` rank-based shift) is rebuilt
from whatever `psites` are passed in. This is wrong: the shift
parameters should be CALIBRATED ONCE on data with strong periodicity
(the full real dataset, like offsets), then applied unchanged to any
test set (novel ORFs, held-out transcripts, null data, etc.).
Recomputing `shiftdf` on small or non-periodic subsets gives
essentially arbitrary rank assignments because (phase, readlen)
counts are near-uniform, introducing frame biases the F-test picks
up. Empirical test confirms: using `shiftdf` fixed from real data
and applied (with Bug 1 also fixed) to randomized null reads gives
p<0.05 at 4.5% (perfectly calibrated) with a flat 33.6/33.3/33.1
frame distribution. The fix is to treat `shiftdf` as a model
parameter alongside `offsets_df`, computed once and passed in.

**Bug 3 — `ftest_orfs` returns spectral power labeled as p-value
(`periodicity_tests.R:50`).**
```r
return(c(Fmax_3nt, spect_3nt))      # F statistic + spectral power

pval <- pf(q = vals[1], df1 = 2, df2 = (2*24) - 2, lower.tail = FALSE)  # DEAD CODE
return(c(vals[2], pval))
```
The early `return()` makes the pval computation unreachable. In
`ftest_orfs`, the output is labeled `c("spec_coef", "p.value", "orf_id")`
but `p.value` is actually the raw spectral power at 1/3 Hz. This means
`periodicity_filter_uORFs(..., remove = filter(p.value < 0.05))` is
filtering on spectral power, not statistical significance.

**Status (Apr 2026): all three bugs fixed.** The installed
`get_psite_gr` now uses an order-preserving `left_join`, calibrates
`shiftdf` on annotated-CDS psites only, and `ftestvect` now returns
the actual F-test p-value alongside the F statistic and spectral
power. Empirically the null is calibrated at 3.9% p<0.05 with a
33/33/33 frame distribution; real-data periodicity detection is
preserved (19% of tested ORFs at p<0.05, 48/34/18 frame mix). Three
regression tests in `tests/testthat/test-regression.R` guard against
reintroducing any of these bugs.

**Short-ORF & sparse-coverage robustness.** Simulations confirm the
F-test is well-calibrated across ORF lengths (15-600 nt) and under
bursty/clustered nulls (5% false positives for smooth Poisson, 0-3%
for concentrated peaks). On real chr22 data, uORFs across all length
and P-site-count bins have null-like p-value distributions (median
0.4-0.65) while annotated CDS show strong signal especially with
high P-site counts. The zero-padding to n=50 for short vectors does
not create spurious 1/3 Hz leakage.

**Open suggestion.** `periodicity_filter_uORFs` uses a raw p<0.05
threshold without multiple-testing correction. With 5000+ candidate
uORFs this produces ~250 expected false positives even under a
perfectly calibrated test. Adding Benjamini-Hochberg correction
would give tighter control at the expected FDR.

See `scripts/debug_real_pvalues.R`, `scripts/debug_fixed_shifts.R`,
`scripts/debug_short_orf_calibration.R`, `scripts/debug_bursty_null.R`
for supporting analyses.



`get_psite_gr()` performs three steps that together make `ftest_orfs()`
structurally incapable of distinguishing a translated ORF from a random
one:

```r
# get_ritpms.R:220
rpfs$phase <- (start(rpfs) - start(orfs)[ov]) %% 3   # phase relative to assigned ORF

# get_ritpms.R:272-274
psites <- rpfs %>% resize(1,"start") %>% shift(.$p_offset)
# offsets are all %% 3 == 0, so P-site frame == read phase (relative to ORF)

# get_ritpms.R:293-295
psites <- psites %>% shift(-.$phase) %>% shift(phaseshifts$shft)
# step 3a (shift(-phase)) forces every P-site to frame 0 of its assigned ORF
# step 3b (shift(shft)) scatters P-sites to frames {0,1,2} by GLOBAL rank of
#   (readlen, phase) -- i.e. by CDS-dominated phase mix, not per-ORF evidence.
```

**Empirical proof (scripts/debug_null_periodicity.R):** randomize the
5' positions of all RPFs uniformly within their host transcripts (pure
noise), then run the full pipeline:

| metric | real reads | randomized reads (null) |
|---|---|---|
| median `ftest_orfs` p-value | 0.25 | **0.045** |
| frac p < 0.05 | 26% | **53%** |
| frame 0 fraction of psites | 48% | 43% |

Under a calibrated test the null should give uniform p-values (~5% at
p<0.05). Instead the pipeline calls **53%** of ORFs periodic on pure
noise, a ~10x false positive rate. The random-reads data even gets
*more* "periodic" calls than real data (52.9% vs 26.3%).

Additionally, there is a separate bug in the phaseshift step: `merge()`
at line 287 reorders rows by the (phase, readlen) keys, so the
`phaseshifts$shft` vector on line 295 is no longer aligned with
`psites` — each P-site gets an effectively random shift value from
another P-site.

**Fix (proposed):** `ftest_orfs()` should operate on psites computed
with a phase-independent offset (one offset per readlength, no
`shift(-phase)` step). This preserves the raw 3-nt signal of read 5'
ends relative to the candidate ORF, which is exactly what a
periodicity test needs.

### Bug: `get_psite_gr()` silently drops reads overlapping multiple ORFs (FIXED)

In `get_ritpms.R:268-269`, reads overlapping >1 ORF had their
randomly-chosen candidate's `p_offset` set to `NULL` before being put
back into `rpfs`. The subsequent `subset(!is.na(p_offset))` then
dropped them entirely. Fixed by removing the `uniq_mov_rpfs$p_offset
<- NULL` line so the general downstream `shift(p_offset)` applies
correctly. Regression test in `tests/testthat/test-regression.R`.

### Bug: `ftest_orfs()` discards psite filter (`periodicity_tests.R:81`)

```r
# CURRENT (wrong) — line 80's filtered `orfs` is immediately overwritten:
orfs <- intersect(unique(psites$orf), names(anno$trspacecds))
orfs <- anno$trspacecds[]          # <-- drops the filter, tests ALL CDSs

# SHOULD BE:
orfs <- anno$trspacecds[orfs]
```

This means `ftest_orfs()` always tests every ORF in the annotation, ignoring which ORFs actually have P-sites. The function still works (coverage is just zero for unpopulated ORFs) but is wasteful and produces misleading ORF-length normalisation.

### Typo in `get_metacodon_profs()` parallelisation check

```r
# Current — 'paralell' will never match an installed package name:
if ('paralell' %in% installed.packages()) ...
# Fix:
if ('parallel' %in% installed.packages()) ...
```

Consequence: metacodon profiling always runs single-threaded even when `parallel` is available.

### `testthat` in `Imports`

`R/anno.R` has `@import testthat` and `DESCRIPTION` lists `testthat` in `Imports`. It should live in `Suggests` only. This causes testthat to be a hard runtime dependency.

### Test files that download from AnnotationHub

`test-annotation.R`, `test-kl_div.R`, and `test-get_codon_occs.R` always re-download and re-export the chr22 GTF even when the local file already exists (the `if(!file.exists(...))` guards are commented out). On CI without a cache this is very slow; on an offline machine they simply fail.

---

## Test Suite Status

| Test file | Data source | Status (Apr 2026) |
|-----------|-------------|-----------------|
| `test-bamread.R` | `data(chr22_anno)` + `extdata/nchr22.bam` | **PASS** (4 assertions) |
| `test-psites.R` | `data(chr22_anno/rpfs/offsets_df)` | **PASS** |
| `test-offsets_work.R` | same | **PASS** |
| `test-multitaper.R` | same | **PASS** |
| `test-quant.R` | same + `data(ms_df)` | **PASS / SKIP** — ms_df not shipped; mass-spec correlation block skipped, `ftest_orfs` columns check passes |
| `test-regression.R` | same | **PASS (13 assertions)** — guards against the multi-ORF silent-drop bug, the merge-reorder bug, the dead-`return()` p-value bug, and the null-calibration regression |
| `test-annotation.R` | AnnotationHub + BSgenome download | Slow; needs network; guards commented out (always re-downloads) |
| `test-kl_div.R` | AnnotationHub + BSgenome download | Slow; needs network |
| `test-get_codon_occs.R` | AnnotationHub + BSgenome download | Slow; needs network |

Run only the fast tests (no network required, ~50 s):

```r
library(testthat); library(Ribostan)
library(GenomicRanges); library(dplyr); library(stringr)
test_dir("tests/testthat",
  filter = "bamread|psites|offsets|multitaper|quant",
  package = "Ribostan")
```

Note: tests must be run with the dependent package namespaces explicitly loaded (`library(GenomicRanges)` etc.) because S4 generics like `seqnames()`/`mcols()` are imported but not re-exported. Using `devtools::test()` from within the package directory handles this automatically.

---

## Analysis: CDS Outside Annotated CDS

This is the main development goal: using Ribostan to detect and quantify translation of ORFs that are **not** in the standard gene annotation (no GTF CDS entry).

### What already works

1. **uORFs** — `load_annotation(add_uorfs=TRUE)` uses `ORFik::findUORFs()` to detect ORFs in 5′ UTRs. These are handled identically to annotated CDS throughout the pipeline.

2. **Custom ORF input** — `load_annotation()` ultimately builds a `cdsgrl` (`GRangesList`) plus matching `exonsgrl`. If you can supply these for novel ORFs, the downstream pipeline is largely agnostic.

3. **`get_ribofasta_anno()`** — parses a Gencode-style extended FASTA (as produced by `make_ext_fasta()`) to create a minimal annotation for a custom ORF set. Useful when reads were aligned to a custom transcript FASTA.

### What would be needed for other non-annotated ORFs

The categories and approaches:

| ORF category | Detection strategy | Integration point |
|---|---|---|
| **uORFs** (5′ UTR) | `ORFik::findUORFs` (already done) | `load_annotation(add_uorfs=TRUE)` |
| **dORFs** (3′ UTR / downstream) | `ORFik::findORFs` on 3′ UTRs | Add to `load_annotation()` analogously to uORFs |
| **Internal / frameshifted ORFs** | Scan CDS in alt frames | Same |
| **Intergenic ORFs / smORFs** | Whole-transcriptome ORF scan | Build custom `cdsgrl` + `exonsgrl`, pass to pipeline |
| **lncRNA-embedded ORFs** | `ORFik::findORFs` on lncRNA exons | Same |

### Minimal recipe for a custom ORF set

```r
library(Ribostan)

# 1. Build annotation normally
anno <- load_annotation(gtf, fafile, add_uorfs = FALSE)

# 2. Define novel ORFs as a GRangesList (one element per ORF,
#    with transcript_id / gene_id / gene_name mcols on unlisted ranges,
#    and each ORF's exon structure in exonsgrl):
novel_cdsgrl  <- ...   # GRangesList, named by ORF id
novel_exonsgrl <- ...  # GRangesList, named by parent transcript id

# 3. Merge into the annotation
anno$cdsgrl   <- c(anno$cdsgrl,  novel_cdsgrl)
anno$exonsgrl <- c(anno$exonsgrl, novel_exonsgrl)
anno$uORF     <- c(anno$uORF, setNames(rep(FALSE, length(novel_cdsgrl)),
                                       names(novel_cdsgrl)))
# Recompute trspacecds and derived slots
anno$trspacecds <- Ribostan:::get_trspace_cds(anno$cdsgrl, anno$exonsgrl[
  Ribostan:::fmcols(anno$cdsgrl, transcript_id)])
anno$cdsstarts  <- GenomicRanges::start(anno$trspacecds) |>
  setNames(names(anno$trspacecds))
anno$cds_prestop_st <- GenomicRanges::end(anno$trspacecds) - 2L |>
  setNames(names(anno$trspacecds))

# 4. Proceed normally
psites  <- get_psite_gr(rpfs, offsets_df, anno)
ritpms  <- get_ritpms(psites, anno)
ftests  <- ftest_orfs(psites, anno)
```

A cleaner path would be a new helper function `add_orfs(anno, new_cdsgrl, new_exonsgrl)` that wraps the above — worth adding if this use-case is confirmed.

### Periodicity testing as a discovery screen

`ftest_orfs()` returns a per-ORF spectral coefficient and p-value. After fixing the bug noted above, it can be run on a custom annotation containing candidate ORFs to rank them by translation evidence:

```r
candidate_anno <- load_annotation(gtf, fafile, add_uorfs = FALSE)
# ... add novel ORFs as above ...
psites  <- get_psite_gr(rpfs, offsets_df, candidate_anno)
ftests  <- ftest_orfs(psites, candidate_anno)
# Filter for periodic (likely translated) ORFs
translated <- ftests |> dplyr::filter(p.value < 0.05)
```

---

## File Map

```
R/
  anno.R              — GTF import, CDS filtering, uORF detection, load_annotation()
  get_offsets.R       — CDS-inclusion offset method
  get_ritpms.R        — BAM reading, P-site assignment, Stan quantification
  kl_div.R            — Metacodon profiles, KL divergence, offset verification
  model_TE.R          — Codon occupancy (RUST/GLM), elongation rates
  periodicity_tests.R — Multitaper 3-nt periodicity, uORF filtering
  data.R              — Dataset documentation
  utils-pipe.R        — %>% re-export

vignettes/Ribostan.rmd   — Full worked example (chr22, hg38)
tests/testthat/          — 8 test files (see table above)
data-raw/                — Scripts used to build shipped test data
```
