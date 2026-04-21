test_that("reading bam file works",
{
  suppressPackageStartupMessages(library(GenomicRanges))
  data(chr22_anno)
  testbam <- system.file("extdata", "nchr22.bam", package = "Ribostan", mustWork = TRUE)
  rpfs <- get_readgr(testbam, chr22_anno)
  expect_equal(length(rpfs), 63346)
  expect_true(all(GenomicRanges::seqnames(rpfs) %in% names(chr22_anno$exonsgrl)))
  expect_true(all(colnames(GenomicRanges::mcols(rpfs)) == "readlen"))
  expect_true(all(GenomicRanges::strand(rpfs) == "+"))
}
)
