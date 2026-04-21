test_that("test optimisation of quantification works", 
{
  data(chr22_anno)
  data(rpfs)
  data(offsets_df)

  psites <- get_psite_gr(rpfs, offsets_df, chr22_anno)
  # now, use Stan to estimate normalized p-site densities for our data
  ritpms <- get_ritpms(psites, chr22_anno)
  # and get these at the gene level (ignoring uORFs)
  gritpms <- gene_level_expr(ritpms, chr22_anno)

  # compare these to mass spec data if ms_df is available
  # (ms_df was generated from a private dataset and is not shipped with the package)
  ms_df_available <- tryCatch({ data(ms_df); TRUE }, warning = function(w) FALSE)
  if (ms_df_available && exists("ms_df")) {
    gritpms_vs <- gritpms %>%
      mutate(across("gene_id", ~stringr::str_replace(., "\\.\\d+$", "")))
    compdf <- ms_df %>%
      left_join(gritpms_vs) %>%
      filter(is.finite(log2(expr))) %>%
      filter(ribo > 0)
    expect_true(nrow(compdf) == 39)
    expect_true(cor(compdf$ribo, log2(compdf$expr)) > 0.66)
  } else {
    skip("ms_df not available — skipping mass-spec correlation check")
  }

  ftests <- ftest_orfs(psites %>% head(10000), chr22_anno, n_cores=1)
  expect_equal(colnames(ftests), c("orf_id", "spec_coef", "p.value"))
}
)
