test_that("Testing plotMetadata function", {
  data("dkd_spe_subset")

  expect_silent(plotMetadata(dkd_spe_subset, column2plot = c("SlideName","disease_status","region")))
})
