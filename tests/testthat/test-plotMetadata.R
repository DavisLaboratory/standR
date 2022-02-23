test_that("Testing plotMetadata function", {
  library(ggalluvial)
  data("dkd_spe_subset")

  expect_error(plotMetadata(dkd_spe_subset, column2plot = c("SlideName","disease_status","xyz")))
})
