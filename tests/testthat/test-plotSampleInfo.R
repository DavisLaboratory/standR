test_that("Testing plotSampleInfo function", {
  library(ggalluvial)
  data("dkd_spe_subset")

  expect_error(plotSampleInfo(dkd_spe_subset, column2plot = c("SlideName", "disease_status", "xyz")))
})
