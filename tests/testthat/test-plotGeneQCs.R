test_that("Testing plotGeneQC", {
  library(ggplot2)
  data("dkd_spe_subset")
  spe <- addQCstat(dkd_spe_subset)

  expect_silent(plotGenewiseQC(spe))
  expect_error(plotGenewiseQC(col = xyz))
})
