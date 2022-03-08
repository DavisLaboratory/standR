test_that("Testing plotGeneQC", {
  library(ggplot2)
  data("dkd_spe_subset")
  spe <- addPerROIQC(dkd_spe_subset)

  expect_silent(plotGeneQC(spe))
  expect_error(plotGeneQC(col = xyz))
})
