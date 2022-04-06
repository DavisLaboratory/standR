test_that("Testing plotlibraryQC function", {
  library(ggplot2)
  library(patchwork)
  data("dkd_spe_subset")
  spe <- addPerROIQC(dkd_spe_subset)

  expect_silent(plotROIQC(spe))
  expect_error(plotROIQC(col = xyz))
})
