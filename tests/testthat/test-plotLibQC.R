test_that("Testing plotlibraryQC function", {
  library(ggplot2)
  library(patchwork)
  data("dkd_spe_subset")
  spe <- addPerROIQC(dkd_spe_subset)

  expect_message(plotROIQC(spe))
  expect_error(plotROIQC(spe, col = xyz))
})
