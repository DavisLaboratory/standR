test_that("Testing plotlibraryQC function", {
  library(ggplot2)
  library(patchwork)
  data("dkd_spe_subset")
  spe <- addPerROIQC(dkd_spe_subset)
  
  p <- plotROIQC(spe)
  expect_s3_class(p, "ggplot")
  expect_error(plotROIQC(col = xyz))
})
