test_that("Testing plotGeneQC", {
  library(ggplot2)
  data("dkd_spe_subset")
  spe <- addPerROIQC(dkd_spe_subset)
  
  p <- plotGeneQC(spe)
  expect_silent(print(p))
  expect_error(plotGeneQC(col = xyz))
})
