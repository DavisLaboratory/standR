test_that("Testing plotlibraryQC function", {
  library(ggplot2)
  library(patchwork)
  data("dkd_spe_subset")
  spe <- addQCstat(dkd_spe_subset)

  expect_message(plotLibrarySizeQC(spe))
  expect_error(plotLibrarySizeQC(spe, col = xyz))
})
