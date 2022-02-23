test_that("Testing addQCstat", {
  data("dkd_spe_subset")

  spe <- addQCstat(dkd_spe_subset)

  expect_equal(ncol(colData(spe))-ncol(colData(dkd_spe_subset)),3)
  expect_equal("lib_size" %in% colnames(colData(spe)),TRUE)
  expect_equal("countOfLowEprGene" %in% colnames(colData(spe)),TRUE)
  expect_equal("percentOfLowEprGene" %in% colnames(colData(spe)),TRUE)

  expect_error(addQCstat(dkd_spe_subset, sample_fraction = 1.2))
  expect_error(addQCstat(dkd_spe_subset, sample_fraction = -1))
  expect_error(addQCstat(dkd_spe_subset, min_count = -1))

  expect_true(is.numeric(colData(spe)$lib_size))
  expect_true(is.numeric(colData(spe)$countOfLowEprGene))
  expect_true(is.numeric(colData(spe)$percentOfLowEprGene))
})
