test_that("Testing negative control gene function", {
  data("dkd_spe_subset")

  spe <- findNCGs(dkd_spe_subset, top_n = 100)

  expect_identical(dim(spe), dim(dkd_spe_subset))
  expect_true("mean_zscore" %in% colnames(rowData(spe)))
  expect_true("mean_expr" %in% colnames(rowData(spe)))
  expect_true("NCGs" %in% names(metadata(spe)))
  expect_equal(length(metadata(spe)$NCGs), 100)

  expect_error(findNCGs(dkd_spe_subset, n_assay = 3))
  expect_error(findNCGs(dkd_spe_subset, batch_name = "xyz"))
  expect_error(findNCGs(dkd_spe_subset, top_n = 3001))

})


test_that("Testing ruv4 batch correction function", {
  data("dkd_spe_subset")

  spe <- findNCGs(dkd_spe_subset, top_n = 100)
  spe_ruv <- geomxBatchCorrection(spe, k = 3, factors = c("disease_status","region"),
                                  NCGs = S4Vectors::metadata(spe)$NCGs)

  expect_identical(dim(spe_ruv), dim(dkd_spe_subset))
  expect_gt(ncol(colData(spe_ruv)), ncol(colData(spe)))
  expect_false(assay(spe,2)[1,1] == assay(spe_ruv,2)[1,1])
  expect_equal(ncol(colData(spe_ruv))-ncol(colData(spe)),3)

  expect_error(geomxBatchCorrection(spe, k = -1, factors = c("disease_status","region"),
                                    NCGs = metadata(spe)$NCGs))
  expect_error(geomxBatchCorrection(spe, k = 3, factors = c("xyz"),
                                    NCGs = metadata(spe)$NCGs))
})


test_that("Testing limma batch correction function", {
  data("dkd_spe_subset")

  spe_limmarmb <- geomxBatchCorrection(dkd_spe_subset, batch = colData(dkd_spe_subset)$SlideName, method = "Limma")

  expect_identical(dim(spe_limmarmb), dim(dkd_spe_subset))
  expect_equal(ncol(colData(spe_limmarmb)), ncol(colData(dkd_spe_subset)))

})


