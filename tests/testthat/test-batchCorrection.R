test_that("Testing negative control gene function", {
  data("dkd_spe_subset")

  spe <- chooseNegCTRLgenes(dkd_spe_subset, top_n = 100)

  expect_identical(dim(spe), dim(dkd_spe_subset))
  expect_true("mean_zscore" %in% colnames(rowData(spe)))
  expect_true("mean_expr" %in% colnames(rowData(spe)))
  expect_true("negGenes" %in% names(metadata(spe)))
  expect_equal(length(metadata(spe)$negGenes), 100)

  expect_error(chooseNegCTRLgenes(dkd_spe_subset, n_assay = 3))
  expect_error(chooseNegCTRLgenes(dkd_spe_subset, batch_name = "xyz"))
  expect_error(chooseNegCTRLgenes(dkd_spe_subset, top_n = 3001))

})


test_that("Testing ruv4 batch correction function", {
  data("dkd_spe_subset")

  spe <- chooseNegCTRLgenes(dkd_spe_subset, top_n = 100)
  spe_ruv <- runRUV4(spe, k = 3, factors = c("disease_status","region"),
                     negctrlGenes = metadata(spe)$negGenes)

  expect_identical(dim(spe_ruv), dim(dkd_spe_subset))
  expect_gt(ncol(colData(spe_ruv)), ncol(colData(spe)))
  expect_false(assay(spe,2)[1,1] == assay(spe_ruv,2)[1,1])
  expect_equal(ncol(colData(spe_ruv))-ncol(colData(spe)),3)

  expect_error(runRUV4(spe, k = -1, factors = c("disease_status","region"),
                       negctrlGenes = metadata(spe)$negGenes))
  expect_error(runRUV4(spe, k = 3, factors = c("xyz"),
                       negctrlGenes = metadata(spe)$negGenes))
})


test_that("Testing combat batch correction function", {
  data("dkd_spe_subset")

  spe_combat <- runCombat(dkd_spe_subset, batch = colData(dkd_spe_subset)$SlideName, bio_factor = colData(dkd_spe_subset)$region)

  expect_identical(dim(spe_combat), dim(dkd_spe_subset))
  expect_equal(ncol(colData(spe_combat)), ncol(colData(dkd_spe_subset)))

  expect_error(runCombat(dkd_spe_subset, batch = colData(dkd_spe_subset)$SlideNames, bio_factor = colData(dkd_spe_subset)$region))
})


test_that("Testing limma batch correction function", {
  data("dkd_spe_subset")

  spe_limmarmb <- runLimmarmb(dkd_spe_subset, batch = colData(dkd_spe_subset)$SlideName)

  expect_identical(dim(spe_limmarmb), dim(dkd_spe_subset))
  expect_equal(ncol(colData(spe_limmarmb)), ncol(colData(dkd_spe_subset)))

  expect_error(runCombat(dkd_spe_subset, batch = colData(dkd_spe_subset)$SlideNames))
})


