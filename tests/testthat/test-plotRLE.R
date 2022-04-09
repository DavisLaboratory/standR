test_that("plotRLExpr works with Eset", {
  data("dkd_spe_subset")

  # Default param
  p <- plotRLExpr(dkd_spe_subset)
  expect_silent(print(p))

  # color with column
  p <- plotRLExpr(dkd_spe_subset, color = SlideName)
  expect_silent(print(p))

  # color with expression
  p <- plotRLExpr(dkd_spe_subset, color = SequencingSaturation > 90)
  expect_silent(print(p))

  # multiple aesthetics
  p <- plotRLExpr(dkd_spe_subset, color = SequencingSaturation > 90, shape = disease_status)
  expect_error(print(p), NA)

  # order with one column
  p <- plotRLExpr(dkd_spe_subset, ordannots = "region")
  expect_silent(print(p))
  p <- plotRLExpr(dkd_spe_subset, ordannots = c("disease_status"), color = disease_status)
  expect_silent(print(p))

  # order with multiple columns
  p <- plotRLExpr(dkd_spe_subset, ordannots = c("region", "disease_status"), color = disease_status)
  expect_silent(print(p))

  # error when the column is not in the data
  p <- plotRLExpr(dkd_spe_subset, color = region2)
  expect_error(
    # We need to force eval here to throw the error
    print(p),
    "object 'region2' not found"
  )

  # will not error if the variable is present in the parent env
  region2 <- dkd_spe_subset$region
  p <- plotRLExpr(dkd_spe_subset, color = region2)
  expect_error(print(p), NA)

  # SCE
  sce <- SingleCellExperiment::SingleCellExperiment(list("a" = matrix(rnorm(6e3), 3)))
  expect_silent(print(plotRLExpr(sce[, 1:10])))
  expect_silent(print(plotRLExpr(sce[, 1:100])))
  expect_silent(print(plotRLExpr(sce)))


  # test dgelist

  dge <- edgeR::SE2DGEList(dkd_spe_subset)
  p <- plotRLExpr(dkd_spe_subset)
  expect_silent(print(p))
})
