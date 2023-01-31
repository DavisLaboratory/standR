test_that("drawPCA works", {
  data("dkd_spe_subset")

  # Default param
  p <- drawPCA(dkd_spe_subset)
  expect_silent(print(p))

  # color with column
  p <- drawPCA(dkd_spe_subset, color = region)
  expect_silent(print(p))

  # color with expression
  p <- drawPCA(dkd_spe_subset, color = SequencingSaturation > 90)
  expect_silent(print(p))

  # multiple aesthetics
  p <- drawPCA(dkd_spe_subset, color = SequencingSaturation > 90, shape = region)
  expect_error(print(p), NA)

  # error when the column is not in the data
  p <- drawPCA(dkd_spe_subset, color = xyz)
  expect_error(
    # We need to force eval here to throw the error
    print(p),
    "object 'xyz' not found"
  )

  # will not error if the variable is present in the parent env
  region2 <- dkd_spe_subset$region
  p <- drawPCA(dkd_spe_subset, color = region2)
  expect_error(print(p), NA)

  # test precompute
  library(scater)
  spe <- suppressWarnings(scater::runPCA(dkd_spe_subset))
  pdata <- reducedDim(spe, "PCA")

  p <- drawPCA(dkd_spe_subset, precomputed = pdata)
  expect_silent(print(p))

  expect_error(drawPCA(dkd_spe_subset, precomputed = "pdata"))

  # test dgelist

  dge <- edgeR::SE2DGEList(spe)

  p <- drawPCA(dge)
  expect_silent(print(p))
})

test_that("plotMDS works", {
  data("dkd_spe_subset")

  # Default param
  p <- standR::plotMDS(dkd_spe_subset)
  expect_silent(print(p))

  # color with column
  p <- standR::plotMDS(dkd_spe_subset, color = region)
  expect_silent(print(p))

  # color with expression
  p <- standR::plotMDS(dkd_spe_subset, color = SequencingSaturation > 90)
  expect_silent(print(p))

  # multiple aesthetics
  p <- standR::plotMDS(dkd_spe_subset, color = SequencingSaturation > 90, shape = region)
  expect_error(print(p), NA)

  # error when the column is not in the data
  p <- standR::plotMDS(dkd_spe_subset, color = xyz)
  expect_error(
    # We need to force eval here to throw the error
    print(p),
    "object 'xyz' not found"
  )

  # will not error if the variable is present in the parent env
  region2 <- dkd_spe_subset$region
  p <- standR::plotMDS(dkd_spe_subset, color = region2)
  expect_error(print(p), NA)

  # test dgelist

  dge <- edgeR::SE2DGEList(dkd_spe_subset)

  p <- standR::plotMDS(dge)
  expect_silent(print(p))
})


test_that("plotDR works", {
  library(scater)
  data("dkd_spe_subset")

  dkd_spe_subset <- suppressWarnings(scater::runPCA(dkd_spe_subset))

  # Default param
  p <- plotDR(dkd_spe_subset, dimred = "PCA")
  expect_silent(print(p))

  # color with column
  p <- plotDR(dkd_spe_subset, color = region, dimred = "PCA")
  expect_silent(print(p))

  # color with expression
  p <- plotDR(dkd_spe_subset, color = SequencingSaturation > 90, dimred = "PCA")
  expect_silent(print(p))

  # multiple aesthetics
  p <- plotDR(dkd_spe_subset, color = SequencingSaturation > 90, shape = region, dimred = "PCA")
  expect_error(print(p), NA)

  # error when the column is not in the data
  p <- plotDR(dkd_spe_subset, color = xyz, dimred = "PCA")
  expect_error(
    # We need to force eval here to throw the error
    print(p),
    "object 'xyz' not found"
  )

  # will not error if the variable is present in the parent env
  region2 <- dkd_spe_subset$region
  p <- plotDR(dkd_spe_subset, color = region2, dimred = "PCA")
  expect_error(print(p), NA)
})
