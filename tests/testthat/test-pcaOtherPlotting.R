test_that("Testing PCA pair plot function", {
  data("dkd_spe_subset")

  p <- plotPairdimensionPCA(dkd_spe_subset)
  expect_silent(print(p))

  p <- plotPairdimensionPCA(dkd_spe_subset, color=region)
  expect_silent(print(p))

  p <- plotPairdimensionPCA(dkd_spe_subset, color=SequencingSaturation > 90)
  expect_silent(print(p))

  #multiple aesthetics
  p <- plotPairdimensionPCA(dkd_spe_subset, color=SequencingSaturation > 90,
                            shape=region)
  expect_error(print(p), NA)

  #error when the column is not in the data
  expect_error(plotPairdimensionPCA(dkd_spe_subset, color=xyz))

})

test_that("Testing PCA scree plot function", {
  data("dkd_spe_subset")

  p <- plotScreePCA(dkd_spe_subset)
  expect_silent(print(p))

})

test_that("Testing PCA biplot function", {
  data("dkd_spe_subset")

  p <- plotPCAbiplot(dkd_spe_subset)
  expect_silent(print(p))

  p <- plotPCAbiplot(dkd_spe_subset, n_loadings = c("RPL7A","ATP5MC2"))
  expect_silent(print(p))

})
