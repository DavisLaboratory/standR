test_that("plotPCA works with Eset", {
  data("dkd_spe_subset")

  # Default param
  p <- plotPCA(dkd_spe_subset)
  expect_silent(print(p))

  # color with column
  p <- plotPCA(dkd_spe_subset, color = region)
  expect_silent(print(p))

  # color with expression
  p <- plotPCA(dkd_spe_subset, color = SequencingSaturation > 90)
  expect_silent(print(p))

  # multiple aesthetics
  p <- plotPCA(dkd_spe_subset, color = SequencingSaturation > 90, shape = region)
  expect_error(print(p), NA)

  # error when the column is not in the data
  p <- plotPCA(dkd_spe_subset, color = xyz)
  expect_error(
    # We need to force eval here to throw the error
    print(p),
    "object 'xyz' not found"
  )

  # will not error if the variable is present in the parent env
  region2 <- dkd_spe_subset$region
  p <- plotPCA(dkd_spe_subset, color = region2)
  expect_error(print(p), NA)
})
