test_that("prepareSpatialDecon uses negative probes retained in the assay", {
  countFile <- test_path("fixtures", "geomx-counts.tsv")
  sampleAnnoFile <- test_path("fixtures", "geomx-samples.tsv")
  spe <- readGeoMx(countFile, sampleAnnoFile, rmNegProbe = FALSE)

  out <- prepareSpatialDecon(spe)

  expected_counts <- matrix(
    c(10, 20, 1, 30, 40, 2),
    nrow = 3,
    dimnames = list(c("GeneA", "GeneB", "NegProbe-WTX"), c("ROI1", "ROI2"))
  )
  expected_norm <- edgeR::cpm(expected_counts, log = TRUE)
  expected_background <- matrix(
    rep(expected_norm["NegProbe-WTX", ], each = 3),
    nrow = 3,
    dimnames = dimnames(expected_norm)
  )

  expect_equal(names(out), c("normCount", "backGround"))
  expect_equal(out$normCount, expected_norm)
  expect_equal(out$backGround, expected_background)
})

test_that("prepareSpatialDecon uses negative probes stored in metadata", {
  countFile <- test_path("fixtures", "geomx-counts.tsv")
  sampleAnnoFile <- test_path("fixtures", "geomx-samples.tsv")
  spe <- readGeoMx(countFile, sampleAnnoFile, rmNegProbe = TRUE)

  out <- prepareSpatialDecon(spe)

  expected_counts <- matrix(
    c(10, 20, 1, 30, 40, 2),
    nrow = 3,
    dimnames = list(c("GeneA", "GeneB", "NegProbe-WTX"), c("ROI1", "ROI2"))
  )
  expected_norm <- edgeR::cpm(expected_counts[1:2, , drop = FALSE], log = TRUE)
  expected_neg <- edgeR::cpm(expected_counts, log = TRUE)["NegProbe-WTX", ]
  expected_background <- matrix(
    rep(expected_neg, each = 2),
    nrow = 2,
    dimnames = dimnames(expected_norm)
  )

  expect_equal(names(out), c("normCount", "backGround"))
  expect_equal(out$normCount, expected_norm)
  expect_equal(out$backGround, expected_background)
})

