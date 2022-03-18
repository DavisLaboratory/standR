test_that("Testing import data from path", {
  library(SpatialExperiment)
  dirPath <- system.file("extdata", package = "standR")
  countFile <- file.path(dirPath, "dkd_subset_TargetCountMatrix.txt")
  sampleAnnoFile <- file.path(dirPath, "dkd_subset_Sample_Annotations.txt")

  # default
  spe <- readGeoMx(countFile, sampleAnnoFile, hasNegProbe = FALSE)
  expect_equal(nrow(spe), 3000)
  expect_equal(ncol(spe), 70)
  expect_equal(length(assayNames(spe)), 2)

  # change colnames to use
  expect_error(readGeoMx(countFile, sampleAnnoFile,
    hasNegProbe = FALSE, colnames.as.rownames = c("x", "y")
  ))
  expect_error(readGeoMx(countFile, sampleAnnoFile,
    hasNegProbe = FALSE, colnames.as.rownames = c("x", "y", "z")
  ))
  expect_error(readGeoMx(countFile, sampleAnnoFile,
    hasNegProbe = FALSE, colnames.as.rownames = c("TargetName", "SegmentDisplayName")
  ))
  expect_identical(readGeoMx(countFile, sampleAnnoFile,
    hasNegProbe = FALSE, colnames.as.rownames = c("TargetName", "SegmentDisplayName", "z")
  ), spe)
})
