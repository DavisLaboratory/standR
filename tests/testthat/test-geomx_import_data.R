test_that("Testing import data from path", {
  library(SpatialExperiment)
  url <- "http://nanostring-public-share.s3-website-us-west-2.amazonaws.com/"
  countFile <- paste0(url, "GeoScriptHub/KidneyDataset/Kidney_Raw_TargetCountMatrix.txt")
  sampleAnnoFile <- paste0(url, "GeoScriptHub/KidneyDataset/Kidney_Sample_Annotations.txt")

  # default
  spe <- readGeoMx(countFile, sampleAnnoFile, rmNegProbe = FALSE)
  expect_equal(nrow(spe), 18504)
  expect_equal(ncol(spe), 231)
  expect_equal(length(assayNames(spe)), 2)

  # change colnames to use
  expect_error(readGeoMx(countFile, sampleAnnoFile,
                         rmNegProbe = FALSE, colnames.as.rownames = c("x", "y")
  ))
  expect_error(readGeoMx(countFile, sampleAnnoFile,
                         rmNegProbe = FALSE, colnames.as.rownames = c("x", "y", "z")
  ))
  expect_error(readGeoMx(countFile, sampleAnnoFile,
                         rmNegProbe = FALSE, colnames.as.rownames = c("TargetName", "SegmentDisplayName")
  ))
  expect_identical(readGeoMx(countFile, sampleAnnoFile,
                             rmNegProbe = FALSE, colnames.as.rownames = c("TargetName", "SegmentDisplayName", "z")
  ), spe)
})
