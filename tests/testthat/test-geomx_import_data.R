test_that("Testing import data from path", {
  library(SpatialExperiment)
  dirPath <- ""
  countFile <- system.file("extdata", "Kidney_Raw_TargetCountMatrix.txt.gz", package = "standR")
  sampleAnnoFile <- system.file("extdata", "Kidney_Sample_Annotations.txt.gz", package = "standR")

  # default
  spe <- geomx_import_from_path(dirPath, countFile, sampleAnnoFile, hasNegProbe = TRUE)
  expect_equal(nrow(spe),18503)
  expect_equal(ncol(spe),231)
  expect_equal(length(assayNames(spe)),2)

  # change colnames to use
  expect_error(geomx_import_from_path(dirPath, countFile, sampleAnnoFile,
                                      hasNegProbe = TRUE, colnames.as.rownames = c("x","y")))
  expect_error(geomx_import_from_path(dirPath, countFile, sampleAnnoFile,
                                      hasNegProbe = TRUE, colnames.as.rownames = c("x","y","z")))
  expect_error(geomx_import_from_path(dirPath, countFile, sampleAnnoFile,
                                      hasNegProbe = TRUE, colnames.as.rownames = c("TargetName", "SegmentDisplayName")))
  expect_identical(geomx_import_from_path(dirPath, countFile, sampleAnnoFile,
                                        hasNegProbe = TRUE, colnames.as.rownames = c("TargetName", "SegmentDisplayName", "z")),spe)
})
