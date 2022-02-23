test_that("Testing import data from path", {
  library(SpatialExperiment)
  dirPath <- ""
  countFile <- system.file("extdata", "dkd_subset_TargetCountMatrix.txt", package = "standR")
  sampleAnnoFile <- system.file("extdata", "dkd_subset_Sample_Annotations.txt", package = "standR")

  # default
  spe <- geomx_import_from_path(dirPath, countFile, sampleAnnoFile, hasNegProbe = FALSE)
  expect_equal(nrow(spe),3000)
  expect_equal(ncol(spe),70)
  expect_equal(length(assayNames(spe)),2)

  # change colnames to use
  expect_error(geomx_import_from_path(dirPath, countFile, sampleAnnoFile,
                                      hasNegProbe = FALSE, colnames.as.rownames = c("x","y")))
  expect_error(geomx_import_from_path(dirPath, countFile, sampleAnnoFile,
                                      hasNegProbe = FALSE, colnames.as.rownames = c("x","y","z")))
  expect_error(geomx_import_from_path(dirPath, countFile, sampleAnnoFile,
                                      hasNegProbe = FALSE, colnames.as.rownames = c("TargetName", "SegmentDisplayName")))
  expect_identical(geomx_import_from_path(dirPath, countFile, sampleAnnoFile,
                                        hasNegProbe = FALSE, colnames.as.rownames = c("TargetName", "SegmentDisplayName", "z")),spe)
})
