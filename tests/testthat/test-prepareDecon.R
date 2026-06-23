test_that("Testing the function of preparation of spatialDecon", {
  library(ExperimentHub)
  eh <- ExperimentHub()
   
  query(eh, "standR")
  countFile <- eh[["EH7364"]]
  sampleAnnoFile <- eh[["EH7365"]]
   
  spe <- readGeoMx(countFile, sampleAnnoFile, rmNegProbe = FALSE)
 
  out <- prepareSpatialDecon(spe)
  
  expect_equal(length(out), 2)

})

test_that("prepareSpatialDecon uses negative probes stored in metadata", {
  countFile <- data.frame(
    TargetName = c("GeneA", "GeneB", "NegProbe-WTX"),
    ROI1 = c(10, 20, 1),
    ROI2 = c(30, 40, 2),
    check.names = FALSE
  )
  sampleAnnoFile <- data.frame(
    SegmentDisplayName = c("ROI1", "ROI2"),
    ROICoordinateX = c(1, 2),
    ROICoordinateY = c(3, 4)
  )
  spe <- readGeoMx(countFile, sampleAnnoFile, rmNegProbe = TRUE)

  out <- prepareSpatialDecon(spe)

  expect_equal(names(out), c("normCount", "backGround"))
  expect_equal(dim(out$normCount), c(2, 2))
  expect_equal(dim(out$backGround), c(2, 2))
})

