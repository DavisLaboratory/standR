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


