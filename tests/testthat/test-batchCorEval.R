test_that("Testing the function of computing evaluation statistics for batch correction", {
  library(scater)
  data("dkd_spe_subset")
  spe <- scater::runPCA(dkd_spe_subset)
  eval_out <- computeClusterEvalStats(spe, "SlideName")

  expect_equal(nrow(eval_out), 3)
  expect_equal(ncol(eval_out), 2)

  expect_error(computeClusterEvalStats(spe, "xyz"))
  expect_error(computeClusterEvalStats(dkd_spe_subset, "SlideName"))
})

test_that("Testing the function of plotting and comparing different batch-corrected data", {
  library(scater)
  data("dkd_spe_subset")
  spe <- scater::runPCA(dkd_spe_subset)
  spe2 <- spe
  spe3 <- spe

  expect_silent(plotClusterEvalStats(list(spe, spe2, spe3), bio_feature_name = "region",
                       batch_feature_name = "SlideName", c("test1","test2","test3")))

  expect_error(plotClusterEvalStats(list(spe, spe2, spe3), bio_feature_name = "xyz",
                                     batch_feature_name = "SlideName", c("test1","test2","test3")))

  expect_error(plotClusterEvalStats(list(spe, spe2, spe3), bio_feature_name = "region",
                                     batch_feature_name = "xyz", c("test1","test2","test3")))

})


