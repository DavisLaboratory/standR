test_that("Testing normalisation function", {
  data("dkd_spe_subset")


  # check output dimensions
  spe_tmm <- geomxNorm(dkd_spe_subset, method = "TMM")
  expect_equal(dim(SummarizedExperiment::assay(spe_tmm, 2)), dim(dkd_spe_subset))
  expect_equal(typeof(SummarizedExperiment::assay(spe_tmm, 2)), "double")
  spe_upq <- geomxNorm(dkd_spe_subset, method = "upperquartile")
  expect_equal(dim(SummarizedExperiment::assay(spe_upq, 2)), dim(dkd_spe_subset))
  expect_equal(typeof(SummarizedExperiment::assay(spe_upq, 2)), "double")
  spe_deseqnorm <- geomxNorm(dkd_spe_subset, method = "sizefactor")
  expect_equal(dim(SummarizedExperiment::assay(spe_deseqnorm, 2)), dim(dkd_spe_subset))
  expect_equal(typeof(SummarizedExperiment::assay(spe_deseqnorm, 2)), "double")

  # check output dimensions of tpm and rpkm
  spe <- dkd_spe_subset
  rowData(spe)$GeneLength <- sample(100:10000, nrow(spe))
  spe_rpkm <- geomxNorm(spe, method = "RPKM")
  expect_equal(dim(SummarizedExperiment::assay(spe_rpkm, 2)), dim(spe))
  expect_equal(typeof(SummarizedExperiment::assay(spe_rpkm, 2)), "double")
  spe_tpm <- geomxNorm(spe, method = "TPM")
  expect_equal(dim(SummarizedExperiment::assay(spe_tpm, 2)), dim(spe))
  expect_equal(typeof(SummarizedExperiment::assay(spe_tpm, 2)), "double")

  # check weird parameter input
  expect_error(geomxNorm(dkd_spe_subset, method = "x"))

})
