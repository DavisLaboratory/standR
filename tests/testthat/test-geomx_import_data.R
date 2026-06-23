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

test_that("NanoStringGeoMxSet-like objects can be converted to SpatialExperiment", {
  library(Biobase)

  expr <- matrix(
    c(10, 20, 30, 1, 5, 6, 7, 2),
    nrow = 4,
    dimnames = list(
      c("ProbeA1", "ProbeA2", "ProbeB1", "NegProbe1"),
      c("dcc1", "dcc2")
    )
  )
  pheno <- data.frame(
    group = c("A", "B"),
    ROICoordinateX = c(1, 2),
    ROICoordinateY = c(3, 4),
    row.names = c("dcc1", "dcc2")
  )
  feature <- data.frame(
    TargetName = c("GeneA", "GeneA", "GeneB", "NegProbe-WTX"),
    CodeClass = c("Endogenous", "Endogenous", "Endogenous", "Negative"),
    row.names = rownames(expr)
  )
  geomx_set <- Biobase::ExpressionSet(
    assayData = expr,
    phenoData = Biobase::AnnotatedDataFrame(pheno),
    featureData = Biobase::AnnotatedDataFrame(feature)
  )

  spe <- readGeoMxFromNanoStringGeoMxSet(geomx_set, rmNegProbe = TRUE)

  expect_s4_class(spe, "SpatialExperiment")
  expect_equal(rownames(spe), c("GeneA", "GeneB"))
  expect_equal(colnames(spe), c("dcc1", "dcc2"))
  expect_equal(
    unlist(SummarizedExperiment::assay(spe, "counts")["GeneA", ]),
    c(dcc1 = 30, dcc2 = 11)
  )
  expect_equal(
    rownames(S4Vectors::metadata(spe)$NegProbes),
    "NegProbe-WTX"
  )
})
