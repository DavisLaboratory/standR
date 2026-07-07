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

test_that("NanoStringGeoMxSet-like import supports one-column feature metadata", {
  library(Biobase)

  expr <- matrix(
    c(10, 1),
    nrow = 2,
    dimnames = list(c("GeneA", "NegProbe-WTX"), "ROI1")
  )
  pheno <- data.frame(
    ROICoordinateX = 1,
    ROICoordinateY = 1,
    row.names = "ROI1"
  )
  feature <- data.frame(
    TargetName = rownames(expr),
    row.names = rownames(expr)
  )
  geomx_set <- Biobase::ExpressionSet(
    assayData = expr,
    phenoData = Biobase::AnnotatedDataFrame(pheno),
    featureData = Biobase::AnnotatedDataFrame(feature)
  )

  spe <- readGeoMxFromNanoStringGeoMxSet(geomx_set)

  expect_equal(dim(spe), c(1, 1))
  expect_equal(rownames(spe), "GeneA")
  expect_equal(colnames(spe), "ROI1")
})

test_that("DCC and PKC files can be imported without GeomxTools", {
  dcc1 <- tempfile(fileext = ".dcc")
  dcc2 <- tempfile(fileext = ".dcc")
  pkc <- tempfile(fileext = ".pkc")

  writeLines(c(
    "<Header>",
    "FileVersion,0.01",
    "</Header>",
    "",
    "<Scan_Attributes>",
    "ID,ROI1",
    "Plate_ID,PlateA",
    "Well,A01",
    "</Scan_Attributes>",
    "",
    "<NGS_Processing_Attributes>",
    "Raw,10",
    "</NGS_Processing_Attributes>",
    "",
    "<Code_Summary>",
    "RTS001,10",
    "RTS002,20",
    "RTSNEG,1",
    "</Code_Summary>"
  ), dcc1)
  writeLines(c(
    "<Header>",
    "FileVersion,0.01",
    "</Header>",
    "",
    "<Scan_Attributes>",
    "ID,ROI2",
    "Plate_ID,PlateA",
    "Well,A02",
    "</Scan_Attributes>",
    "",
    "<NGS_Processing_Attributes>",
    "Raw,20",
    "</NGS_Processing_Attributes>",
    "",
    "<Code_Summary>",
    "RTS001,5",
    "RTS002,6",
    "RTSNEG,2",
    "</Code_Summary>"
  ), dcc2)
  writeLines(c(
    "{",
    '  "Name": "Tiny PKC",',
    '  "Version": 1.0,',
    '  "AnalyteType": "RNA",',
    '  "Targets": [',
    '    {"DisplayName": "GeneA", "CodeClass": "Endogenous01", "Probes": [',
    '      {"RTS_ID": "RTS001", "ProbeID": "P1", "GeneID": ["1"], "SystematicName": ["GeneA"]},',
    '      {"RTS_ID": "RTS002", "ProbeID": "P2", "GeneID": ["1"], "SystematicName": ["GeneA"]}',
    "    ]},",
    '    {"DisplayName": "NegProbe-WTX", "CodeClass": "Negative", "Probes": [',
    '      {"RTS_ID": "RTSNEG", "ProbeID": "N1", "GeneID": [], "SystematicName": []}',
    "    ]}",
    "  ]",
    "}"
  ), pkc)

  sampleAnno <- data.frame(
    Sample_ID = basename(c(dcc1, dcc2)),
    SegmentDisplayName = c("ROI1", "ROI2"),
    ROICoordinateX = c(1, 2),
    ROICoordinateY = c(3, 4)
  )

  spe <- readGeoMxFromDcc(c(dcc1, dcc2), pkcFiles = pkc,
                          sampleAnnoFile = sampleAnno)

  expect_s4_class(spe, "SpatialExperiment")
  expect_equal(rownames(spe), "GeneA")
  expect_equal(colnames(spe), c("ROI1", "ROI2"))
  expect_equal(
    unlist(SummarizedExperiment::assay(spe, "counts")["GeneA", ]),
    c(ROI1 = 30, ROI2 = 11)
  )
  expect_equal(
    rownames(S4Vectors::metadata(spe)$NegProbes),
    "NegProbe-WTX"
  )
})

test_that("GeomxTools is not required as a suggested package", {
  description <- read.dcf(test_path("../../DESCRIPTION"))

  expect_false(grepl("GeomxTools", description[1, "Suggests"], fixed = TRUE))
})

test_that("DCC import supports one sample after negative probe removal", {
  dcc <- tempfile(fileext = ".dcc")
  pkc <- tempfile(fileext = ".pkc")

  writeLines(c(
    "<Header>",
    "FileVersion,0.01",
    "</Header>",
    "<Scan_Attributes>",
    "ID,ROI1",
    "</Scan_Attributes>",
    "<Code_Summary>",
    "RTS001,10",
    "RTSNEG,1",
    "</Code_Summary>"
  ), dcc)
  writeLines(c(
    "{",
    '  "Name": "Tiny PKC", "Version": 1.0, "AnalyteType": "RNA",',
    '  "Targets": [',
    '    {"DisplayName": "GeneA", "CodeClass": "Endogenous", "Probes": [',
    '      {"RTS_ID": "RTS001", "ProbeID": "P1"}',
    "    ]},",
    '    {"DisplayName": "NegProbe-WTX", "CodeClass": "Negative", "Probes": [',
    '      {"RTS_ID": "RTSNEG", "ProbeID": "N1"}',
    "    ]}",
    "  ]",
    "}"
  ), pkc)

  sampleAnno <- data.frame(
    Sample_ID = basename(dcc),
    SegmentDisplayName = "ROI1",
    ROICoordinateX = 1,
    ROICoordinateY = 1
  )

  spe <- readGeoMxFromDcc(dcc, pkcFiles = pkc, sampleAnnoFile = sampleAnno)

  expect_equal(dim(spe), c(1, 1))
  expect_equal(colnames(spe), "ROI1")
  expect_equal(
    as.numeric(SummarizedExperiment::assay(spe, "counts")["GeneA", ]),
    10
  )
})

test_that("DCC import creates default sample annotation when none is supplied", {
  dcc <- tempfile(fileext = ".dcc")

  writeLines(c(
    "<Header>",
    "FileVersion,0.01",
    "</Header>",
    "<Scan_Attributes>",
    "ID,ROI1",
    "</Scan_Attributes>",
    "<Code_Summary>",
    "RTS001,10",
    "</Code_Summary>"
  ), dcc)

  spe <- readGeoMxFromDcc(dcc, rmNegProbe = FALSE)

  expect_equal(dim(spe), c(1, 1))
  expect_equal(colnames(spe), basename(dcc))
})
