#' Import GeoMX DSP data into a saptial experiment object from file paths
#'
#' @param countFile tsv file or a dataframe object. Count matrix, with samples in columns and features/genes in rows. The first column is gene names/ids.
#' @param sampleAnnoFile tsv file or a dataframe object. Sample annotations.
#' @param featureAnnoFile tsv file or a dataframe object. Feature/Gene annotations.
#' @param rmNegProbe Logical. Default is TRUE, indicating there are negative probe genes in the data.
#' @param NegProbeName Character. Name of negative probe genes, default is NegProbe-WTX.
#' @param colnames.as.rownames Vector of characters, length of 3. Column names used to capture gene names, sample names and gene names in countFile, sampleAnnoFile and featureAnnoFile, respectively.
#' @param coord.colnames Vector of characters, length of 2. Column names used to capture ROI coordinates.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
#' library(ExperimentHub)
#' 
#' eh <- ExperimentHub()
#' query(eh, "standR")
#' countFile <- eh[["EH7364"]]
#' sampleAnnoFile <- eh[["EH7365"]]
#'
#' spe <- readGeoMx(countFile, sampleAnnoFile, rmNegProbe = FALSE)
#'
readGeoMx <- function(countFile, sampleAnnoFile, featureAnnoFile = NA,
                      rmNegProbe = TRUE, NegProbeName = "NegProbe-WTX",
                      colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),
                      coord.colnames = c("ROICoordinateX", "ROICoordinateY")) {
  stopifnot(is.character(NegProbeName))
  stopifnot(length(colnames.as.rownames) == 3)
  stopifnot(length(coord.colnames) == 2)
  spe <- geomx_import_fun(
    countFile, sampleAnnoFile, featureAnnoFile,
    rmNegProbe, NegProbeName, colnames.as.rownames, coord.colnames
  )
  return(spe)
}


#' Import GeoMX DSP data from a NanoStringGeoMxSet-like object
#'
#' Converts a `Biobase::ExpressionSet`-derived object into a
#' `SpatialExperiment` object.
#'
#' @param geomxSet A `Biobase::ExpressionSet`-derived object.
#' @param assay2use Name of the assay data element to use. Default is `exprs`.
#' @param sampleIDCol Column in the sample annotation used as sample IDs. If the
#'   column is absent it is created from the assay column names.
#' @param featureIDCol Column in the feature annotation used as feature IDs. If
#'   duplicated, counts are summed by this column.
#' @param rmNegProbe Logical. Default is TRUE, indicating negative probe genes
#'   are removed from the assay and stored in metadata.
#' @param NegProbeName Character. Name of negative probe genes.
#' @param coord.colnames Vector of characters, length of 2. Column names used to
#'   capture ROI coordinates.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
#' expr <- matrix(c(10, 1), nrow = 2,
#'                dimnames = list(c("GeneA", "NegProbe-WTX"), "ROI1"))
#' pheno <- data.frame(ROICoordinateX = 1, ROICoordinateY = 1, row.names = "ROI1")
#' feature <- data.frame(TargetName = rownames(expr), row.names = rownames(expr))
#' eset <- Biobase::ExpressionSet(
#'   assayData = expr,
#'   phenoData = Biobase::AnnotatedDataFrame(pheno),
#'   featureData = Biobase::AnnotatedDataFrame(feature)
#' )
#' spe <- readGeoMxFromNanoStringGeoMxSet(eset)
readGeoMxFromNanoStringGeoMxSet <- function(geomxSet, assay2use = "exprs",
                                            sampleIDCol = "SegmentDisplayName",
                                            featureIDCol = "TargetName",
                                            rmNegProbe = TRUE,
                                            NegProbeName = "NegProbe-WTX",
                                            coord.colnames = c("ROICoordinateX", "ROICoordinateY")) {
  if (!methods::is(geomxSet, "ExpressionSet")) {
    stop("geomxSet must be a Biobase::ExpressionSet-derived object.")
  }

  count_matrix <- Biobase::assayDataElement(geomxSet, assay2use)
  count_matrix <- as.matrix(count_matrix)
  samplemeta <- Biobase::pData(geomxSet)
  featuremeta <- Biobase::fData(geomxSet)

  if (!featureIDCol %in% colnames(featuremeta)) {
    featuremeta[[featureIDCol]] <- rownames(count_matrix)
  }
  feature_ids <- as.character(featuremeta[[featureIDCol]])

  if (anyDuplicated(feature_ids)) {
    count_matrix <- rowsum(count_matrix, group = feature_ids, reorder = FALSE)
    featuremeta <- featuremeta[!duplicated(feature_ids), , drop = FALSE]
    featuremeta[[featureIDCol]] <- unique(feature_ids)
    rownames(featuremeta) <- featuremeta[[featureIDCol]]
  } else {
    rownames(count_matrix) <- feature_ids
    rownames(featuremeta) <- feature_ids
  }

  countdata <- data.frame(
    stats::setNames(list(rownames(count_matrix)), featureIDCol),
    as.data.frame(count_matrix, check.names = FALSE),
    check.names = FALSE
  )

  if (!sampleIDCol %in% colnames(samplemeta)) {
    samplemeta[[sampleIDCol]] <- colnames(count_matrix)
  }
  if (!all(colnames(count_matrix) %in% samplemeta[[sampleIDCol]])) {
    stop("sampleIDCol must identify all columns in the selected assay.")
  }

  readGeoMx(
    countFile = countdata,
    sampleAnnoFile = samplemeta,
    featureAnnoFile = featuremeta,
    rmNegProbe = rmNegProbe,
    NegProbeName = NegProbeName,
    colnames.as.rownames = c(featureIDCol, sampleIDCol, featureIDCol),
    coord.colnames = coord.colnames
  )
}


#' Import GeoMX DSP data from DCC files
#'
#' Reads NanoString GeoMX DCC count files directly and optionally uses PKC files
#' to map probe RTS IDs to target names. Multiple probes targeting the same gene
#' are summed before creating the `SpatialExperiment`.
#'
#' @param dccFiles Character vector of DCC file paths.
#' @param pkcFiles Optional character vector of PKC file paths.
#' @param sampleAnnoFile Optional sample annotation file path or data.frame.
#'   Delimited files are read as CSV when the extension is `.csv`, otherwise as
#'   TSV.
#' @param sampleIDCol Column in the sample annotation used as final sample IDs.
#' @param dccColName Column in the sample annotation containing DCC file names.
#' @param rmNegProbe Logical. Default is TRUE, indicating negative probe genes
#'   are removed from the assay and stored in metadata.
#' @param NegProbeName Character. Name of negative probe genes. If `NULL`, names
#'   are inferred from PKC rows with negative code class.
#' @param coord.colnames Vector of characters, length of 2. Column names used to
#'   capture ROI coordinates.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
#' dccFile <- tempfile(fileext = ".dcc")
#' writeLines(c(
#'   "<Header>", "FileVersion,0.01", "</Header>",
#'   "<Scan_Attributes>", "ID,ROI1", "</Scan_Attributes>",
#'   "<Code_Summary>", "RTS001,10", "RTSNEG,1", "</Code_Summary>"
#' ), dccFile)
#' pkcFile <- tempfile(fileext = ".pkc")
#' writeLines(c(
#'   "{",
#'   '  "Name": "Tiny PKC", "Version": 1.0, "AnalyteType": "RNA",',
#'   '  "Targets": [',
#'   '    {"DisplayName": "GeneA", "CodeClass": "Endogenous", "Probes": [',
#'   '      {"RTS_ID": "RTS001", "ProbeID": "P1"}',
#'   "    ]},",
#'   '    {"DisplayName": "NegProbe-WTX", "CodeClass": "Negative", "Probes": [',
#'   '      {"RTS_ID": "RTSNEG", "ProbeID": "N1"}',
#'   "    ]}",
#'   "  ]",
#'   "}"
#' ), pkcFile)
#' sampleAnno <- data.frame(
#'   Sample_ID = basename(dccFile),
#'   SegmentDisplayName = "ROI1",
#'   ROICoordinateX = 1,
#'   ROICoordinateY = 1
#' )
#' spe <- readGeoMxFromDcc(dccFile, pkcFile, sampleAnno)
readGeoMxFromDcc <- function(dccFiles, pkcFiles = NULL, sampleAnnoFile = NULL,
                             sampleIDCol = "SegmentDisplayName",
                             dccColName = "Sample_ID",
                             rmNegProbe = TRUE, NegProbeName = NULL,
                             coord.colnames = c("ROICoordinateX", "ROICoordinateY")) {
  if (!is.character(dccFiles) || length(dccFiles) == 0) {
    stop("dccFiles must be a non-empty character vector.")
  }
  if (!all(file.exists(dccFiles))) {
    stop("All dccFiles must exist.")
  }
  if (!is.null(pkcFiles) && (!is.character(pkcFiles) || !all(file.exists(pkcFiles)))) {
    stop("pkcFiles must be NULL or existing file paths.")
  }
  stopifnot(length(coord.colnames) == 2)

  dcc_ids <- basename(dccFiles)
  dcc_data <- stats::setNames(lapply(dccFiles, .read_dcc_file), dcc_ids)
  counts <- .dcc_count_matrix(dcc_data)
  featuremeta <- .dcc_feature_metadata(rownames(counts), pkcFiles)
  featuremeta <- featuremeta[match(rownames(counts), featuremeta$RTS_ID), , drop = FALSE]

  target_names <- featuremeta$TargetName
  count_by_target <- rowsum(counts, group = target_names, reorder = FALSE)
  feature_by_target <- featuremeta[!duplicated(target_names), , drop = FALSE]
  feature_by_target$TargetName <- unique(target_names)
  rownames(feature_by_target) <- feature_by_target$TargetName

  sample_info <- .dcc_sample_metadata(
    dcc_data = dcc_data,
    sampleAnnoFile = sampleAnnoFile,
    sampleIDCol = sampleIDCol,
    dccColName = dccColName,
    coord.colnames = coord.colnames
  )
  colnames(count_by_target) <- sample_info$sample_names

  if (is.null(NegProbeName)) {
    neg_idx <- grepl("negative", feature_by_target$CodeClass, ignore.case = TRUE) |
      grepl("negative", feature_by_target$TargetName, ignore.case = TRUE)
    NegProbeName <- unique(feature_by_target$TargetName[neg_idx])
  }
  if (rmNegProbe && length(NegProbeName) == 0) {
    stop("Negative probes could not be inferred; supply NegProbeName or set rmNegProbe = FALSE.")
  }
  if (is.null(NegProbeName)) {
    NegProbeName <- character()
  }

  countdata <- data.frame(
    TargetName = rownames(count_by_target),
    as.data.frame(count_by_target, check.names = FALSE),
    check.names = FALSE
  )

  readGeoMx(
    countFile = countdata,
    sampleAnnoFile = sample_info$samplemeta,
    featureAnnoFile = feature_by_target,
    rmNegProbe = rmNegProbe,
    NegProbeName = NegProbeName,
    colnames.as.rownames = c("TargetName", sampleIDCol, "TargetName"),
    coord.colnames = coord.colnames
  )
}


# the importing function itself
geomx_import_fun <- function(countFile, sampleAnnoFile, featureAnnoFile,
                             rmNegProbe, NegProbeName,
                             colnames.as.rownames,
                             coord.colnames) {

  # remove the NegProbe gene from the count matrix and save it in the metadata
  if (rmNegProbe) {
    if(is.data.frame(countFile)){
      countdata <- countFile
    } else {
      countdata <- as.data.frame(readr::read_tsv(countFile), optional = TRUE)
    }

    # raw count without negprobes
    # make sure count data have the gene column name as pre-defined, such as TargetName.
    if (!colnames.as.rownames[1] %in% colnames(countdata)) {
      stop("colnames.as.rownames[1] is not in the column names of your count file.")
    }
    # make sure the name of negprobe is in the gene column of count data.
    if (!all(NegProbeName %in% as.matrix(countdata[, colnames.as.rownames[1]]))) {
      stop("NegProbeName is not found in your count file.")
    }

    # filter the count data, remove the negprobe.
    countdata_filtered0 <- countdata[!countdata[, colnames.as.rownames[1]] %in% NegProbeName, , drop = FALSE]
    countdata_filtered <- countdata_filtered0[, !colnames(countdata_filtered0) %in%
      colnames.as.rownames[1], drop = FALSE]
    rownames(countdata_filtered) <- as.vector(as.matrix(countdata_filtered0[, colnames.as.rownames[1]]))


    # gene meta without negprobes
    if (!all(is.na(featureAnnoFile))) {
      if(is.data.frame(featureAnnoFile)){
        genemeta <- featureAnnoFile
      } else {
        genemeta <- as.data.frame(readr::read_tsv(featureAnnoFile), optional = TRUE)
      }

      stopifnot(colnames.as.rownames[3] %in% colnames(genemeta)) # make sure column name is there in the gene meta.

      genemeta_filtered0 <- genemeta[!genemeta[, colnames.as.rownames[3]] %in% NegProbeName, , drop = FALSE]
      genemeta_filtered <- genemeta_filtered0[, !colnames(genemeta_filtered0) %in%
        colnames.as.rownames[3], drop = FALSE]
      rownames(genemeta_filtered) <- genemeta_filtered0[, colnames.as.rownames[3]]
      genemeta_filtered <- genemeta_filtered[rownames(countdata_filtered), , drop = FALSE]
      # arrange the gene meta, as the same order as count table.
    } else {
      genemeta_filtered <- data.frame(Type = rep("gene", nrow(countdata_filtered)))
    }

    # sample meta
    if(is.data.frame(sampleAnnoFile)){
      samplemeta <- sampleAnnoFile
    } else {
      samplemeta <- as.data.frame(readr::read_tsv(sampleAnnoFile), optional = TRUE)
    }

    stopifnot(colnames.as.rownames[2] %in% colnames(samplemeta)) # make sure column name is there.


    samplemeta_filtered <- samplemeta[, !colnames(samplemeta) %in%
      colnames.as.rownames[2], drop = FALSE]
    rownames(samplemeta_filtered) <- samplemeta[, colnames.as.rownames[2]]
    samplemeta_filtered <- samplemeta_filtered[colnames(countdata_filtered), , drop = FALSE]
    # arrange according to count table.
    # negprobe raw count
    negprobecount <- countdata[countdata[, colnames.as.rownames[1]] %in% 
                                 NegProbeName, , drop = FALSE]
    nprobename <- as.vector(as.matrix(negprobecount[,colnames.as.rownames[1]]))
    if(length(nprobename) != length(unique(nprobename))){
      nprobename <- paste0(nprobename,"_",seq(length(nprobename)))
    }
    negprobecount <- negprobecount[, !colnames(negprobecount) %in% 
                                     colnames.as.rownames[1], drop = FALSE]
    rownames(negprobecount) <- nprobename


    # logCPM count
    countdata_filtered_lcpm <- edgeR::cpm(countdata_filtered, log = TRUE)


    # output spe
    spe <- SpatialExperiment::SpatialExperiment(
      assay = list(
        counts = countdata_filtered,
        logcounts = countdata_filtered_lcpm
      ),
      colData = samplemeta_filtered,
      rowData = genemeta_filtered,
      metadata = list(NegProbes = negprobecount),
      spatialCoords = as.matrix(samplemeta_filtered[, coord.colnames])
    )
  } else {
    # it doesn't remove the NegProbe genes, leave them in the count matrix
    # raw count
    if (is.data.frame(countFile)) {
      countdata0 <- countFile
    }
    else {
      countdata0 <- as.data.frame(readr::read_tsv(countFile), 
                                 optional = TRUE)
    }
    stopifnot(colnames.as.rownames[1] %in% colnames(countdata0))
    
    countdata <- countdata0[, !colnames(countdata0) %in% colnames.as.rownames[1], drop = FALSE]
    rownames(countdata) <- as.vector(as.matrix(countdata0[, colnames.as.rownames[1]]))
    
    if (!all(is.na(featureAnnoFile))) {
      if (is.data.frame(featureAnnoFile)) {
        genemeta0 <- featureAnnoFile
      }
      else {
        genemeta0 <- as.data.frame(readr::read_tsv(featureAnnoFile), 
                                  optional = TRUE)
      }
      stopifnot(colnames.as.rownames[3] %in% colnames(genemeta0))
      
      genemeta <- genemeta0[, !colnames(genemeta0) %in% colnames.as.rownames[3], drop = FALSE]
      rownames(genemeta) <- genemeta0[, colnames.as.rownames[3]]
      genemeta <- genemeta[rownames(countdata), , drop = FALSE]
    }
    else {
      genemeta <- data.frame(Type = rep("gene", nrow(countdata)))
    }
    
    if (is.data.frame(sampleAnnoFile)) {
      samplemeta0 <- sampleAnnoFile
    }
    else {
      samplemeta0 <- as.data.frame(readr::read_tsv(sampleAnnoFile), 
                                  optional = TRUE)
    }
    stopifnot(colnames.as.rownames[2] %in% colnames(samplemeta0))
    
    samplemeta <- samplemeta0[, !colnames(samplemeta0) %in% colnames.as.rownames[2], drop = FALSE]
    rownames(samplemeta) <- samplemeta0[, colnames.as.rownames[2]]
    samplemeta <- samplemeta[colnames(countdata), , drop = FALSE]
    
    countdata_lcpm <- edgeR::cpm(countdata, log = TRUE)
    spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = countdata, 
                                                              logcounts = countdata_lcpm), colData = samplemeta, 
                                                rowData = genemeta, spatialCoords = as.matrix(samplemeta[, 
                                                                                                          coord.colnames]))
  }
  return(spe)
}



#' Import GeoMX DSP data into a spatial experiment object from DGEList object
#'
#' @param dge_object a DGEList object (created using edgeR::DGEList).
#' @param spatialCoord a matrix with coordinates of samples, rowname must be cosistent with the colnames of dge_object.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#' @examples
#' # making a simple DGEList object
#' ng <- 1000
#' ns <- 10
#' Counts <- matrix(rnbinom(ng * ns, mu = 5, size = 2), ng, ns)
#' rownames(Counts) <- seq(ng)
#' y <- edgeR::DGEList(counts = Counts, group = rep(seq(2), each = 5))
#'
#' # transfer into spatial experiment object
#' coords <- matrix(rnorm(2 * ns), 10, 2)
#' spe <- readGeoMxFromDGE(dge_object = y, spatialCoord = coords)
#' spe
#'
readGeoMxFromDGE <- function(dge_object, spatialCoord = NULL) {
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(
      counts = dge_object$counts,
      logcounts = edgeR::cpm(dge_object, log = TRUE)
    ),
    colData = dge_object$samples,
    rowData = dge_object$genes,
    spatialCoords = spatialCoord
  )
}


.read_dcc_file <- function(file) {
  lines <- readLines(file, warn = FALSE)
  list(
    Header = .read_dcc_key_value_section(lines, "Header"),
    Scan_Attributes = .read_dcc_key_value_section(lines, "Scan_Attributes"),
    NGS_Processing_Attributes = .read_dcc_key_value_section(lines, "NGS_Processing_Attributes"),
    Code_Summary = .read_dcc_counts(lines, file)
  )
}

.read_dcc_section <- function(lines, section) {
  lines_trimmed <- trimws(lines)
  start <- match(paste0("<", section, ">"), lines_trimmed)
  end <- match(paste0("</", section, ">"), lines_trimmed)
  if (is.na(start) || is.na(end) || end <= start) {
    return(character())
  }
  section_lines <- lines[(start + 1):(end - 1)]
  section_lines[nzchar(trimws(section_lines))]
}

.read_dcc_key_value_section <- function(lines, section) {
  section_lines <- .read_dcc_section(lines, section)
  if (length(section_lines) == 0) {
    return(data.frame(check.names = FALSE))
  }

  values <- utils::read.csv(
    text = paste(section_lines, collapse = "\n"),
    header = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  keys <- make.unique(as.character(values[[1]]))
  out <- as.data.frame(as.list(as.character(values[[2]])), check.names = FALSE)
  colnames(out) <- keys
  out
}

.read_dcc_counts <- function(lines, file) {
  section_lines <- .read_dcc_section(lines, "Code_Summary")
  if (length(section_lines) == 0) {
    stop("No Code_Summary section found in ", basename(file), ".")
  }

  counts <- utils::read.csv(
    text = paste(section_lines, collapse = "\n"),
    header = FALSE,
    col.names = c("RTS_ID", "Count"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  counts$Count <- suppressWarnings(as.numeric(counts$Count))
  if (anyNA(counts$RTS_ID) || anyNA(counts$Count)) {
    stop("Invalid Code_Summary counts in ", basename(file), ".")
  }
  counts <- stats::aggregate(Count ~ RTS_ID, data = counts, sum)
  stats::setNames(counts$Count, counts$RTS_ID)
}

.dcc_count_matrix <- function(dcc_data) {
  feature_ids <- unique(unlist(lapply(dcc_data, function(x) names(x$Code_Summary)), use.names = FALSE))
  counts <- matrix(
    0,
    nrow = length(feature_ids),
    ncol = length(dcc_data),
    dimnames = list(feature_ids, names(dcc_data))
  )

  for (sample_id in names(dcc_data)) {
    sample_counts <- dcc_data[[sample_id]]$Code_Summary
    counts[names(sample_counts), sample_id] <- sample_counts
  }

  counts
}

.dcc_feature_metadata <- function(feature_ids, pkcFiles) {
  if (is.null(pkcFiles)) {
    return(data.frame(
      RTS_ID = feature_ids,
      TargetName = feature_ids,
      CodeClass = "Unknown",
      stringsAsFactors = FALSE
    ))
  }

  pkc_features <- do.call(rbind, lapply(pkcFiles, .read_pkc_features))
  pkc_features <- pkc_features[!duplicated(pkc_features$RTS_ID), , drop = FALSE]
  missing_ids <- setdiff(feature_ids, pkc_features$RTS_ID)
  if (length(missing_ids) > 0) {
    pkc_features <- rbind(
      pkc_features,
      data.frame(
        RTS_ID = missing_ids,
        TargetName = missing_ids,
        Module = NA_character_,
        CodeClass = "Unknown",
        ProbeID = NA_character_,
        GeneID = NA_character_,
        SystematicName = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }
  pkc_features
}

.read_pkc_features <- function(pkcFile) {
  pkc <- rjson::fromJSON(file = pkcFile)
  module <- sub("\\.pkc$", "", basename(pkcFile), ignore.case = TRUE)
  analyte <- pkc[["AnalyteType"]]

  rows <- lapply(pkc[["Targets"]], function(target) {
    code_class <- gsub("\\d+$", "", target[["CodeClass"]])
    target_name <- target[["DisplayName"]]
    probes <- target[["Probes"]]

    do.call(rbind, lapply(probes, function(probe) {
      rts_id <- if (identical(tolower(analyte), "protein")) {
        target[["RTS_ID"]]
      } else {
        probe[["RTS_ID"]]
      }
      data.frame(
        RTS_ID = rts_id,
        TargetName = target_name,
        Module = module,
        CodeClass = code_class,
        ProbeID = .pkc_collapse(probe[["ProbeID"]]),
        GeneID = .pkc_collapse(probe[["GeneID"]]),
        SystematicName = .pkc_collapse(probe[["SystematicName"]]),
        stringsAsFactors = FALSE
      )
    }))
  })

  do.call(rbind, rows)
}

.pkc_collapse <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NA_character_)
  }
  paste(unlist(value), collapse = ", ")
}

.dcc_sample_metadata <- function(dcc_data, sampleAnnoFile, sampleIDCol,
                                 dccColName, coord.colnames) {
  dcc_ids <- names(dcc_data)
  dcc_meta <- do.call(rbind, lapply(dcc_data, function(x) {
    .dcc_metadata_row(x)
  }))
  rownames(dcc_meta) <- dcc_ids

  if (is.null(sampleAnnoFile)) {
    samplemeta <- data.frame(dcc_ids, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(samplemeta) <- sampleIDCol
  } else {
    samplemeta <- .read_sample_annotation(sampleAnnoFile)
    if (dccColName %in% colnames(samplemeta)) {
      dcc_keys <- as.character(samplemeta[[dccColName]])
      dcc_keys <- ifelse(grepl("\\.dcc$", dcc_keys, ignore.case = TRUE),
                         dcc_keys, paste0(dcc_keys, ".dcc"))
      idx <- match(dcc_ids, dcc_keys)
      if (anyNA(idx)) {
        stop("sampleAnnoFile is missing annotations for: ",
             paste(dcc_ids[is.na(idx)], collapse = ", "))
      }
      samplemeta <- samplemeta[idx, , drop = FALSE]
    } else if (sampleIDCol %in% colnames(samplemeta)) {
      idx <- match(dcc_ids, as.character(samplemeta[[sampleIDCol]]))
      if (!anyNA(idx)) {
        samplemeta <- samplemeta[idx, , drop = FALSE]
      }
    } else {
      stop("sampleAnnoFile must contain ", dccColName, " or ", sampleIDCol, ".")
    }
    if (!sampleIDCol %in% colnames(samplemeta)) {
      samplemeta[[sampleIDCol]] <- dcc_ids
    }
  }

  sample_names <- as.character(samplemeta[[sampleIDCol]])
  if (anyNA(sample_names) || anyDuplicated(sample_names)) {
    stop(sampleIDCol, " must contain non-missing unique sample names.")
  }

  for (col in colnames(dcc_meta)) {
    if (!col %in% colnames(samplemeta)) {
      samplemeta[[col]] <- dcc_meta[[col]]
    }
  }
  for (coord_col in coord.colnames) {
    if (!coord_col %in% colnames(samplemeta)) {
      samplemeta[[coord_col]] <- NA_real_
    }
  }

  list(samplemeta = samplemeta, sample_names = sample_names)
}

.dcc_metadata_row <- function(dcc) {
  sections <- list(dcc$Header, dcc$Scan_Attributes, dcc$NGS_Processing_Attributes)
  sections <- sections[vapply(sections, ncol, integer(1)) > 0]
  if (length(sections) == 0) {
    return(data.frame(check.names = FALSE))
  }
  do.call(cbind, sections)
}

.read_sample_annotation <- function(sampleAnnoFile) {
  if (is.data.frame(sampleAnnoFile)) {
    return(as.data.frame(sampleAnnoFile, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (!is.character(sampleAnnoFile) || length(sampleAnnoFile) != 1 || !file.exists(sampleAnnoFile)) {
    stop("sampleAnnoFile must be NULL, a data.frame, or an existing delimited file.")
  }
  if (grepl("\\.csv$", sampleAnnoFile, ignore.case = TRUE)) {
    return(as.data.frame(readr::read_csv(sampleAnnoFile, show_col_types = FALSE), check.names = FALSE))
  }
  as.data.frame(readr::read_tsv(sampleAnnoFile, show_col_types = FALSE), check.names = FALSE)
}

utils::globalVariables(c("."))
