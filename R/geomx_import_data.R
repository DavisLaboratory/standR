#' Import GeoMX DSP data into a saptial experiment object from file paths
#'
#' @param countFile tsv file or a dataframe object. Count matrix, with samples in columns and features/genes in rows. The first column is gene names/ids.
#' @param sampleAnnoFile tsv file or a dataframe object. Sample annotations.
#' @param featureAnnoFile tsv file or a dataframe object. Feature/Gene annotations.
#' @param hasNegProbe Logical. Default is TRUE, indicating there are negative probe genes in the data.
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
#' spe <- readGeoMx(countFile, sampleAnnoFile, hasNegProbe = FALSE)
#'
readGeoMx <- function(countFile, sampleAnnoFile, featureAnnoFile = NA,
                      hasNegProbe = TRUE, NegProbeName = "NegProbe-WTX",
                      colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),
                      coord.colnames = c("ROICoordinateX", "ROICoordinateY")) {
  stopifnot(is.character(NegProbeName))
  stopifnot(length(colnames.as.rownames) == 3)
  stopifnot(length(coord.colnames) == 2)
  spe <- geomx_import_fun(
    countFile, sampleAnnoFile, featureAnnoFile,
    hasNegProbe, NegProbeName, colnames.as.rownames, coord.colnames
  )
  return(spe)
}


# the importing function itself
geomx_import_fun <- function(countFile, sampleAnnoFile, featureAnnoFile,
                             hasNegProbe, NegProbeName,
                             colnames.as.rownames,
                             coord.colnames) {

  # remove the NegProbe gene from the count matrix and save it in the metadata
  if (hasNegProbe) {
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
    if (!NegProbeName %in% as.matrix(countdata[, colnames.as.rownames[1]])) {
      stop("NegProbeName is not found in your count file.")
    }

    # filter the count data, remove the negprobe.
    countdata_filtered0 <- countdata[countdata[, colnames.as.rownames[1]] != NegProbeName, ]
    countdata_filtered <- countdata_filtered0[, !colnames(countdata_filtered0) %in%
      colnames.as.rownames[1]]
    rownames(countdata_filtered) <- as.vector(as.matrix(countdata_filtered0[, colnames.as.rownames[1]]))


    # gene meta without negprobes
    if (!is.na(featureAnnoFile)) {
      if(is.data.frame(featureAnnoFile)){
        genemeta <- featureAnnoFile
      } else {
        genemeta <- as.data.frame(readr::read_tsv(featureAnnoFile), optional = TRUE)
      }

      stopifnot(colnames.as.rownames[3] %in% colnames(genemeta)) # make sure column name is there in the gene meta.

      genemeta_filtered0 <- genemeta[genemeta[, colnames.as.rownames[3]] != NegProbeName, ]
      genemeta_filtered <- genemeta_filtered0[, !colnames(genemeta_filtered0) %in%
        colnames.as.rownames[3]]
      rownames(genemeta_filtered) <- genemeta_filtered0[, colnames.as.rownames[3]]
      genemeta_filtered <- genemeta_filtered[rownames(countdata_filtered), ]
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
      colnames.as.rownames[2]]
    rownames(samplemeta_filtered) <- samplemeta[, colnames.as.rownames[2]]
    samplemeta_filtered <- samplemeta_filtered[colnames(countdata_filtered), ]
    # arrange according to count table.

    # negprobe raw count
    negprobecount <- countdata[countdata[, colnames.as.rownames[1]] == NegProbeName, ]
    negprobecount <- negprobecount[, !colnames(negprobecount) %in%
      colnames.as.rownames[1]]
    rownames(negprobecount) <- NegProbeName


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
    countdata <- as.data.frame(readr::read_tsv(countFile), optional = TRUE)

    stopifnot(colnames.as.rownames[1] %in% colnames(countdata))

    countdata <- tibble::column_to_rownames(countdata, colnames.as.rownames[1])
    # gene meta check
    if (!is.na(featureAnnoFile)) {
      gene_meta <- as.data.frame(readr::read_tsv(featureAnnoFile),
        optional = TRUE
      )
      gene_meta <- tibble::column_to_rownames(gene_meta, colnames.as.rownames[3])
      gene_meta <- gene_meta[rownames(countdata), ]
    } else {
      gene_meta <- data.frame(Type = rep("gene", nrow(countdata)))
    }

    # sample meta
    samplemeta <- as.data.frame(readr::read_tsv(sampleAnnoFile),
      optional = TRUE
    )

    stopifnot(colnames.as.rownames[2] %in% colnames(samplemeta))

    samplemeta <- tibble::column_to_rownames(samplemeta, colnames.as.rownames[2])
    samplemeta <- samplemeta[colnames(countdata), ]

    # negprobe raw count
    countdata_lcpm <- edgeR::cpm(countdata, log = TRUE)

    spe <- SpatialExperiment::SpatialExperiment(
      assays = list(
        counts = countdata,
        logcounts = countdata_lcpm
      ),
      colData = samplemeta,
      rowData = gene_meta,
      spatialCoords = as.matrix(samplemeta[, coord.colnames])
    )
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


utils::globalVariables(c("."))
