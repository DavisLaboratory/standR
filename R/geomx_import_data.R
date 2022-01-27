`%ni` <- Negate(`%in%`)


#' Import geomx data into spe object
#'
#' @param dirPath The folder path to all necessary tsv files.
#' @param countFile tsv file. Count matrix, with samples in columns and features/genes in rows. The first column is gene names/ids.
#' @param sampleAnnoFile tsv file. Sample annotations.
#' @param featureAnnoFile tsv file. Feature/Gene annotations.
#' @param hasNegProbe Logical. Indicating there are negative probe genes in the data.
#' @param NegProbeName Character. Name of negative probe genes.
#' @param colnames.as.rownames Vector of characters, length of 3. Column names used to capture gene names, sample names and gene names in countFile, sampleAnnoFile and featureAnnoFile, respectively.
#' @param coord.colnames Vector of characters, length of 2. Column names used to capture ROI coordinates.
#'
#' @return A Spatial Experiment object.
#' @export
#'
#'
geomx_import_fun <- function(dirPath, countFile, sampleAnnoFile, featureAnnoFile,
                         hasNegProbe, NegProbeName,
                         colnames.as.rownames,
                         coord.colnames){

  . = NULL
  datalist <- list()

  stopifnot(file.exists(dirPath))
  stopifnot(file.exists(file.path(dirPath, countFile)))
  stopifnot(file.exists(file.path(dirPath, sampleAnnoFile)))
  stopifnot(file.exists(file.path(dirPath, featureAnnoFile)))
  stopifnot(is.character(NegProbeName))
  stopifnot(length(colnames.as.rownames) == 3)
  stopifnot(length(coord.colnames) == 2)

  if(hasNegProbe == TRUE){
    c <- readr::read_tsv(file.path(dirPath, countFile))
    # raw count without negprobes
    datalist[[1]] <- c %>% dplyr::filter(!!rlang::sym(colnames.as.rownames[1]) != NegProbeName) %>%
      tibble::column_to_rownames(colnames.as.rownames[1])
    # gene meta without negprobes
    datalist[[2]] <- readr::read_tsv(file.path(dirPath, featureAnnoFile)) %>%
      dplyr::filter(!!rlang::sym(colnames.as.rownames[3]) != NegProbeName) %>%
      tibble::column_to_rownames(colnames.as.rownames[3])
    # sample meta
    datalist[[3]] <- readr::read_tsv(file.path(dirPath, sampleAnnoFile))
    datalist[[3]] <- datalist[[3]] %>%
      tibble::column_to_rownames(colnames.as.rownames[2]) %>%
      .[colnames(datalist[[1]]),]
    # negprobe raw count
    datalist[[4]] <- c %>% dplyr::filter(!!rlang::sym(colnames.as.rownames[1]) == NegProbeName) %>%
      tibble::column_to_rownames(colnames.as.rownames[1])
    # logCPM count
    datalist[[5]] <- datalist[[1]] %>%
      edgeR::cpm(log = TRUE)

    spe <- SpatialExperiment::SpatialExperiment(
      assay = list(counts = datalist[[1]],
                   logcounts = datalist[[5]]),
      colData = datalist[[3]],
      rowData = datalist[[2]] %>%
        .[rownames(datalist[[1]]),],
      metadata = list(NegProbes = datalist[[4]]),
      spatialCoords = datalist[[3]] %>%
        dplyr::select(coord.colnames) %>%
        as.matrix())
  } else {
    datalist <- list(countFile, featureAnnoFile, sampleAnnoFile) %>%
      lapply(function(x) readr::read_tsv(file.path(dirPath, x)))
    datalist[[1]] <- datalist[[1]] %>%
      tibble::column_to_rownames(colnames.as.rownames[1])
    datalist[[3]] <- datalist[[3]] %>%
      tibble::column_to_rownames(colnames.as.rownames[2]) %>%
      .[colnames(datalist[[1]]),]
    datalist[[4]] <- 1
    datalist[[5]] <- datalist[[1]] %>%
      edgeR::cpm(log = TRUE)


    spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = datalist[[1]],
                                           logcounts = datalist[[5]]),
                             colData = datalist[[3]],
                             rowData = datalist[[2]] %>%
                               tibble::column_to_rownames(colnames.as.rownames[3]) %>%
                               .[rownames(datalist[[1]]),],
                             spatialCoords = datalist[[3]] %>%
                               dplyr::select(c(coord.colnames)) %>%
                               as.matrix())
  }
  return(spe)
}


#' Import GeoMX DSP data into a saptial experiment object from path
#'
#' @param dirPath The folder path to all necessary tsv files.
#' @param countFile tsv file. Count matrix, with samples in columns and features/genes in rows. The first column is gene names/ids.
#' @param sampleAnnoFile tsv file. Sample annotations.
#' @param featureAnnoFile tsv file. Feature/Gene annotations.
#' @param hasNegProbe Logical. Default is TRUE, indicating there are negative probe genes in the data.
#' @param NegProbeName Character. Name of negative probe genes, default is NegProbe-WTX.
#' @param colnames.as.rownames Vector of characters, length of 3. Column names used to capture gene names, sample names and gene names in countFile, sampleAnnoFile and featureAnnoFile, respectively.
#' @param coord.colnames Vector of characters, length of 2. Column names used to capture ROI coordinates.
#'
#' @return A SpatialExperiment object.
#' @export
#'
#'
geomx_import_from_path <- function(dirPath, countFile, sampleAnnoFile, featureAnnoFile,
                              hasNegProbe = TRUE, NegProbeName = "NegProbe-WTX",
                              colnames.as.rownames = c("TargetName","SegmentDisplayName","TargetName"),
                              coord.colnames = c("ROICoordinateX", "ROICoordinateY")){
  stopifnot(file.exists(dirPath))
  stopifnot(file.exists(file.path(dirPath, countFile)))
  stopifnot(file.exists(file.path(dirPath, sampleAnnoFile)))
  stopifnot(file.exists(file.path(dirPath, featureAnnoFile)))
  stopifnot(is.character(NegProbeName))
  stopifnot(length(colnames.as.rownames) == 3)
  stopifnot(length(coord.colnames) == 2)
  spe <- suppressMessages(geomx_import_fun(dirPath, countFile, sampleAnnoFile, featureAnnoFile,
                                       hasNegProbe, NegProbeName, colnames.as.rownames, coord.colnames))
  return(spe)
}



#' Import GeoMX DSP data into a saptial experiment object from DGEList object
#'
#' @param dge_object a DGEList object.
#' @param spatialCoord a dataframe with coordinates of samples, rowname must be cosistent with the colnames of dge_object.
#'
#' @return A SpatialExperiment object.
#' @export
geomx_import_from_dge <- function(dge_object, spatialCoord = NULL){
  spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = dge_object$counts,
                                                            logcounts = edgeR::cpm(dge_object, log = TRUE)),
                                              colData = dge_object$samples,
                                              rowData = dge_object$genes,
                                              spatialCoords = spatialCoord)
}
