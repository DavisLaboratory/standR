#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames

# a handy function opposite to %in%
`%ni` <- Negate(`%in%`)


#' Import GeoMX DSP data into a saptial experiment object from file paths
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
#' @examples
#' dirPath <- ""
#' countFile <- system.file("extdata", "dkd_subset_TargetCountMatrix.txt", package = "standR")
#' sampleAnnoFile <- system.file("extdata", "dkd_subset_Sample_Annotations.txt", package = "standR")
#'
#' spe <- geomx_import_from_path(dirPath, countFile, sampleAnnoFile, hasNegProbe = FALSE)
#'
geomx_import_from_path <- function(dirPath, countFile, sampleAnnoFile, featureAnnoFile = NA,
                              hasNegProbe = TRUE, NegProbeName = "NegProbe-WTX",
                              colnames.as.rownames = c("TargetName","SegmentDisplayName","TargetName"),
                              coord.colnames = c("ROICoordinateX", "ROICoordinateY")){
  stopifnot(file.exists(file.path(dirPath, countFile)))
  stopifnot(file.exists(file.path(dirPath, sampleAnnoFile)))
  if(!is.na(featureAnnoFile)){
    stopifnot(file.exists(file.path(dirPath, featureAnnoFile)))
  }
  stopifnot(is.character(NegProbeName))
  stopifnot(length(colnames.as.rownames) == 3)
  stopifnot(length(coord.colnames) == 2)
  spe <- suppressMessages(geomx_import_fun(dirPath, countFile, sampleAnnoFile, featureAnnoFile,
                                       hasNegProbe, NegProbeName, colnames.as.rownames, coord.colnames))
  return(spe)
}


# the importing function itself
geomx_import_fun <- function(dirPath, countFile, sampleAnnoFile, featureAnnoFile,
                             hasNegProbe, NegProbeName,
                             colnames.as.rownames,
                             coord.colnames){

  . = NULL
  datalist <- list()

  # remove the NegProbe gene from the count matrix and save it in the metadata
  if(hasNegProbe == TRUE){
    countdata <- readr::read_tsv(file.path(dirPath, countFile))
    # raw count without negprobes
    stopifnot(colnames.as.rownames[1] %in% colnames(countdata)) # make sure count data have the gene column name as pre-defined, such as TargetName.
    stopifnot(NegProbeName %in% as.matrix(countdata[,colnames.as.rownames[1]])) # make sure the name of negprobe is in the gene column of count data.
    countdata_filtered <- countdata %>% # filter the count data, remove the negprobe.
      dplyr::filter(!!rlang::sym(colnames.as.rownames[1]) != NegProbeName) %>%
      tibble::column_to_rownames(colnames.as.rownames[1])

    # gene meta without negprobes
    if(!is.na(featureAnnoFile)){
      genemeta <- readr::read_tsv(file.path(dirPath, featureAnnoFile))
      stopifnot(colnames.as.rownames[3] %in% colnames(genemeta)) # make sure column name is there in the gene meta.
      genemeta_filtered <- genemeta %>%
        tibble::column_to_rownames(colnames.as.rownames[3]) %>%
        .[rownames(countdata_filtered),] # arrange the gene meta, as the same order as count table.
    } else {
      genemeta_filtered <- data.frame(Type = rep("gene",nrow(countdata_filtered)))
    }

    # sample meta
    samplemeta <- readr::read_tsv(file.path(dirPath, sampleAnnoFile))
    stopifnot(colnames.as.rownames[2] %in% colnames(samplemeta)) # make sure column name is there.
    samplemeta_filtered <- samplemeta %>%
      tibble::column_to_rownames(colnames.as.rownames[2]) %>%
      .[colnames(countdata_filtered),] # arrange acoording to count table.

    # negprobe raw count
    negprobecount <- countdata %>%
      dplyr::filter(!!rlang::sym(colnames.as.rownames[1]) == NegProbeName) %>%
      tibble::column_to_rownames(colnames.as.rownames[1])

    # logCPM count
    countdata_filtered_lcpm <- countdata_filtered %>%
      edgeR::cpm(log = TRUE)

    # output spe
    spe <- SpatialExperiment::SpatialExperiment(
      assay = list(counts = countdata_filtered,
                   logcounts = countdata_filtered_lcpm),
      colData = samplemeta_filtered,
      rowData = genemeta_filtered,
      metadata = list(NegProbes = negprobecount),
      spatialCoords = samplemeta_filtered %>%
        dplyr::select(coord.colnames) %>%
        as.matrix())
  } else {
    # it doesn't remove the NegProbe genes, leave them in the count matrix
    # raw count
    countdata <- readr::read_tsv(file.path(dirPath, countFile))
    stopifnot(colnames.as.rownames[1] %in% colnames(countdata))
    countdata <- countdata %>%
      tibble::column_to_rownames(colnames.as.rownames[1])
    # gene meta check
    if(!is.na(featureAnnoFile)){
      gene_meta <-  readr::read_tsv(file.path(dirPath, featureAnnoFile))%>%
        tibble::column_to_rownames(colnames.as.rownames[3]) %>%
        .[rownames(countdata),]
    } else {
      gene_meta <- data.frame(Type = rep("gene",nrow(countdata)))
    }

    # sample meta
    samplemeta <- readr::read_tsv(file.path(dirPath, sampleAnnoFile))
    stopifnot(colnames.as.rownames[2] %in% colnames(samplemeta))
    samplemeta<- samplemeta %>%
      tibble::column_to_rownames(colnames.as.rownames[2]) %>%
      .[colnames(countdata),]
    # negprobe raw count
    countdata_lcpm <- countdata %>%
      edgeR::cpm(log = TRUE)


    spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = countdata,
                                                              logcounts = countdata_lcpm),
                                                colData = samplemeta,
                                                rowData = gene_meta,
                                                spatialCoords = samplemeta %>%
                                                  dplyr::select(c(coord.colnames)) %>%
                                                  as.matrix())
  }
  return(spe)
}



#' Import GeoMX DSP data into a saptial experiment object from DGEList object
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
#' Counts <- matrix(rnbinom(ng*ns,mu=5,size=2),ng,ns)
#' rownames(Counts) <- 1:ng
#' y <- edgeR::DGEList(counts=Counts, group=rep(1:2,each=5))
#'
#' # transfer into spatial experiment object
#' coords <- matrix(rnorm(2*ns),10,2)
#' spe <- geomx_import_from_dge(dge_object = y, spatialCoord = coords)
#' spe
#'
geomx_import_from_dge <- function(dge_object, spatialCoord = NULL){
  spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = dge_object$counts,
                                                            logcounts = edgeR::cpm(dge_object, log = TRUE)),
                                              colData = dge_object$samples,
                                              rowData = dge_object$genes,
                                              spatialCoords = spatialCoord)
}

