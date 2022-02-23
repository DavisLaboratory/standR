calDESeq2NormCount <- function(spe_object, log = TRUE){
  # get raw count
  count_df <- SummarizedExperiment::assay(spe_object,1)

  # compute geometric mean
  loggeomeans <- rowMeans(log(count_df))

  # compute size factor
  sf <- apply(count_df, 2, function(x) {
    exp(stats::median((log(x) - loggeomeans)[is.finite(loggeomeans) & x > 0]))
  })

  # normalise by size factor
  norm_count <- t(t(count_df)/sf)

  if(isTRUE(log)){
    lognorm_count <- log2(norm_count + 1)
  } else {
    lognorm_count <- norm_count
  }

  return(lognorm_count)
}


rpkm2tpm <- function(x){
  tpm <- t(t(x)/colSums(x)) * 1e6
  return(tpm)
}

#' Perform normalisation to GeoMX data
#'
#' @param spe_object A spatial experiment object.
#' @param method Normalisation method to use. Options: TMM, RPKM, TPM, CPM, upperquartile, deseq2norm
#' @param log Log-transformed or not.
#'
#' @return A spatial experiment object.
#' @export
#'
geomx_normalize <- function(spe_object, method = "TMM", log = TRUE){

  if(!(method %in% c("TMM","RPKM","TPM","CPM","upperquartile","deseq2norm"))){
    stop("Please make sure method mathced one of the following strings: TMM,RPKM,TPM,CPM,upperquartile,deseq2")
  }

  . <- NULL
  # spe object to dgelist
  spe <- spe_object

  y <- edgeR::SE2DGEList(spe)

  ## TMM

  if(method == "TMM"){
    if(isTRUE(log)){
      spe@assays@data$logcounts <- edgeR::calcNormFactors(y) %>%
        edgeR::cpm(., log = TRUE)
    } else {
      spe@assays@data$logcounts <- edgeR::calcNormFactors(y) %>%
        edgeR::cpm(.)
    }
  }

  ## RPKM

  if(method == "RPKM"){
    if(isTRUE(log)){
      spe@assays@data$logcounts <- edgeR::rpkm(y, log = TRUE)
    } else {
      spe@assays@data$logcounts <- edgeR::rpkm(y, log = FALSE)
    }
  }


  ## TPM

  if(method == "TPM"){
    if(isTRUE(log)){
      spe@assays@data$logcounts <- log2(rpkm2tpm(2^edgeR::rpkm(y)))
    } else {
      spe@assays@data$logcounts <- rpkm2tpm(2^edgeR::rpkm(y))
    }
  }


  ## upper quantile

  if(method == "upperquartile"){
    if(isTRUE(log)){
      spe@assays@data$logcounts <- edgeR::calcNormFactors(y, method = "upperquartile") %>%
        edgeR::cpm(., log = TRUE)
    } else {
      spe@assays@data$logcounts <- edgeR::calcNormFactors(y, method = "upperquartile") %>%
        edgeR::cpm(.)
    }
  }


  ## normalisation methods from deseq2

  if(method == "deseq2norm"){
    if(isTRUE(log)){
      spe@assays@data$logcounts <- calDESeq2NormCount(spe, log = TRUE)
    } else {
      spe@assays@data$logcounts <- calDESeq2NormCount(spe, log = FALSE)
    }
  }

  return(spe)
}
