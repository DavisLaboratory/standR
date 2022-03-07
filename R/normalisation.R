# DEseq2 normalisation
calNormCount <- function(spe_object, log = TRUE){
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



rpkm2tpm <- function(x) {
  colSumMat <- edgeR::expandAsMatrix(colSums(x, na.rm = TRUE),
                                     byrow = TRUE,
                                     dim = dim(x))
  tpm <- x / colSumMat * 1e6
  return(tpm)
}


#' Perform normalisation to GeoMX data
#'
#' @param spe_object A spatial experiment object.
#' @param method Normalisation method to use. Options: TMM, RPKM, TPM, CPM, upperquartile, sizefactor. RPKM and TPM require gene length information, which should be added into rowData(spe). Note that TMM here is TMM + CPM.
#' @param log Log-transformed or not.
#'
#' @return A spatial experiment object, with the second assay being the normalised count matrix.
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#'
#' spe_tmm <- geomxNorm(dkd_spe_subset, method = "TMM")
#' head(SummarizedExperiment::assay(spe_tmm, 2))
#' spe_upq <- geomxNorm(dkd_spe_subset, method = "upperquartile")
#' head(SummarizedExperiment::assay(spe_upq, 2))
#' spe_deseqnorm <- geomxNorm(dkd_spe_subset, method = "sizefactor")
#' head(SummarizedExperiment::assay(spe_deseqnorm, 2))
#'
geomxNorm <- function(spe_object, method = "TMM", log = TRUE){

  if(!(method %in% c("TMM","RPKM","TPM","CPM","upperquartile","sizefactor"))){
    stop("Please make sure method mathced one of the following strings: TMM,RPKM,TPM,CPM,upperquartile,sizefactor")
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
      spe@assays@data$logcounts <- edgeR::rpkm(y, log = TRUE, prior.count = 0)
    } else {
      spe@assays@data$logcounts <- edgeR::rpkm(y, log = FALSE, prior.count = 0)
    }
  }


  ## TPM

  if(method == "TPM"){
    if(isTRUE(log)){
      spe@assays@data$logcounts <- log(rpkm2tpm(edgeR::rpkm(y))+1, 2)
    } else {
      spe@assays@data$logcounts <- rpkm2tpm(edgeR::rpkm(y))
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


  ## calculating size factor based on geomean

  if(method == "sizefactor"){
    if(isTRUE(log)){
      spe@assays@data$logcounts <- calNormCount(spe, log = TRUE)
    } else {
      spe@assays@data$logcounts <- calNormCount(spe, log = FALSE)
    }
  }

  return(spe)
}

