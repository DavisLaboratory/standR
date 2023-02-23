# DEseq2 normalization
calNormCount <- function(spe_object, log = TRUE) {
  # get raw count
  count_df <- assay(spe_object, 1)
  # compute geometric mean
  loggeomeans <- rowMeans(log(count_df))
  # compute size factor
  sf <- apply(count_df, 2, function(x) {
    exp(stats::median((log(x) - loggeomeans)[is.finite(loggeomeans) & x > 0]))
  })

  # normalise by size factor
  norm_count <- t(t(count_df) / sf)

  if (isTRUE(log)) {
    lognorm_count <- log2(norm_count + 1)
  } else {
    lognorm_count <- norm_count
  }
  return(list(lognorm_count, sf))
}



rpkm2tpm <- function(x) {
  colSumMat <- edgeR::expandAsMatrix(colSums(x, na.rm = TRUE),
    byrow = TRUE,
    dim = dim(x)
  )
  tpm <- x / colSumMat * 1e6
  return(tpm)
}


#' Perform normalization to GeoMX data
#'
#' @param spe_object A SpatialExperiment object.
#' @param method Normalization method to use. Options: TMM, RPKM, TPM, CPM, upperquartile, sizefactor. RPKM and TPM require gene length information, which should be added into rowData(spe). Note that TMM here is TMM + CPM.
#' @param log Log-transformed or not.
#'
#' @return A SpatialExperiment object, with the second assay being the normalized count matrix. The normalised count is stored in the assay slot called "logcounts" by default. 
#' With method TMM and sizefactor, the norm.factor will be saved in the metadata of the SpatialExperiment object.
#' @export
#'
#' @references Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139-140.
#' @references Love, M., Anders, S., & Huber, W. (2014). Differential analysis of count data–the DESeq2 package. Genome Biol, 15(550), 10-1186.
#'
#' @note The normalised count is not intended to be used directly for linear modelling. For linear modelling, it is better to include the normalized factors in the "norm.factors" column of the DGEList object.
#' 
#' @examples
#' data("dkd_spe_subset")
#'
#' spe_tmm <- geomxNorm(dkd_spe_subset, method = "TMM")
#' spe_upq <- geomxNorm(dkd_spe_subset, method = "upperquartile")
#' spe_deseqnorm <- geomxNorm(dkd_spe_subset, method = "sizefactor")
#'
geomxNorm <- function(spe_object, method = c(
                        "TMM", "RPKM", "TPM", "CPM",
                        "upperquartile", "sizefactor"
                      ),
                      log = TRUE) {

  if (!(method %in% c("TMM", "RPKM", "TPM", "CPM", "upperquartile", "sizefactor"))) {
    stop("Please make sure method matched one of the following strings: 
         TMM, 
         RPKM, 
         TPM, 
         CPM, 
         upperquartile, 
         sizefactor")
  }
  
  method <- match.arg(method)

  # spe object to dgelist
  spe <- spe_object

  y <- edgeR::SE2DGEList(spe)

  ## TMM

  if (method == "TMM") {
    S4Vectors::metadata(spe)$norm.factor <- edgeR::calcNormFactors(y)$samples$norm.factors
    if (isTRUE(log)) {
      assay(spe, "logcounts") <- edgeR::calcNormFactors(y) |>
        (\(.) edgeR::cpm(., log = TRUE))()
    } else {
      assay(spe, "logcounts") <- edgeR::calcNormFactors(y) |>
        (\(.) edgeR::cpm(.))()
    }
  }

  ## RPKM

  if (method == "RPKM") {
    if (isTRUE(log)) {
      assay(spe, "logcounts") <- edgeR::rpkm(y, log = TRUE, prior.count = 0)
    } else {
      assay(spe, "logcounts") <- edgeR::rpkm(y, log = FALSE, prior.count = 0)
    }
  }


  ## TPM

  if (method == "TPM") {
    if (isTRUE(log)) {
      assay(spe, "logcounts") <- log(rpkm2tpm(edgeR::rpkm(y)) + 1, 2)
    } else {
      assay(spe, "logcounts") <- rpkm2tpm(edgeR::rpkm(y))
    }
  }


  ## upper quantile

  if (method == "upperquartile") {
    S4Vectors::metadata(spe)$norm.factor <- edgeR::calcNormFactors(y, method = "upperquartile")$samples$norm.factors
    if (isTRUE(log)) {
      assay(spe, "logcounts") <- edgeR::calcNormFactors(y, method = "upperquartile") |>
        (\(.) edgeR::cpm(., log = TRUE))()
    } else {
      assay(spe, "logcounts") <- edgeR::calcNormFactors(y, method = "upperquartile") |>
        (\(.) edgeR::cpm(.))()
    }
  }


  ## calculating size factor based on geomean

  if (method == "sizefactor") {
    S4Vectors::metadata(spe)$norm.factor <- calNormCount(spe, log = TRUE)[[2]]
    if (isTRUE(log)) {
      assay(spe, "logcounts") <- calNormCount(spe, log = TRUE)[[1]]
    } else {
      assay(spe, "logcounts") <- calNormCount(spe, log = FALSE)[[1]]
    }
  }

  return(spe)
}

#' Transfer SpatialExperiment object into DGEList object for DE analysis
#'
#' @param spe SpatialExperiment object.
#' @param norm Prior normalization method to be used. Norm factor calculated from the previous normalization step to be used in the DGEList.
#'
#' @return A DGEList.
#' @export 
#'
#' @examples 
#' data("dkd_spe_subset")
#' 
#' dge <- spe2dge(dkd_spe_subset)
#' 
#' spe_tmm <- geomxNorm(dkd_spe_subset, method = "TMM")
#' dge <- spe2dge(spe_tmm, norm = "TMM")
#' 
spe2dge <- function(spe, norm = NULL){
  if (!is(spe, "SummarizedExperiment")) 
    stop("spe is not of the SummarizedExperiment class")
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) 
    stop("SummarizedExperiment package required but is not installed (or can't be loaded)")
  if (!("counts" %in% SummarizedExperiment::assayNames(spe))) 
    stop("spe doesn't contain counts")
  counts <- SummarizedExperiment::assay(spe, "counts")
  if (!is.null(rownames(spe))) 
    rownames(counts) <- rownames(spe)
  if (!is.null(colnames(spe))) 
    colnames(counts) <- colnames(spe)
  genes <- samples <- NULL
  if (ncol(SummarizedExperiment::colData(spe))) {
    samples <- as.data.frame(SummarizedExperiment::colData(spe))
  }
  if (is(SummarizedExperiment::rowRanges(spe), "GRanges")) 
    genes <- as.data.frame(SummarizedExperiment::rowRanges(spe))
  else if (ncol(SummarizedExperiment::rowData(spe))) 
    genes <- as.data.frame(SummarizedExperiment::rowData(spe))
  dge <- edgeR::DGEList(counts = counts, samples = samples, genes = genes)
  
  if(!is.null(norm)){
    if (!(norm %in% c("TMM", "upperquartile", "sizefactor"))) {
      stop("Please make sure norm matched one of the following strings: 
         TMM, 
         upperquartile, 
         sizefactor")
    } else {
      dge$samples$norm.factors <- S4Vectors::metadata(spe)$norm.factor
    }
  }
  
  return(dge)
}



utils::globalVariables(c("."))
