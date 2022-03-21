# DEseq2 normalization
calNormCount <- function(spe_object, log = TRUE) {
  # get raw count
  count_df <- SummarizedExperiment::assay(spe_object, 1)
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
#' @param spe_object A spatial experiment object.
#' @param method Normalization method to use. Options: TMM, RPKM, TPM, CPM, upperquartile, sizefactor. RPKM and TPM require gene length information, which should be added into rowData(spe). Note that TMM here is TMM + CPM.
#' @param log Log-transformed or not.
#'
#' @return A spatial experiment object, with the second assay being the normalized count matrix.
#' @export
#'
#' @references Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139-140.
#' @references Love, M., Anders, S., & Huber, W. (2014). Differential analysis of count data–the DESeq2 package. Genome Biol, 15(550), 10-1186.
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
geomxNorm <- function(spe_object, method = c(
                        "TMM", "RPKM", "TPM", "CPM",
                        "upperquartile", "sizefactor"
                      ),
                      log = TRUE) {
  method <- match.arg(method)

  if (!(method %in% c("TMM", "RPKM", "TPM", "CPM", "upperquartile", "sizefactor"))) {
    stop("Please make sure method mathced one of the following strings: TMM,RPKM,TPM,CPM,upperquartile,sizefactor")
  }

  # spe object to dgelist
  spe <- spe_object

  y <- edgeR::SE2DGEList(spe)

  ## TMM

  if (method == "TMM") {
    S4Vectors::metadata(spe)$norm.factor <- edgeR::calcNormFactors(y)$samples$norm.factors
    if (isTRUE(log)) {
      SummarizedExperiment::assay(spe, 2) <- edgeR::calcNormFactors(y) %>%
        edgeR::cpm(., log = TRUE)
    } else {
      SummarizedExperiment::assay(spe, 2) <- edgeR::calcNormFactors(y) %>%
        edgeR::cpm(.)
    }
  }

  ## RPKM

  if (method == "RPKM") {
    if (isTRUE(log)) {
      SummarizedExperiment::assay(spe, 2) <- edgeR::rpkm(y, log = TRUE, prior.count = 0)
    } else {
      SummarizedExperiment::assay(spe, 2) <- edgeR::rpkm(y, log = FALSE, prior.count = 0)
    }
  }


  ## TPM

  if (method == "TPM") {
    if (isTRUE(log)) {
      SummarizedExperiment::assay(spe, 2) <- log(rpkm2tpm(edgeR::rpkm(y)) + 1, 2)
    } else {
      SummarizedExperiment::assay(spe, 2) <- rpkm2tpm(edgeR::rpkm(y))
    }
  }


  ## upper quantile

  if (method == "upperquartile") {
    S4Vectors::metadata(spe)$norm.factor <- edgeR::calcNormFactors(y, method = "upperquartile")$samples$norm.factors
    if (isTRUE(log)) {
      SummarizedExperiment::assay(spe, 2) <- edgeR::calcNormFactors(y, method = "upperquartile") %>%
        edgeR::cpm(., log = TRUE)
    } else {
      SummarizedExperiment::assay(spe, 2) <- edgeR::calcNormFactors(y, method = "upperquartile") %>%
        edgeR::cpm(.)
    }
  }


  ## calculating size factor based on geomean

  if (method == "sizefactor") {
    S4Vectors::metadata(spe)$norm.factor <- calNormCount(spe, log = TRUE)[[2]]
    if (isTRUE(log)) {
      SummarizedExperiment::assay(spe, 2) <- calNormCount(spe, log = TRUE)[[1]]
    } else {
      SummarizedExperiment::assay(spe, 2) <- calNormCount(spe, log = FALSE)[[1]]
    }
  }

  return(spe)
}


utils::globalVariables(c("."))
