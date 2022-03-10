#' Add QC statistics to the Spatial Experiment object
#'
#' @param spe_object A spatial experiment object
#' @param sample_fraction Double. Genes with low count in more than this threshold of the samples will be removed. Default is 0.9
#' @param rm_genes Logical. Decide whether genes with low count in more than sample_fraction of the samples are removed from the dataset. Default is TRUE.
#' @param min_count Interger. Minimum read count to calculate count threshold. Default is 5.
#' @param design Generate using `model.matrix`, if this is specify, `edgeR::filterByExpr` will be used to filter genes.
#'
#' @return A spatial experiment object
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#' spe_filtered <- addPerROIQC(dkd_spe_subset)
#' spe_filtered
#'
addPerROIQC <- function(spe_object, sample_fraction = 0.9, rm_genes = TRUE, min_count = 5, design = NULL) {
  spe <- spe_object

  stopifnot(nrow(spe) > 0)
  stopifnot(ncol(spe) > 0)
  stopifnot(length(SummarizedExperiment::assayNames(spe)) > 1)
  stopifnot(sample_fraction <= 1 & sample_fraction >= 0)
  stopifnot(min_count >= 0)

  ## add library size to metadata
  SummarizedExperiment::colData(spe)$lib_size <- SummarizedExperiment::assay(spe, 1) %>% colSums()


  ## calculate lcpm threshold for a gene to be expressed
  L <- mean(SummarizedExperiment::colData(spe)$lib_size) * 1e-6
  M <- stats::median(SummarizedExperiment::colData(spe)$lib_size) * 1e-6
  lcpm_threshold <- log2(min_count / M + 2 / L)
  S4Vectors::metadata(spe)$lcpm_threshold <- lcpm_threshold

  ## Feature-wise QC & Filter
  ###### Genes with low count in more than the threshold of the samples will be removed
  if (is.null(design)) {
    genes_lowCount_overNsamples <- rowSums(SummarizedExperiment::assay(spe, 2) <= lcpm_threshold) >= round(ncol(SummarizedExperiment::assay(spe, 2)) * sample_fraction)
    SummarizedExperiment::rowData(spe)$genes_lowCount_overNsamples <- genes_lowCount_overNsamples
  } else {
    SummarizedExperiment::rowData(spe)$genes_lowCount_overNsamples <- !edgeR::filterByExpr(edgeR::SE2DGEList(spe), design)
  }


  # remove genes and store the removed genes to metadata
  if (rm_genes == TRUE) {
    S4Vectors::metadata(spe)$genes_rm_rawCount <- SummarizedExperiment::assay(spe, 1)[SummarizedExperiment::rowData(spe)$genes_lowCount_overNsamples, ]
    S4Vectors::metadata(spe)$genes_rm_logCPM <- SummarizedExperiment::assay(spe, 2)[SummarizedExperiment::rowData(spe)$genes_lowCount_overNsamples, ]
    spe <- spe[!SummarizedExperiment::rowData(spe)$genes_lowCount_overNsamples, ]
  }

  ## Sample-wise QC & Filter
  SummarizedExperiment::colData(spe)$countOfLowEprGene <- colSums(SummarizedExperiment::assay(spe, 2) <= lcpm_threshold)
  SummarizedExperiment::colData(spe)$percentOfLowEprGene <- colSums(SummarizedExperiment::assay(spe, 2) <= lcpm_threshold) / nrow(SummarizedExperiment::assay(spe, 2)) * 100

  return(spe)
}

utils::globalVariables(c("."))
