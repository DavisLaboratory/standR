#' Add QC statistics to the Spatial Experiment object
#'
#' @param spe_object A SpatialExperiment object
#' @param sample_fraction Double. Genes with low count in more than this threshold of the samples will be removed. Default is 0.9
#' @param rm_genes Logical. Decide whether genes with low count in more than sample_fraction of the samples are removed from the dataset. Default is TRUE.
#' @param min_count Integer. Minimum read count to calculate count threshold. Default is 5.
#' @param design Generate using `model.matrix`, if this is specify, `edgeR::filterByExpr` will be used to filter genes.
#'
#' @return A SpatialExperiment object
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
  stopifnot(length(assayNames(spe)) > 1)
  stopifnot(sample_fraction <= 1 & sample_fraction >= 0)
  stopifnot(min_count >= 0)

  ## add library size to metadata
  colData(spe)$lib_size <- assay(spe, "counts") |> colSums()


  ## calculate lcpm threshold for a gene to be expressed
  L <- mean(colData(spe)$lib_size) * 1e-6
  M <- stats::median(colData(spe)$lib_size) * 1e-6
  lcpm_threshold <- log2(min_count / M + 2 / L)
  S4Vectors::metadata(spe)$lcpm_threshold <- lcpm_threshold

  ## Feature-wise QC & Filter
  ###### Genes with low count in more than the threshold of the samples will be removed
  if (is.null(design)) {
    genes_lowCount_overNsamples <- rowSums(assay(spe, "logcounts") <= lcpm_threshold) >= round(ncol(assay(spe, "logcounts")) * sample_fraction)
    rowData(spe)$genes_lowCount_overNsamples <- genes_lowCount_overNsamples
  } else {
    rowData(spe)$genes_lowCount_overNsamples <- !edgeR::filterByExpr(edgeR::SE2DGEList(spe), design)
  }


  # remove genes and store the removed genes to metadata
  if (rm_genes == TRUE) {
    S4Vectors::metadata(spe)$genes_rm_rawCount <- assay(spe, "counts")[rowData(spe)$genes_lowCount_overNsamples, ]
    S4Vectors::metadata(spe)$genes_rm_logCPM <- assay(spe, "logcounts")[rowData(spe)$genes_lowCount_overNsamples, ]
    spe <- spe[!rowData(spe)$genes_lowCount_overNsamples, ]
  }

  ## Sample-wise QC & Filter
  colData(spe)$countOfLowEprGene <- colSums(assay(spe, "logcounts") <= lcpm_threshold)
  colData(spe)$percentOfLowEprGene <- colSums(assay(spe, "logcounts") <= lcpm_threshold) / nrow(assay(spe, "logcounts")) * 100

  return(spe)
}

utils::globalVariables(c("."))
