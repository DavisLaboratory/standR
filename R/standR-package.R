#' @import ggplot2
#' @importFrom methods is
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @import patchwork
#' @importFrom SingleCellExperiment reducedDim
NULL


#' Tools for nalysing Nanostring's GeoMX spatial transcriptomics data
#'
#' `standR` implements a series of functions to facilitate inspection,
#' analysis and visualization of the Nanostring's GeoMX DSP datasets. `standR`
#' takes the either the csv files from the Nanostring or DGEList object as
#' input, allowing for multiple methods to be analyzed together.
#'
#' `standR` represents the GeoMX DSP data as SpatialExperiment objects, which can
#' easily be integrated with a wide variety of Bioconductor packages.`standR`
#' generates various plots, such as QC distribution plots, dimension reduction
#' plots and RLE plots, for quality control of genes and region of interest
#' (ROI) samples. Multiple normalisation and batch correction methods are also
#' provided in the package as well, with the ability to compute statistics for
#' assessing the normalisation/batch correction results.
#'
#'
#' @author Ning Liu \email{liu.n@@wehi.edu.au}
#' @name standR-package
#' @docType package
#' @aliases standR standR-package
#' @keywords internal
#'
NULL

#' Description of the standR example datasets
#'
#' standR-package has 1 datasets: \itemize{
#'   \item dkd_spe_subset  Example subset of a GeoMX DSP WTA dataset,
#'  }
#'
#' @docType data
#' @name dkd_spe_subset
#' @usage data("dkd_spe_subset")
#' @keywords internal
#' @format A spatial experiment object with 3000 rows and 70 samples:
#' @source \url{http://nanostring-public-share.s3-website-us-west-2.amazonaws.com/GeoScriptHub/KidneyDataset/}
#' @examples
#' data(dkd_spe_subset)
"dkd_spe_subset"
