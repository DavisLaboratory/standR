#' @import ggplot2
NULL

#' Plot the user-defined meta data using alluvium plot
#'
#' @param spe_object A spatial experiment object.
#' @param column2plot Which columns to plot.
#' @param textsize text size.
#'
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' library(ggalluvial)
#' data("dkd_spe_subset")
#' plotMetadata(dkd_spe_subset, column2plot = c("SlideName","disease_status","region"))
#'
plotMetadata <- function(spe_object, column2plot, textsize = 3){

  checkPackages("ggalluvial")

  stopifnot(all(column2plot %in% colnames(SummarizedExperiment::colData(spe_object))))
  x = stratum = alluvium = NULL

  #using alluvium plot to visualise user-defined metadata in the data.
  SummarizedExperiment::colData(spe_object) %>%
    as.data.frame() %>%
    dplyr::select(column2plot) %>%
    ggalluvial::to_lodes_form() %>%
    ggplot(aes(x = x, stratum = stratum,
               alluvium = alluvium, fill = stratum, label = stratum)) +
    ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "frontback",
              color = "NA") +
    ggalluvial::geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = textsize) +
    theme_light() +
    theme(legend.position = "none") +
    xlab("Metadata") +
    ylab("Frequency")
}
