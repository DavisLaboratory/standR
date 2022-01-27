

#' Plot the removed genes with their expressions in all samples from a spatial experiment object after QC
#'
#' @param spe A spatial experiment object.
#' @param top_n Integer. Indicating the top n genes will be plotted.
#' @param point_size Numeric. Point size.
#' @param line_type Character. Line types for ggplot.
#' @param line_col Color for line.
#' @param line_cex Cex for line.
#' @param text_size Text size.
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#'
#' @return A ggplot object.
#' @export
#'
plotRmGenes <- function(spe, top_n, point_size, line_type, line_col, line_cex, text_size){

  . = sample = lcpm = rowname = NULL

  data <- spe@metadata$genes_rm_logCPM %>%
    as.data.frame() %>%
    mutate(m = rowMeans(.)) %>%
    arrange(-.$m)

  top_n <- min(nrow(data), top_n)

  p1 <- data[1:top_n,] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tidyr::gather(sample, lcpm, -rowname) %>%
    ggplot(aes(.$sample, .$lcpm)) +
    geom_point(size = point_size) +
    facet_wrap(~rowname) +
    geom_hline(yintercept = spe@metadata$lcpm_threshold, linetype = line_type,
               cex = line_cex, col = line_col, ) +
    theme_test() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          text = element_text(size = text_size)) +
    xlab(paste0("Samples (n=", ncol(spe),")")) +
    ylab("Expression (logCPM)") +
    ggtitle(paste0("Removed genes (top",as.integer(top_n),")"))

  return(p1)
}



#' Plot the distribution of the percentage of low expressed genes in samples.
#'
#' @param spe A spatial experiment object.
#' @param hist_col Color for histogram.
#' @param hist_fill Fill for histogram.
#' @param bin_num Bin numbers for histogram.
#' @param text_size Text size.
#'
#' @return A ggplot object.
#' @export
#'
plotNEGpercentHist <- function(spe, hist_col, hist_fill, bin_num, text_size) {
  . = NULL
  p2 <- SummarizedExperiment::colData(spe)$percentOfLowEprGene %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    magrittr::set_colnames(c("sample","percent")) %>%
    ggplot(aes(.$percent)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    theme_test() +
    xlab("Percentage of non-expressed genes in each sample (%)") +
    ylab("Frequency") +
    theme(text = element_text(size = text_size)) +
    ggtitle("Distribution")

  return(p2)
}


#' Plot gene-wise QC plot
#'
#' @param spe A spatial experiment object.
#' @param top_n Integer. Indicating the top n genes will be plotted. Default is 9.
#' @param point_size Numeric. Point size.
#' @param line_type Character. Line types for ggplot.
#' @param line_col Color for line.
#' @param line_cex Cex for line.
#' @param hist_col Color for histogram.
#' @param hist_fill Fill for histogram.
#' @param bin_num Bin numbers for histogram.
#' @param text_size Text size.
#' @param layout_ncol Integer. Column number for layout. Default is 1.
#' @param layout_nrow Integer. Row number for layout. Default is 2.
#' @param leyout_height Vector of numerics with length of 2. Default is c(1, .4).
#'
#' @return A ggplot object
#' @export
#'
plotGenewiseQC <- function(spe, top_n = 9, point_size = 1,
                           line_type = "dashed", line_col = "darkred", line_cex = 1,
                           hist_col = "black", hist_fill = "skyblue", bin_num = 30,
                           text_size = 13, layout_ncol = 1, layout_nrow = 2,
                           leyout_height = c(1, .4)){
  p1 <- plotRmGenes(spe, top_n, point_size, line_type, line_col, line_cex, text_size)
  p2 <- plotNEGpercentHist(spe, hist_col, hist_fill, bin_num, text_size)

  print(p1 + p2 + patchwork::plot_layout(layout_ncol, layout_nrow,
                                         heights = leyout_height))
}



#' Plot Sample-wise QC plot
#'
#' @param spe_object A spatial experiment object.
#' @param x_axis Numeric feature to plot as x axis.
#' @param y_axis Numeric feature to plot as y axis.
#' @param x_lab Label name for x axis.
#' @param y_lab Label name for y axis.
#' @param x_threshold Threshold to draw.
#' @param y_threshold Threshold to draw.
#' @param regression_col Color for the regression line.
#' @param hist_col Color for the histograms.
#' @param hist_fill Fill for the histograms.
#' @param bin_num Bin numbers for the histograms.
#' @param threshold_col Threshold line color.
#' @param threshold_linetype Threshold line type.
#' @param layout_ncol Column number layout.
#' @param layout_nrow Row number layout.
#' @param leyout_height Height layout.
#' @param layout_width Width layout.
#'
#'
#'
#' @return A ggplot object.
#' @export
#'
plotLibrarySizeQC <- function(spe_object,
                              x_axis = "AOINucleiCount",
                              y_axis = "lib_size",
                              x_lab = "AOINucleiCount",
                              y_lab = "Library size",
                              x_threshold = NULL,
                              y_threshold = NULL,
                              regression_col = "purple",
                              hist_col = "black", hist_fill = "white", bin_num = 50,
                              threshold_col = "red", threshold_linetype = "dashed",
                              layout_ncol = 2, layout_nrow = 2, leyout_height = c(0.8, 2.5), layout_width = c(2.5, 0.8)){

  . = NULL

  stopifnot(x_axis %in% colnames(SummarizedExperiment::colData(spe_object)))
  stopifnot(y_axis %in% colnames(SummarizedExperiment::colData(spe_object)))

  p1 <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame() %>%
    ggplot(aes_string(x_axis, y_axis)) +
    geom_point() +
    geom_smooth(method='loess', se = FALSE, col = regression_col) +
    theme_test() +
    xlab(x_lab) +
    ylab(y_lab)

  p2 <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame() %>%
    ggplot(aes_string(y_axis)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    theme_test() +
    coord_flip() +
    ylab("Frequency") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())

  p_blank <- ggplot() + theme_void()

  p3 <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame() %>%
    ggplot(aes_string(x_axis)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    theme_test() +
    ylab("Frequency") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())

  if(!is.null(x_threshold)){
    stopifnot(is.numeric(x_threshold))
    p1 <- p1 + geom_vline(xintercept = x_threshold, col = threshold_col, linetype = threshold_linetype)
    p3 <- p3 + geom_vline(xintercept = x_threshold, col = threshold_col, linetype = threshold_linetype)
  }

  if(!is.null(y_threshold)){
    stopifnot(is.numeric(y_threshold))
    p1 <- p1 + geom_hline(yintercept = y_threshold, col = threshold_col, linetype = threshold_linetype)
    p2 <- p2 + geom_vline(xintercept = y_threshold, col = threshold_col, linetype = threshold_linetype)
  }

  print(p3 + p_blank + p1 + p2 + patchwork::plot_layout(layout_ncol, layout_nrow,
                                                        widths = layout_width,
                                                        heights = leyout_height))
}
