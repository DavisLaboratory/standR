#' @import ggplot2
#' @import patchwork
NULL

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
#' @param ... aesthetic mappings to pass to `ggplot2::aes()` of the dot plots.
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
                              layout_ncol = 2, layout_nrow = 2,
                              leyout_height = c(0.8, 2.5), layout_width = c(2.5, 0.8),
                              ...){

  . = NULL

  stopifnot(x_axis %in% colnames(SummarizedExperiment::colData(spe_object)))
  stopifnot(y_axis %in% colnames(SummarizedExperiment::colData(spe_object)))

  aesmap = rlang::enquos(...)
  x_axis = rlang::sym(x_axis)
  y_axis = rlang::sym(y_axis)

  # plot dot plot
  p1 <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame() %>%
    ggplot(aes(!!x_axis, !!y_axis, !!!aesmap)) +
    geom_point() +
    geom_smooth(method='loess', se = FALSE, col = regression_col) +
    theme_test() +
    xlab(x_lab) +
    ylab(y_lab)

  # plot distribution of y axis
  p2 <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame() %>%
    ggplot(aes(!!y_axis)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    theme_test() +
    coord_flip() +
    ylab("Frequency") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())

  p_blank <- ggplot() + theme_void()

  # plot distribution of x axis
  p3 <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame() %>%
    ggplot(aes(!!x_axis)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    theme_test() +
    ylab("Frequency") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())

  # plot threshold
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
                                                        heights = leyout_height, guides = "collect"))
}

