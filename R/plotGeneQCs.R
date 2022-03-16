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
#' @param layout_height Vector of numerics with length of 2. Default is c(1, .4).
#' @param ... aesthetic mappings to pass to `ggplot2::aes()` of the dot plots.
#' @param ordannots variables or computations to sort samples by (tidy style).
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#' spe <- addPerROIQC(dkd_spe_subset)
#' plotGeneQC(spe)
plotGeneQC <- function(spe, top_n = 9, ordannots = c(), point_size = 1,
                       line_type = "dashed", line_col = "darkred", line_cex = 1,
                       hist_col = "black", hist_fill = "skyblue", bin_num = 30,
                       text_size = 13, layout_ncol = 1, layout_nrow = 2,
                       layout_height = c(1, .4), ...) {
  p1 <- suppressWarnings(plotRmGenes(spe, top_n, ordannots, point_size, line_type, line_col, line_cex, text_size, ...))
  p2 <- suppressWarnings(plotNEGpercentHist(spe, hist_col, hist_fill, bin_num, text_size))

  suppressWarnings(print(p1 + p2 + patchwork::plot_layout(layout_ncol, layout_nrow,
    heights = layout_height
  )))
}


# plot removed genes
plotRmGenes <- function(spe, top_n, ordannots, point_size, line_type,
                        line_col, line_cex, text_size, ...) {
  # get order

  sampleorder <- SummarizedExperiment::colData(spe) %>%
    as.data.frame(optional = TRUE) %>%
    orderSamples(., ordannots)

  # get removed genes and order by mean expression
  data <- S4Vectors::metadata(spe)$genes_rm_logCPM %>%
    as.data.frame() %>%
    mutate(m = rowMeans(.)) %>%
    arrange(m) %>%
    .[, sampleorder]

  # top N genes to plot
  top_n <- min(nrow(data), top_n)

  # plotting

  aesmap <- rlang::enquos(...)

  data[1:top_n, ] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tidyr::gather(sample, lcpm, -rowname) %>%
    left_join(as.data.frame(SummarizedExperiment::colData(spe)) %>%
      rownames_to_column(), by = c("sample" = "rowname")) %>%
    mutate(sample = factor(sample, levels = colnames(data))) %>%
    ggplot(aes(sample, lcpm, !!!aesmap)) +
    geom_point(size = point_size, alpha = .5) +
    scale_colour_discrete(na.translate = F) +
    scale_shape_discrete(na.translate = F) +
    facet_wrap(~rowname) +
    geom_hline(
      yintercept = S4Vectors::metadata(spe)$lcpm_threshold, linetype = line_type,
      cex = line_cex, col = line_col,
    ) +
    theme_test() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size = text_size)
    ) +
    xlab(paste0("Samples (n=", ncol(spe), ")")) +
    ylab("Expression (logCPM)") +
    ggtitle(paste0("Removed genes (top", as.integer(top_n), ")"))

}



# plot non-expressed gene percentage histogram
plotNEGpercentHist <- function(spe, hist_col, hist_fill, bin_num, text_size) {
  SummarizedExperiment::colData(spe)$percentOfLowEprGene %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    magrittr::set_colnames(c("sample", "percent")) %>%
    ggplot(aes(percent)) +
    geom_histogram(col = hist_col, fill = hist_fill, bins = bin_num) +
    theme_test() +
    xlab("Percentage of non-expressed genes in each sample (%)") +
    ylab("Frequency") +
    theme(text = element_text(size = text_size)) +
    ggtitle("Distribution")

}

utils::globalVariables(c(".", "sample", "lcpm", "rowname", "m", "percent"))
