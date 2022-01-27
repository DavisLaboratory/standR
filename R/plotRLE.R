#' Plot Relative log expression
#'
#' @param spe_object A spatial experiment object
#' @param n_assay integer. the nth assay in the spe object to plot.
#' @param colorby feature to color by.
#' @param sortby feature to sort by.
#' @param font_size font size.
#' @param line_cor line color.
#'
#' @return A ggplot object
#' @export
#'
plotRLE <- function(spe_object, n_assay = 1, colorby = NA, sortby = colorby, font_size = 20, line_cor = "red"){

  . = med = sample = count = NULL

  stopifnot(is.numeric(n_assay))
  stopifnot(n_assay <= length(spe_object@assays))

  pre_rle <- SummarizedExperiment::assay(spe_object, n_assay) %>%
    as_tibble() %>%
    mutate(med = rowMeans(.)) %>%
    mutate_at(vars(-med), funs(. - med)) %>%
    dplyr::select(-med) %>%
    tidyr::gather(sample, count) %>%
    left_join(SummarizedExperiment::colData(spe_object) %>% as.data.frame() %>% rownames_to_column(), by = c("sample"="rowname"))

  p <- graphics::boxplot(count ~ sample, data = pre_rle, outline = F, plot = FALSE)

  boxstat <- p$stats

  if(!is.na(colorby)){
    stopifnot(colorby %in% colnames(pre_rle))
    stopifnot(sortby %in% colnames(pre_rle))

    p1 <- pre_rle %>%
      arrange((!!sym(sortby))) %>%
      mutate(sample = factor(sample, levels = unique(.$sample))) %>%
      ggplot(aes(sample, count)) +
      geom_boxplot(outlier.shape = NA, aes_string(fill = colorby), col = "black") +
      #geom_jitter(aes_string(col = colorby)) +
      #stat_boxplot(geom = "errorbar", linetype = 2, width = .5) +
      #stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
      ylim(c(min(boxstat[1,])-1, max(boxstat[5,])+1)) +
      xlab('Samples') +
      ylab("Relative log expression") +
      theme_test() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size = font_size)) +
      ggtitle(names(spe_object@assays@data)[n_assay]) +
      scale_fill_brewer(palette="Dark2") +
      geom_hline(yintercept = 0, col = line_cor, linetype = "dashed", cex = 1)
  } else {
    p1 <- pre_rle %>%
      ggplot(aes(sample, count)) +
      geom_boxplot(outlier.shape = NA, col = "black") +
      #geom_jitter(aes_string(col = colorby)) +
      #stat_boxplot(geom = "errorbar", linetype = 2, width = .5) +
      #stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
      ylim(c(min(boxstat[1,]), max(boxstat[5,]))) +
      xlab('Samples') +
      ylab("Relative log expression") +
      theme_test() +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size = font_size)) +
      ggtitle(names(spe_object@assays@data)[n_assay]) +
      scale_fill_brewer(palette="Dark2") +
      geom_hline(yintercept = 0, col = line_cor, linetype = "dashed", cex = 1)
  }

  return(p1)
}
