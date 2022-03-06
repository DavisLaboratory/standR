# expand matrix from vector
expand.grid.rmdup <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)

  y <- unique(y)

  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])

    if(length(z)) cbind(x[i], z, deparse.level=0)
  }

  do.call(rbind, lapply(seq_along(x), g))
}

#' Plot pair-wise PCA plots for multiple dimensions
#'
#' @param spe_object A spatial experiment object.
#' @param n_dimension The top n dimensions to be plotted
#' @param precomputed a dimensional reduction results from `stats::prcomp`.
#'   result in `reducedDims(object)` to plot. Default is NULL,
#'   we will compute for you.
#' @param assay a numeric or character, specifying the assay to use (for
#'   `SummarizedExperiment` and its derivative classes).
#' @param ... aesthetic mappings to pass to `ggplot2::aes()`.
#' @param title Character vector, title to put at the top.
#' @param title.size Numeric vector, size of the title.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#' plotPairdimensionPCA(dkd_spe_subset)
plotPairdimensionPCA <- function(spe_object, n_dimension = 3,
                                 precomputed = NULL,
                                 assay = 1, title = NA, title.size = 14, ...){

  set.seed(44)

  . = NULL

  stopifnot(is.numeric(n_dimension))

  #compute PCA
  if (is.null(precomputed)) {
    pcdata = calcPCA(SummarizedExperiment::assay(spe_object, assay), n_dimension)
  } else {
    pcdata = checkPrecomputedPCA(spe_object, precomputed)
  }

  pca_object <- pcdata

  n = n_dimension

  PCsToPlot <- expand.grid.rmdup(1:n, 1:n) # get pairs of PCs

  m <- matrix(seq((n-1)^2),n-1,n-1, byrow = TRUE) # make position matrix

  index_realPlots <- m[upper.tri(m,diag = TRUE)] %>% sort() # these are the positions of real plots

  index_emptyPlots <- seq((n-1)^2) %>% setdiff(index_realPlots) # these are blank spaces

  plot_blank <- ggplot() + theme_void()

  realplots <- list()
  for(i in seq(nrow(PCsToPlot))){
    realplots[[i]] <- plotPCA(spe_object,
                              dims = PCsToPlot[i,],
                              precomputed = precomputed,
                              assay = assay, ...)
  }

  plotting_list <- list()
  k = 1
  for(j in seq((n-1)^2)){
    if(j %in% index_realPlots){
      plotting_list[[j]] <- realplots[[k]]
      k = k+1
    } else if(j %in% index_emptyPlots){
      plotting_list[[j]] <- plot_blank
    }
  }

  checkPackages("ggpubr")

  PCApairplot <- ggpubr::ggarrange(plotlist = plotting_list,
                                   ncol = n-1,nrow = n-1,
                                   common.legend = TRUE, legend = "top")

  if(!is.na(title)){
    PCApairplot <- ggpubr::annotate_figure(PCApairplot,
                                   top = ggpubr::text_grob(title,
                                                   color = "black",
                                                   face = "bold",
                                                   size = title.size))
  }


  return(PCApairplot)
}




#' Plot the PCA scree plot.
#'
#' @param spe_object A spatial experiment object.
#' @param dims The top n dimensions to be plotted
#' @param precomputed a dimensional reduction results from `stats::prcomp`.
#'   result in `reducedDims(object)` to plot. Default is NULL,
#'   we will compute for you.
#' @param assay a numeric or character, specifying the assay to use (for
#'   `SummarizedExperiment` and its derivative classes).
#' @param bar_color Color for bar.
#' @param bar_fill Fill for bar.
#' @param bar_width Bar width.
#' @param point_col Color for point.
#' @param line_col Color for line.
#' @param point_size Point size.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#' plotScreePCA(dkd_spe_subset, dims = 10)
plotScreePCA <- function(spe_object, dims = ncol(spe_object), precomputed = NULL, assay = 1,
                         bar_color = "black",bar_fill = "royalblue", bar_width = .8,
                         point_col = "tomato3", line_col = "tomato3",
                         point_size = 2){

  . = rowname = ev = csum = NULL

  #compute PCA
  if (is.null(precomputed)) {
    pcdata = calcPCA(SummarizedExperiment::assay(spe_object, assay), dims)
  } else {
    pcdata = checkPrecomputedPCA(spe_object, precomputed)
  }

  pca_object <- pcdata

  var_exp <- attr(pca_object, "percentVar") %>%
    round(., 2)
  names(var_exp) <- colnames(pca_object)

  scree_plot <- var_exp %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("ev")) %>%
    rownames_to_column() %>%
    mutate(rowname = factor(.$rowname, levels = .$rowname)) %>%
    mutate(csum = cumsum(.$ev)) %>%
    ggplot(aes(rowname, ev)) +
    geom_bar(stat = "identity", col = bar_color, fill = bar_fill, width = bar_width) +
    geom_point(aes(x = rowname, y = csum), size = point_size, col = point_col) +
    geom_line(aes(x = rowname, y = csum), group = 1, col = line_col) +
    theme_test() +
    xlab("PCs") +
    ylab("Explained variance(%)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylim(0,100)

  return(scree_plot)
}



caculate_coor_end <- function(x, r, f = 1.5){
  xend <- x * r * f
  return(xend)
}


#' Plot PCA bi plot
#'
#' @param spe_object A spatial experiment object.
#' @param n_loadings Plot the top n gene loadings
#' @param dims The top n dimensions to be plotted
#' @param precomputed a dimensional reduction results from `stats::prcomp`.
#'   result in `reducedDims(object)` to plot. Default is NULL,
#'   we will compute for you.
#' @param assay a numeric or character, specifying the assay to use (for
#'   `SummarizedExperiment` and its derivative classes).
#' @param arrow_x a numeric, indicating the x coordinate of the base of the arrow.
#' @param arrow_y a numeric, indicating the y coordinate of the base of the arrow.
#' @param ... aesthetic mappings to pass to `ggplot2::aes()`.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#' plotPCAbiplot(dkd_spe_subset)
plotPCAbiplot <- function(spe_object, n_loadings = 10,
                          dims = c(1,2), precomputed = NULL, assay = 1,
                          arrow_x = 0, arrow_y = 0,
                          ...){

  . = rowname = x_end = y_end = NULL

  #compute PCA
  if (is.null(precomputed)) {
    pcdata = calcPCA(SummarizedExperiment::assay(spe_object, assay), dims)
  } else {
    pcdata = checkPrecomputedPCA(spe_object, precomputed)
  }

  pca_object <- pcdata

  loadings <- attr(pca_object, "rotation")

  if(is(n_loadings, "numeric")){
    pc1_genes <- loadings %>%
      as.data.frame() %>%
      arrange(-abs(!!sym(paste0("PC",dims[1])))) %>%
      .[1:n_loadings,] %>%
      rownames()
    pc2_genes <- loadings %>%
      as.data.frame() %>%
      arrange(-abs(!!sym(paste0("PC",dims[2])))) %>%
      .[1:n_loadings,] %>%
      rownames()
  } else if(is(n_loadings, "character")){
    stopifnot(n_loadings %in% rownames(loadings))
    pc1_genes <- loadings[n_loadings,] %>%
      rownames()
    pc2_genes <- pc1_genes
  }

  genes2plot <- unique(pc1_genes, pc2_genes)

  r <- min(
    (max(pca_object[,dims[1]]) - min(pca_object[,dims[1]]) /
       (max(loadings[,dims[1]]) - min(loadings[,dims[1]]))),
    (max(pca_object[,dims[2]]) - min(pca_object[,dims[2]]) /
       (max(loadings[,dims[2]]) - min(loadings[,dims[2]]))))

  loadings2plot <- loadings %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    dplyr::select(c(rowname,paste0("PC",dims[1]),paste0("PC",dims[2]))) %>%
    filter(rowname %in% genes2plot) %>%
    mutate(x_end = caculate_coor_end((!!sym(paste0("PC",dims[1]))), r),
           y_end = caculate_coor_end((!!sym(paste0("PC",dims[2]))), r))


  plotPCA(spe_object, dims = dims,
          precomputed = precomputed, assay = assay, ...) +
    geom_segment(data = loadings2plot,
                 aes(x = arrow_x, y = arrow_y,
                     xend = x_end,
                     yend = y_end),
                 arrow = arrow(length = unit(1/2, 'picas'), ends = 'last'),
                 color = "black",
                 size = .5,
                 alpha = 1,
                 show.legend = NA) +
    ggrepel::geom_text_repel(data = loadings2plot,
                             aes(label = rowname,
                                 x = x_end,
                                 y = y_end,
                                 hjust = 0),
                             color = "black",
                             size = 3)
}

