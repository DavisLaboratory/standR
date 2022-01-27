
#' Plot PCA with two dimensions
#'
#' @param spe_object A spatial experiment object.
#' @param featureToPlot Feature to color by.
#' @param featureToShape Feature to shape by.
#' @param mainPCs Which PCs to plot.
#' @param color Color values.
#' @param legend.title Legend title.
#' @param font.size Font size.
#' @param shape.size Point size.
#' @param plotName if TRUE, plot the name of the samples.
#'
#' @importFrom SingleCellExperiment reducedDim
#'
#' @return A ggplot object.
#' @export
#'
plot2dimensionPCA <- function(spe_object, featureToPlot, featureToShape = NA,
                              mainPCs = c("PC1", "PC2"), color = NA, legend.title = NA,
                              font.size = 15, shape.size = 6, plotName = FALSE){

  . = NULL
  pca_object <- reducedDim(spe_object, "PCA")
  var_exp <- attr(pca_object, "percentVar") %>%
    round(., 2)
  names(var_exp) <- colnames(pca_object)

  if(length(mainPCs)!=2){
    stop("'length(mainPCs)' is not equal to 2.")
  }

  if(!startsWith(mainPCs[1],"PC") | !startsWith(mainPCs[2],"PC")){
    stop("Each element in mainPCs should be start with PC.")
  }

  pca_plots <- pca_object %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    left_join(SummarizedExperiment::colData(spe_object) %>%
                as.data.frame() %>%
                rownames_to_column())
  if(!is.na(featureToShape)){
    stopifnot(featureToPlot %in% colnames(pca_plots))
    stopifnot(featureToShape %in% colnames(pca_plots))
    pca_plots <- pca_plots %>%
      ggplot(aes_string(mainPCs[1], mainPCs[2], col = featureToPlot, shape = featureToShape))
  } else{
    stopifnot(featureToPlot %in% colnames(pca_plots))
    pca_plots <- pca_plots %>%
      ggplot(aes_string(mainPCs[1], mainPCs[2], col = featureToPlot))
  }
  pca_plots <- pca_plots +
    geom_point(size = shape.size, alpha = .7) +
    theme_test() +
    xlab(paste0(mainPCs[1],"(", var_exp[mainPCs[1]],"%)")) +
    ylab(paste0(mainPCs[2],"(", var_exp[mainPCs[2]],"%)")) +
    theme(text = element_text(size = font.size))

  if(isTRUE(plotName)){
    pca_plots <- pca_plots +
      geom_text(label = rownames(pca_object))
  }

  if(!is.na(color[1])){
    if(length(color)!=length(unique(SummarizedExperiment::colData(spe_object)[,featureToPlot]))){
      stop("length of color should be the same as feature")
    }
    pca_plots <- pca_plots +
      scale_color_manual(values=color)
  }

  if(!all(is.na(legend.title))){
    if(length(legend.title)==1){
      pca_plots <- pca_plots +
        guides(col=guide_legend(title=legend.title))
    } else if(length(legend.title)==2){
      pca_plots <- pca_plots +
        guides(col=guide_legend(title=legend.title[1]),
               shape=guide_legend(title=legend.title[2]))
    }
  }

  return(pca_plots)
}



#' Set up the grid for multiple-panel PCA plots
#'
#' @param x x
#' @param y y
#' @param include.equals Logical
#'
#' @export
#'
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
#' @param featureToPlot Feature to color by.
#' @param featureToShape Feature to shape by.
#' @param n_dimension The top n dimensions to be plotted
#' @param color Color values.
#' @param legend.title Legend title.
#' @param font.size Font size.
#' @param shape.size Dot size.
#'
#' @return A ggplot object.
#' @export
#'
plotPairdimensionPCA <- function(spe_object, featureToPlot, featureToShape = NA,
                                 n_dimension = 3, color = NA, legend.title = NA,
                                 font.size = 15, shape.size = 6){

  . = NULL

  stopifnot(is.numeric(n_dimension))

  pca_object <- reducedDim(spe_object, "PCA")

  n = n_dimension

  PCnames <- colnames(pca_object)[1:n]

  PCsToPlot <- expand.grid.rmdup(PCnames, PCnames) # get pairs of PCs

  m <- matrix(seq((n-1)^2),n-1,n-1, byrow = TRUE) # make position matrix

  index_realPlots <- m[upper.tri(m,diag = TRUE)] %>% sort() # these are the positions of real plots

  index_emptyPlots <- seq((n-1)^2) %>% setdiff(index_realPlots) # these are blank spaces

  plot_blank <- ggplot() + theme_void()

  realplots <- list()
  for(i in seq(nrow(PCsToPlot))){
    realplots[[i]] <- plot2dimensionPCA(spe_object, featureToPlot,
                                        featureToShape = featureToShape,
                                        mainPCs = PCsToPlot[i,], color = color,
                                        legend.title = legend.title,
                                        font.size = font.size,
                                        shape.size = shape.size)
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

  PCApairplot <- ggpubr::ggarrange(plotlist = plotting_list,
                                   ncol = n-1,nrow = n-1,
                                   common.legend = TRUE, legend = "top")

  return(PCApairplot)
}




#' Plot the PCA scree plot.
#'
#' @param spe_object A spatial experiment object.
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
plot_screePCA <- function(spe_object, bar_color = "black", bar_fill = "royalblue", bar_width = .8, point_col = "tomato3", line_col = "tomato3", point_size = 2){

  . = NULL
  pca_object <- reducedDim(spe_object, "PCA")
  var_exp <- attr(pca_object, "percentVar") %>%
    round(., 2)
  names(var_exp) <- colnames(pca_object)

  scree_plot <- var_exp %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("ev")) %>%
    rownames_to_column() %>%
    mutate(rowname = factor(.$rowname, levels = .$rowname)) %>%
    mutate(csum = cumsum(.$ev)) %>%
    ggplot(aes(.$rowname, .$ev)) +
    geom_bar(stat = "identity", col = bar_color, fill = bar_fill, width = bar_width) +
    geom_point(aes(x = .$rowname, y = .$csum), size = point_size, col = point_col) +
    geom_line(aes(x = .$rowname, y = .$csum), group = 1, col = line_col) +
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
#' @param featureToPlot Feautre to color by.
#' @param n_loadings Plot the top n gene loadings
#' @param mainPCs PCs to plot.
#' @param featureToShape Feature to shape by.
#' @param color Color values
#' @param legend.title Legend title.
#' @param font.size Font size.
#' @param shape.size Point size.
#'
#' @return A ggplot object.
#' @export
#'
plotPCAbiplot <- function(spe_object, featureToPlot, featureToShape = NA, n_loadings = 10, mainPCs = c("PC1", "PC2"), color = NA, legend.title = NA,
                          font.size = 15, shape.size = 6){

  . = NULL
  pca_object <- reducedDim(spe_object, "PCA")
  loadings <- attr(pca_object, "rotation")

  pc1_genes <- loadings %>%
    as.data.frame() %>%
    arrange(-abs(!!sym(mainPCs[1]))) %>%
    .[1:n_loadings,] %>%
    rownames()

  pc2_genes <- loadings %>%
    as.data.frame() %>%
    arrange(-abs(!!sym(mainPCs[2]))) %>%
    .[1:n_loadings,] %>%
    rownames()

  genes2plot <- unique(pc1_genes, pc2_genes)

  r <- min(
    (max(pca_object[,mainPCs[1]]) - min(pca_object[,mainPCs[1]]) /
       (max(loadings[,mainPCs[1]]) - min(loadings[,mainPCs[1]]))),
    (max(pca_object[,mainPCs[2]]) - min(pca_object[,mainPCs[2]]) /
       (max(loadings[,mainPCs[2]]) - min(loadings[,mainPCs[2]]))))

  loadings2plot <- loadings %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    dplyr::select(c(.$rowname,mainPCs[1],mainPCs[2])) %>%
    filter(.$rowname %in% genes2plot) %>%
    mutate(x_end = caculate_coor_end((!!sym(mainPCs[1])), r),
           y_end = caculate_coor_end((!!sym(mainPCs[2])), r))


  plot2dimensionPCA(spe_object, featureToPlot, featureToShape,
                    mainPCs = mainPCs, color = color,
                    legend.title = legend.title,
                    font.size = font.size,
                    shape.size = shape.size) +
    geom_segment(data = loadings2plot,
                 aes(x = 0, y = 0,
                     xend = .$x_end,
                     yend = .$y_end),
                 arrow = arrow(length = unit(1/2, 'picas'), ends = 'last'),
                 color = "black",
                 size = .5,
                 alpha = 1,
                 show.legend = NA) +
    ggrepel::geom_text_repel(data = loadings2plot,
                             aes(label = loadings2plot$rowname,
                                 x = loadings2plot$x_end,
                                 y = loadings2plot$y_end,
                                 hjust = 0),
                             color = "black",
                             size = 3)
}
