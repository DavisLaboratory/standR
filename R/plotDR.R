#' @import ggplot2
NULL

checkPackages <- function(...) {
  pkgs = list(...)
  names(pkgs) = pkgs
  pkgs = sapply(pkgs, requireNamespace, quietly = TRUE)
  if (!all(pkgs)) {
    pkgs = names(pkgs)[!pkgs]
    pkgs = paste(paste0("'", pkgs, "'"), collapse = ', ')
    stop(sprintf(
      'The following packages need to be installed for this feature to work: %s',
      pkgs
    ))
  }
}

bhuvad_theme <- function (rl = 1.1) {
  stopifnot(rl > 0)
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = element_rect(colour = 'black', fill = NA),
      panel.grid = element_blank(),
      axis.title = element_text(size = rel(rl) * 1.1),
      axis.text = element_text(size = rel(rl)),
      plot.title = element_text(size = rel(rl) * 1.2),
      strip.background = element_rect(fill = NA, colour = 'black'),
      strip.text = element_text(size = rel(rl)),
      legend.text = element_text(size = rel(rl)),
      legend.title = element_text(size = rel(rl), face = 'italic')
    )
}

#' Compute and plot the results of a PCA analysis on gene expression data
#'
#' @param edata a DGEList, SummarizedExperiment or ExpressionSet object
#'   containing gene expression data.
#' @param dims a numeric, containing 2 values specifying the dimensions to plot.
#' @param assay a numeric or character, specifying the assay to use (for
#'   `SummarizedExperiment` and its derivative classes).
#' @param precomputed a dimensional reduction results from `stats::prcomp`.
#'   result in `reducedDims(object)` to plot.
#' @param rl a numeric, specifying the relative scale factor to apply to text on
#'   the plot.
#' @param ... aesthetic mappings to pass to `ggplot2::aes_string()`.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' se = emtdata::cursons2018_se()
#' dge = emtdata::asDGEList(se)
#' plotPCA(dge, colour = Subline)
#'
setGeneric("plotPCA",
           function(edata,
                    dims = c(1, 2),
                    ...) standardGeneric("plotPCA"))

#' @rdname plotPCA
setMethod("plotPCA",
          signature('DGEList','ANY'),
          function(edata, dims, precomputed = NULL, rl = 1, ...){
            #compute PCA
            if (is.null(precomputed)) {
              pcdata = calcPCA(edgeR::cpm(edata, log = TRUE), dims)
            } else {
              pcdata = checkPrecomputedPCA(edata, precomputed)
            }

            #extract sample data
            sdata = edata$samples
            #create data structure
            drdf = pdataPC_intl(pcdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotPCA
setMethod("plotPCA",
          signature('ExpressionSet','ANY'),
          function(edata, dims, precomputed = NULL, rl = 1, ...){
            #compute PCA
            if (is.null(precomputed)) {
              pcdata = calcPCA(Biobase::exprs(edata), dims)
            } else {
              pcdata = checkPrecomputedPCA(edata, precomputed)
            }

            #extract sample data
            sdata = Biobase::pData(edata)
            #create data structure
            drdf = pdataPC_intl(pcdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotPCA
setMethod("plotPCA",
          signature('SummarizedExperiment', 'ANY'),
          function(edata, dims, assay = 1, precomputed = NULL, rl = 1, ...){
            #compute PCA
            if (is.null(precomputed)) {
              pcdata = calcPCA(SummarizedExperiment::assay(edata, assay), dims)
            } else {
              pcdata = checkPrecomputedPCA(edata, precomputed)
            }

            #extract sample data
            if (is(edata, 'ExperimentList')) {
              colFun = ExperimentList::colWithExperimentData
            } else {
              colFun = SummarizedExperiment::colData
            }
            sdata = BiocGenerics::as.data.frame(colFun(edata), optional = TRUE)
            #create data structure
            drdf = pdataPC_intl(pcdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })


#' Compute and plot the results of a PCA analysis on gene expression data
#'
#' @param dimred a string or integer scalar indicating the reduced dimension
#'   result in `reducedDims(object)` to plot.
#'
#' @inheritParams plotPCA
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' se = emtdata::cursons2018_se()
#' dge = emtdata::asDGEList(se)
#' plotPCA(dge, colour = Subline)
#'
setGeneric("plotDR",
           function(edata,
                    dims = c(1, 2),
                    ...) standardGeneric("plotDR"))
#' @rdname plotDR
setMethod("plotDR",
          signature('SingleCellExperiment', 'ANY'),
          function(edata, dims, dimred = 'PCA', rl = 1, ...){
            #compute PCA
            pcdata = SingleCellExperiment::reducedDim(edata, type = dimred)

            #extract sample data
            if (is(edata, 'ExperimentList')) {
              colFun = ExperimentList::colWithExperimentData
            } else {
              colFun = SummarizedExperiment::colData
            }
            sdata = BiocGenerics::as.data.frame(colFun(edata), optional = TRUE)
            #create data structure
            drdf = pdataPC_intl(pcdata, dims, relabel = FALSE)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotDR
setMethod("plotDR",
          signature('SpatialExperiment', 'ANY'),
          function(edata, dims, dimred = 'PCA', rl = 1, ...){
            #compute PCA
            pcdata = SingleCellExperiment::reducedDim(edata, type = dimred)

            #extract sample data
            if (is(edata, 'ExperimentList')) {
              colFun = ExperimentList::colWithExperimentData
            } else {
              colFun = SummarizedExperiment::colData
            }
            sdata = BiocGenerics::as.data.frame(colFun(edata), optional = TRUE)
            #create data structure
            drdf = pdataPC_intl(pcdata, dims, relabel = FALSE)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' Compute and plot the results of a MDS analysis on gene expression data
#'
#' @param precomputed a dimensional reduction results from either
#'   `limma::plotMDS`.
#'
#' @inheritParams plotPCA
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' se = emtdata::cursons2018_se()
#' dge = emtdata::asDGEList(se)
#' plotMDS(dge, colour = Subline)
#'
setGeneric("plotMDS",
           function(edata,
                    dims = c(1, 2),
                    precomputed = NULL,
                    rl = 1,
                    ...) standardGeneric("plotMDS"))

#' @rdname plotMDS
setMethod("plotMDS",
          signature('DGEList','ANY', 'ANY', 'ANY'),
          function(edata, dims, precomputed, rl, ...){
            #compute PCA
            if (is.null(precomputed)) {
              mdsdata = limma::plotMDS(edata, plot = FALSE)
            } else {
              mdsdata = checkPrecomputedMDS(edata, precomputed)
            }

            #extract sample data
            sdata = edata$samples
            #create data structure
            drdf = pdataMDS_intl(mdsdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotMDS
setMethod("plotMDS",
          signature('ExpressionSet','ANY', 'ANY', 'ANY'),
          function(edata, dims, precomputed, rl, ...){
            #compute PCA
            if (is.null(precomputed)) {
              mdsdata = limma::plotMDS(edata, plot = FALSE)
            } else {
              mdsdata = checkPrecomputedMDS(edata, precomputed)
            }

            #extract sample data
            sdata = Biobase::pData(edata)
            #create data structure
            drdf = pdataMDS_intl(mdsdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

#' @rdname plotMDS
setMethod("plotMDS",
          signature('SummarizedExperiment','ANY', 'ANY', 'ANY'),
          function(edata, dims, precomputed, rl, ...){
            #compute PCA
            if (is.null(precomputed)) {
              mdsdata = limma::plotMDS(SummarizedExperiment::assay(edata), plot = FALSE)
            } else {
              mdsdata = checkPrecomputedMDS(edata, precomputed)
            }

            #extract sample data
            sdata = BiocGenerics::as.data.frame(SummarizedExperiment::colData(edata), optional = TRUE)
            #create data structure
            drdf = pdataMDS_intl(mdsdata, dims)
            p1 = plotDR_intl(drdf, sdata, rl, ...)

            return(p1)
          })

calcPCA <- function(edata, dims) {
  maxdim = max(dims)
  if (maxdim < ncol(edata)) {
    pcdata = scater::calculatePCA(edata, ncomponents = maxdim)
  } else {
    pcdata = stats::prcomp(t(edata))
    pcdata = checkPrecomputedPCA(edata, pcdata)
  }

  return(pcdata)
}

checkPrecomputedPCA <- function(edata, pcdata) {
  if (is(pcdata, 'prcomp')) {
    stopifnot(all(rownames(pcdata$x) == colnames(edata)))

    #prepare results
    pvar = pcdata$sdev ^ 2 / sum(pcdata$sdev ^ 2) * 100
    rot = pcdata$rotation
    pcdata = pcdata$x
    attr(pcdata, 'percentVar') = pvar
    attr(pcdata, 'rotation') = rot
  } else if(is(pcdata, 'matrix')) {
    stopifnot(all(rownames(pcdata) == colnames(edata)))
    stopifnot(c('dim', 'dimnames', 'varExplained', 'percentVar', 'rotation') %in% names(attributes(pcdata)))
  } else {
    stop('provide results from prcomp or scater::calculatePCA')
  }

  return(pcdata)
}

checkPrecomputedMDS <- function(edata, mdsdata) {
  stopifnot(
    is(mdsdata, 'MDS'),
    all(rownames(mdsdata$distance.matrix.squared) == colnames(edata))
  )
  return(mdsdata)
}

#plot data preparation using PCA results
pdataPC_intl <- function(pcdata, dims, relabel = TRUE) {
  stopifnot(length(dims) == 2)

  #check sample names
  rnames = rownames(pcdata)
  if (is.null(rnames))
    rnames = 1:nrow(pcdata)

  plotdf = data.frame(
    'RestoolsMtchID' = rnames,
    x = pcdata[, dims[1]],
    y = pcdata[, dims[2]]
  )
  if (relabel) {
    pca_prop = round(attr(pcdata, 'percentVar'), digits = 2)
    pca_labs = paste0('PC', 1:length(pca_prop), ' (', pca_prop, '%)')
    colnames(plotdf)[-1] = pca_labs[dims]
  }

  return(plotdf)
}

#plot data preparation using MDS results
pdataMDS_intl <- function(mdsdata, dims) {
  stopifnot(length(dims) == 2)

  lambda = pmax(mdsdata$eigen.values, 0)
  plotdf = data.frame(
    'RestoolsMtchID' = rownames(mdsdata$distance.matrix.squared),
    x = mdsdata$eigen.vectors[, dims[1]] * lambda[dims[1]],
    y = mdsdata$eigen.vectors[, dims[2]] * lambda[dims[2]]
  )
  pca_labs = paste0('Leading logFC dim ', 1:ncol(mdsdata$eigen.vectors))
  pca_labs = paste0(pca_labs, ' (', round(mdsdata$var.explained * 100, digits = 2) , '%)')
  colnames(plotdf)[-1] = pca_labs[dims]

  return(plotdf)
}

#add colour annotation
addSampleAnnot <- function(plotdf, sdata) {
  plotdf = cbind(plotdf, sdata[plotdf$RestoolsMtchID, , drop = FALSE])
  return(plotdf)
}

plotDR_intl <- function(drdf, sdata, rl, ...) {
  #annotate samples
  plotdf = addSampleAnnot(drdf, sdata)

  #extract aes
  aesmap = rlang::enquos(...)
  #compute plot
  aesmap = aesmap[!names(aesmap) %in% c('x', 'y')] #remove x,y mappings if present

  #split aes params into those that are not aes i.e. static parametrisation
  if (length(aesmap) > 0) {
    is_aes = sapply(aesmap, rlang::quo_is_symbolic)
    defaultmap = lapply(aesmap[!is_aes], rlang::eval_tidy)
    aesmap = aesmap[is_aes]
  } else {
    defaultmap = list()
  }

  # aes requires x & y to be explicit: https://github.com/tidyverse/ggplot2/issues/3176
  x = rlang::sym(colnames(plotdf)[[2]])
  y = rlang::sym(colnames(plotdf)[[3]])

  #add size if not present
  if (is.null(defaultmap$size))
    defaultmap$size = 2

  # tidystyle recommends no explicit return statements at end of functions
  ggplot2::ggplot(plotdf, ggplot2::aes(!!x, !!y, !!!aesmap)) +
    ggplot2::geom_point() +
    # ggplot2::update_geom_defaults('point', defaultmap) +
    bhuvad_theme(rl)
}
