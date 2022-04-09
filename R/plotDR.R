# compuate PCA
calcPCA <- function(object, dims) {
  maxdim <- max(dims)
  if (requireNamespace("scater") & maxdim < ncol(object)) {
    pcdata <- scater::calculatePCA(object, ncomponents = maxdim)
  } else {
    pcdata <- stats::prcomp(t(object))
    pcdata <- checkPrecomputedPCA(object, pcdata)
  }

  return(pcdata)
}

#' Compute and plot the results of a PCA analysis on gene expression data
#'
#' @param object a DGEList, SummarizedExperiment or ExpressionSet object
#'   containing gene expression data.
#' @param dims a numeric, containing 2 values specifying the dimensions to plot.
#' @param assay a numeric or character, specifying the assay to use (for
#'   `SummarizedExperiment` and its derivative classes).
#' @param precomputed a dimensional reduction results from `stats::prcomp`.
#'   result in `reducedDims(object)` to plot.
#' @param textScale a numeric, specifying the relative scale factor to apply to text on
#'   the plot.
#' @param ... aesthetic mappings to pass to `ggplot2::aes_string()`.
#'
#' @return a ggplot2 object
#' @export
#'
#' @docType methods
#' @name drawPCA
#' @rdname drawPCA
#' @aliases drawPCA drawPCA,DGEList-method drawPCA,ExpressionSet-method
#' @aliases drawPCA,SummarizedExperiment-method drawPCA,SingleCellExperiment-method drawPCA,SpatialExperiment-method
#'
#' @examples
#' data("dkd_spe_subset")
#' drawPCA(dkd_spe_subset)
#'
setGeneric("drawPCA",
           function(object,
                    dims = c(1, 2),
                    ...) standardGeneric("drawPCA"))

setMethod(
  "drawPCA",
  signature("DGEList"),
  function(object, dims = c(1, 2), precomputed = NULL, textScale = 1, ...) {
    # compute PCA
    if (is.null(precomputed)) {
      pcdata <- calcPCA(edgeR::cpm(object, log = TRUE), dims)
    } else {
      pcdata <- checkPrecomputedPCA(object, precomputed)
    }

    # extract sample data
    sdata <- object$samples
    # create data structure
    drdf <- pdataPC_intl(pcdata, dims)
    p1 <- plotDR_intl(drdf, sdata, textScale, ...)

    return(p1)
  }
)

#' @rdname drawPCA
setMethod(
  "drawPCA",
  signature("ExpressionSet"),
  function(object, dims = c(1, 2), precomputed = NULL, textScale = 1, ...) {
    # compute PCA
    if (is.null(precomputed)) {
      pcdata <- calcPCA(Biobase::exprs(object), dims)
    } else {
      pcdata <- checkPrecomputedPCA(object, precomputed)
    }

    # extract sample data
    sdata <- Biobase::pData(object)
    # create data structure
    drdf <- pdataPC_intl(pcdata, dims)
    p1 <- plotDR_intl(drdf, sdata, textScale, ...)

    return(p1)
  }
)

.drawPCA_se <- function(object, dims = c(1, 2), assay = 1, precomputed = NULL, textScale = 1, ...) {
  # compute PCA
  if (is.null(precomputed)) {
    pcdata <- calcPCA(assay(object, assay), dims)
  } else {
    pcdata <- checkPrecomputedPCA(object, precomputed)
  }

  # extract sample data
  sdata <- BiocGenerics::as.data.frame(colData(object), optional = TRUE)


  # create data structure
  drdf <- pdataPC_intl(pcdata, dims)
  p1 <- plotDR_intl(drdf, sdata, textScale, ...)

  return(p1)
}

#' @rdname drawPCA
setMethod("drawPCA", signature("SummarizedExperiment"), .drawPCA_se)

#' @rdname drawPCA
setMethod("drawPCA", signature("SingleCellExperiment"), .drawPCA_se)

#' @rdname drawPCA
setMethod("drawPCA", signature("SpatialExperiment"), .drawPCA_se)

#' Compute and plot the results of any dimension reduction methods on gene expression data
#'
#' @param dimred a string or integer scalar indicating the reduced dimension
#'   result in `reducedDims(object)` to plot.
#'
#' @inheritParams drawPCA
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' library(scater)
#' data("dkd_spe_subset")
#' spe <- scater::runPCA(dkd_spe_subset)
#' plotDR(spe, dimred = "PCA")
#'
setGeneric(
  "plotDR",
  function(object,
           dims = c(1, 2),
           ...) {
    standardGeneric("plotDR")
  }
)
#' @rdname plotDR
setMethod(
  "plotDR",
  signature("SingleCellExperiment", "ANY"),
  function(object, dims, dimred = "PCA", textScale = 1, ...) {
    # compute PCA
    pcdata <- SingleCellExperiment::reducedDim(object, type = dimred)

    # extract sample data
    sdata <- BiocGenerics::as.data.frame(colData(object), optional = TRUE)

    # create data structure
    drdf <- pdataPC_intl(pcdata, dims, relabel = FALSE)
    p1 <- plotDR_intl(drdf, sdata, textScale, ...)

    return(p1)
  }
)

#' @rdname plotDR
setMethod(
  "plotDR",
  signature("SpatialExperiment", "ANY"),
  function(object, dims, dimred = "PCA", textScale = 1, ...) {
    # compute PCA
    pcdata <- SingleCellExperiment::reducedDim(object, type = dimred)

    # extract sample data
    sdata <- BiocGenerics::as.data.frame(colData(object), optional = TRUE)

    # create data structure
    drdf <- pdataPC_intl(pcdata, dims, relabel = FALSE)
    p1 <- plotDR_intl(drdf, sdata, textScale, ...)

    return(p1)
  }
)

#' Compute and plot the results of a MDS analysis on gene expression data
#'
#' @param precomputed a dimensional reduction results from either
#'   `limma::plotMDS`.
#' @param assay a numeric or character, specifying the assay to use (for
#'   `SummarizedExperiment` and its derivative classes).
#' @inheritParams drawPCA
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#' standR::plotMDS(dkd_spe_subset)
#'
setGeneric(
  "plotMDS",
  function(object,
           dims = c(1, 2),
           precomputed = NULL,
           textScale = 1, assay = 1,
           ...) {
    standardGeneric("plotMDS")
  }
)

#' @rdname plotMDS
setMethod(
  "plotMDS",
  signature("DGEList", "ANY", "ANY", "ANY"),
  function(object, dims, precomputed, textScale, ...) {
    # compute PCA
    if (is.null(precomputed)) {
      mdsdata <- limma::plotMDS(object, plot = FALSE)
    } else {
      mdsdata <- checkPrecomputedMDS(object, precomputed)
    }

    # extract sample data
    sdata <- object$samples
    # create data structure
    drdf <- pdataMDS_intl(mdsdata, dims)
    p1 <- plotDR_intl(drdf, sdata, textScale, ...)

    return(p1)
  }
)

#' @rdname plotMDS
setMethod(
  "plotMDS",
  signature("ExpressionSet", "ANY", "ANY", "ANY"),
  function(object, dims, precomputed, textScale, ...) {
    # compute PCA
    if (is.null(precomputed)) {
      mdsdata <- limma::plotMDS(object, plot = FALSE)
    } else {
      mdsdata <- checkPrecomputedMDS(object, precomputed)
    }

    # extract sample data
    sdata <- Biobase::pData(object)
    # create data structure
    drdf <- pdataMDS_intl(mdsdata, dims)
    p1 <- plotDR_intl(drdf, sdata, textScale, ...)

    return(p1)
  }
)


.plotMDS_se <- function(object, dims = c(1, 2), precomputed = NULL, textScale = 1, assay = 1, ...) {
  # compute PCA
  if (is.null(precomputed)) {
    mdsdata <- limma::plotMDS(assay(object, assay), plot = FALSE)
  } else {
    mdsdata <- checkPrecomputedMDS(object, precomputed)
  }

  # extract sample data
  sdata <- BiocGenerics::as.data.frame(colData(object), optional = TRUE)
  # create data structure
  drdf <- pdataMDS_intl(mdsdata, dims)
  p1 <- plotDR_intl(drdf, sdata, textScale, ...)

  return(p1)
}

#' @rdname plotMDS
setMethod("plotMDS", signature("SummarizedExperiment"), .plotMDS_se)

#' @rdname plotMDS
setMethod("plotMDS", signature("SingleCellExperiment"), .plotMDS_se)

#' @rdname plotMDS
setMethod("plotMDS", signature("SpatialExperiment"), .plotMDS_se)

checkPrecomputedPCA <- function(object, pcdata) {
  if (is(pcdata, "prcomp")) {
    stopifnot(all(rownames(pcdata$x) == colnames(object)))

    # prepare results
    pvar <- pcdata$sdev^2 / sum(pcdata$sdev^2) * 100
    rot <- pcdata$rotation
    pcdata <- pcdata$x
    attr(pcdata, "percentVar") <- pvar
    attr(pcdata, "rotation") <- rot
  } else if (is.matrix(pcdata)) {
    stopifnot(all(rownames(pcdata) == colnames(object)))
    stopifnot(c("dim", "dimnames", "varExplained", "percentVar", "rotation") %in% names(attributes(pcdata)))
  } else {
    stop("provide results from prcomp or scater::calculatePCA")
  }

  return(pcdata)
}

checkPrecomputedMDS <- function(object, mdsdata) {
  stopifnot(
    is(mdsdata, "MDS"),
    all(rownames(mdsdata$distance.matrix.squared) == colnames(object))
  )
  return(mdsdata)
}

# plot data preparation using PCA results
pdataPC_intl <- function(pcdata, dims, relabel = TRUE) {
  stopifnot(length(dims) == 2)

  # check sample names
  rnames <- rownames(pcdata)
  if (is.null(rnames)) {
    rnames <- seq(nrow(pcdata))
  }

  plotdf <- data.frame(
    "RestoolsMtchID" = rnames,
    x = pcdata[, dims[1]],
    y = pcdata[, dims[2]]
  )
  if (relabel) {
    pca_prop <- round(attr(pcdata, "percentVar"), digits = 2)
    pca_labs <- paste0("PC", seq(length(pca_prop)), " (", pca_prop, "%)")
    colnames(plotdf)[-1] <- pca_labs[dims]
  }

  return(plotdf)
}

# plot data preparation using MDS results
pdataMDS_intl <- function(mdsdata, dims) {
  stopifnot(length(dims) == 2)

  lambda <- pmax(mdsdata$eigen.values, 0)
  plotdf <- data.frame(
    "RestoolsMtchID" = rownames(mdsdata$distance.matrix.squared),
    x = mdsdata$eigen.vectors[, dims[1]] * lambda[dims[1]],
    y = mdsdata$eigen.vectors[, dims[2]] * lambda[dims[2]]
  )
  pca_labs <- paste0("Leading logFC dim ", seq(ncol(mdsdata$eigen.vectors)))
  pca_labs <- paste0(pca_labs, " (", round(mdsdata$var.explained * 100, digits = 2), "%)")
  colnames(plotdf)[-1] <- pca_labs[dims]

  return(plotdf)
}

# add colour annotation
addSampleAnnot <- function(plotdf, sdata) {
  plotdf <- cbind(plotdf, sdata[plotdf$RestoolsMtchID, , drop = FALSE])
  return(plotdf)
}

plotDR_intl <- function(drdf, sdata, textScale, ...) {
  # annotate samples
  plotdf <- addSampleAnnot(drdf, sdata)

  # extract aes
  aesmap <- rlang::enquos(...)
  # compute plot
  aesmap <- aesmap[!names(aesmap) %in% c("x", "y")] # remove x,y mappings if present

  # split aes params into those that are not aes i.e. static parametrisation
  if (length(aesmap) > 0) {
    is_aes <- vapply(aesmap, rlang::quo_is_symbolic, FUN.VALUE = logical(1))
    defaultmap <- lapply(aesmap[!is_aes], rlang::eval_tidy)
    aesmap <- aesmap[is_aes]
  } else {
    defaultmap <- list()
  }

  # aes requires x & y to be explicit: https://github.com/tidyverse/ggplot2/issues/3176
  x <- rlang::sym(colnames(plotdf)[[2]])
  y <- rlang::sym(colnames(plotdf)[[3]])

  # add size if not present
  if (is.null(defaultmap$size)) {
    defaultmap$size <- 2
  }

  # tidystyle recommends no explicit return statements at end of functions
  ggplot2::ggplot(plotdf, ggplot2::aes(!!x, !!y, !!!aesmap)) +
    ggplot2::geom_point() +
    do.call(ggplot2::geom_point,defaultmap) +
    bhuvad_theme(textScale)
}
