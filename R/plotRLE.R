#' Compute and plot relative log expression (RLE) values of gene expression data
#'
#' @param ordannots variables or computations to sort samples by (tidy style).
#'
#' @inheritParams plotPCA
#' @return a ggplot2 object, containing the RLE plot.
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#' plotRLE(dkd_spe_subset)
#'
setGeneric(
  "plotRLE",
  function(object,
           ordannots = c(),
           ...) {
    standardGeneric("plotRLE")
  }
)

#' @rdname plotRLE
setMethod(
  "plotRLE",
  signature("DGEList", "ANY"),
  function(object, ordannots, ...) {
    # extract sample data
    sdata <- object$samples
    # extract expression data (and transform)
    object <- edgeR::cpm(object, log = TRUE)
    # create data structure
    samporder <- orderSamples(sdata, rlang::enquo(ordannots))
    rledf <- pdataRLE_intl(object, samporder)
    p1 <- plotRLE_intl(rledf, sdata, isSCE = FALSE, ...)

    return(p1)
  }
)

#' @rdname plotRLE
setMethod(
  "plotRLE",
  signature("ExpressionSet", "ANY"),
  function(object, ordannots, ...) {
    # extract sample data
    sdata <- Biobase::pData(object)
    # extract expression data (and transform)
    object <- Biobase::exprs(object)
    # create data structure
    samporder <- orderSamples(sdata, rlang::enquo(ordannots))
    rledf <- pdataRLE_intl(object, samporder)
    p1 <- plotRLE_intl(rledf, sdata, isSCE = FALSE, ...)

    return(p1)
  }
)

#' @rdname plotRLE
setMethod(
  "plotRLE",
  signature("SummarizedExperiment", "ANY"),
  function(object, ordannots, assay = 1, ...) {
    isSCE <- is(object, "SingleCellExperiment")

    # extract sample data
    sdata <- BiocGenerics::as.data.frame(SummarizedExperiment::colData(object), optional = TRUE)

    # extract expression data (and transform)
    object <- SummarizedExperiment::assay(object, i = assay)
    # create data structure
    samporder <- orderSamples(sdata, rlang::enquo(ordannots))
    rledf <- pdataRLE_intl(object, samporder)
    p1 <- plotRLE_intl(rledf, sdata, isSCE = isSCE, ...)

    return(p1)
  }
)

# plot data preparation using MDS results
pdataRLE_intl <- function(emat, sampord) {
  # compute RLE
  rle <- emat - Biobase::rowMedians(as.matrix(emat))
  # order samples
  rle <- rle[, sampord]

  # compute boxplot
  rledf <- t(apply(rle, 2, function(x) grDevices::boxplot.stats(x)$stats))
  rledf <- as.data.frame(rledf)
  colnames(rledf) <- c("ymin", "lower", "middle", "upper", "ymax")
  rledf$x <- seq(nrow(rledf))
  rledf$RestoolsMtchID <- rownames(rledf)

  return(rledf)
}

plotRLE_intl <- function(plotdf, sdata, isSCE = FALSE, rl = 1, ...) {

  # constant - sample size at which standard plot becomes dense
  dense_thresh <- 50
  sce_thresh <- 1000

  # extract aes
  aesmap <- rlang::enquos(...)

  # annotate samples
  plotdf <- addSampleAnnot(plotdf, sdata)

  # compute plot
  aesmap <- aesmap[!names(aesmap) %in% c("x", "ymin", "ymax", "upper", "middle", "lower")] # remove fixed mappings if present

  # split aes params into those that are not aes i.e. static parametrisation
  if (length(aesmap) > 0) {
    is_aes <- vapply(aesmap, rlang::quo_is_symbolic, FUN.VALUE = logical(1))
    defaultmap <- lapply(aesmap[!is_aes], rlang::eval_tidy)
    aesmap <- aesmap[is_aes]
  } else {
    defaultmap <- list()
  }

  # build plot
  if (isSCE & nrow(plotdf) > sce_thresh) {
    p1 <- ggplot2::ggplot(plotdf, aes(x = upper - lower, y = middle, !!!aesmap)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, colour = 2, lty = 2) +
      ggplot2::labs(y = "Relative log expression (median)", x = "Relative log expression (IQR)") +
      do.call(ggplot2::geom_point, defaultmap) +
      bhuvad_theme(rl)
  } else {
    p1 <- ggplot2::ggplot(plotdf, aes(x = x, y = middle, group = x, !!!aesmap)) +
      ggplot2::geom_boxplot(
        aes(ymin = ymin, ymax = ymax, upper = upper, middle = middle, lower = lower),
        stat = "identity"
      ) +
      ggplot2::geom_hline(yintercept = 0, colour = 2, lty = 2) +
      ggplot2::ylab("Relative log expression") +
      do.call(ggplot2::geom_boxplot, defaultmap) +
      bhuvad_theme(rl) +
      ggplot2::theme(axis.text.x = element_blank())

    # update plot if too many samples are plot
    if (nrow(plotdf) > dense_thresh) {
      ## geom_point will inherit relevant aesthetics from top `aes`, include y=middle
      p1 <- p1 + ggplot2::geom_point()
    }
  }

  return(p1)
}


utils::globalVariables(c("upper", "lower", "middle", "x", "ymin", "ymax"))
