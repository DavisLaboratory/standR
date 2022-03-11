#' Compute and plot relative log expression (RLE) values of gene expression data
#'
#' @param ordannots variables or computations to sort samples by (tidy style).
#'
#' @inheritParams plotPCA
#' @return a ggplot2 object, containing the RLE plot.
#' @export
#'
#' @examples
#' se <- emtdata::cursons2018_se()
#' dge <- emtdata::asDGEList(se)
#' plotRLE(dge,
#'   colour = Subline, lty = Treatment, lwd = 2, ordannots =
#'     c(Subline, Treatment)
#' )
#'
setGeneric(
  "plotRLE",
  function(edata,
           ordannots = c(),
           ...) {
    standardGeneric("plotRLE")
  }
)

#' @rdname plotRLE
setMethod(
  "plotRLE",
  signature("DGEList", "ANY"),
  function(edata, ordannots, ...) {
    # extract sample data
    sdata <- edata$samples
    # extract expression data (and transform)
    edata <- edgeR::cpm(edata, log = TRUE)
    # create data structure
    samporder <- orderSamples(sdata, rlang::enquo(ordannots))
    rledf <- pdataRLE_intl(edata, samporder)
    p1 <- plotRLE_intl(rledf, sdata, isSCE = FALSE, ...)

    return(p1)
  }
)

#' @rdname plotRLE
setMethod(
  "plotRLE",
  signature("ExpressionSet", "ANY"),
  function(edata, ordannots, ...) {
    # extract sample data
    sdata <- Biobase::pData(edata)
    # extract expression data (and transform)
    edata <- Biobase::exprs(edata)
    # create data structure
    samporder <- orderSamples(sdata, rlang::enquo(ordannots))
    rledf <- pdataRLE_intl(edata, samporder)
    p1 <- plotRLE_intl(rledf, sdata, isSCE = FALSE, ...)

    return(p1)
  }
)

#' @rdname plotRLE
setMethod(
  "plotRLE",
  signature("SummarizedExperiment", "ANY"),
  function(edata, ordannots, assay = 1, ...) {
    isSCE <- is(edata, "SingleCellExperiment")
    
    # extract sample data
    if (is(edata, "ExperimentList")) {
      sdata <- BiocGenerics::as.data.frame(colFun(edata, experimentData = TRUE), optional = TRUE)
    } else {
      sdata <- BiocGenerics::as.data.frame(colFun(edata), optional = TRUE)
    }

    # extract expression data (and transform)
    edata <- SummarizedExperiment::assay(edata, i = assay)
    # create data structure
    samporder <- orderSamples(sdata, rlang::enquo(ordannots))
    rledf <- pdataRLE_intl(edata, samporder)
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
  rledf$x <- 1:nrow(rledf)
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
    is_aes <- sapply(aesmap, rlang::quo_is_symbolic)
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
      # ggplot2::update_geom_defaults('boxplot', defaultmap) +
      bhuvad_theme(rl)
  } else {
    p1 <- ggplot2::ggplot(plotdf, aes(x = x, y = middle, group = x, !!!aesmap)) +
      ggplot2::geom_boxplot(
        aes(ymin = ymin, ymax = ymax, upper = upper, middle = middle, lower = lower),
        stat = "identity"
      ) +
      ggplot2::geom_hline(yintercept = 0, colour = 2, lty = 2) +
      ggplot2::ylab("Relative log expression") +
      # ggplot2::update_geom_defaults('boxplot', defaultmap) +
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
