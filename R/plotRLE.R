#' Compute and plot relative log expression (RLE) values of gene expression data
#'
#' @param ordannots variables or computations to sort samples by (tidy style).
#'
#' @inheritParams plotPCA
#' @return a ggplot2 object, containing the RLE plot.
#' @export
#'
#' @examples
#' se = emtdata::cursons2018_se()
#' dge = emtdata::asDGEList(se)
#' plotRLE(dge, colour = Subline, lty = Treatment, lwd = 2, ordannots =
#' c(Subline, Treatment))
#'
setGeneric("plotRLE",
           function(edata,
                    ordannots = c(),
                    rl = 1,
                    ...) standardGeneric("plotRLE"))

#' @rdname plotRLE
setMethod("plotRLE",
          signature('DGEList','ANY', 'ANY'),
          function(edata, ordannots, rl, ...){
            #extract sample data
            sdata = edata$samples
            #extract expression data (and transform)
            edata = edgeR::cpm(edata, log = TRUE)
            #create data structure
            samporder = orderSamples(sdata, rlang::enquo(ordannots))
            rledf = pdataRLE_intl(edata, samporder)
            p1 = plotRLE_intl(rledf, sdata, rl, FALSE, ...)

            return(p1)
          })

#' @rdname plotRLE
setMethod("plotRLE",
          signature('ExpressionSet','ANY', 'ANY'),
          function(edata, ordannots, rl, ...){
            #extract sample data
            sdata = Biobase::pData(edata)
            #extract expression data (and transform)
            edata = Biobase::exprs(edata)
            #create data structure
            samporder = orderSamples(sdata, rlang::enquo(ordannots))
            rledf = pdataRLE_intl(edata, samporder)
            p1 = plotRLE_intl(rledf, sdata, rl, FALSE, ...)

            return(p1)
          })

#' @rdname plotRLE
setMethod("plotRLE",
          signature('SummarizedExperiment','ANY', 'ANY'),
          function(edata, ordannots, rl, ...){
            #extract sample data
            sdata = BiocGenerics::as.data.frame(SummarizedExperiment::colData(edata), optional = TRUE)
            #extract expression data (and transform)
            edata = SummarizedExperiment::assay(edata)
            #create data structure
            samporder = orderSamples(sdata, rlang::enquo(ordannots))
            rledf = pdataRLE_intl(edata, samporder)
            p1 = plotRLE_intl(rledf, sdata, rl, FALSE, ...)

            return(p1)
          })

#plot data preparation using MDS results
pdataRLE_intl <- function(emat, sampord) {
  #compute RLE
  rle = emat - Biobase::rowMedians(emat)
  #order samples
  rle = rle[, sampord]

  #compute boxplot
  rledf = t(apply(rle, 2, function(x) grDevices::boxplot.stats(x)$stats))
  rledf = as.data.frame(rledf)
  colnames(rledf) = c('ymin', 'lower', 'middle', 'upper', 'ymax')
  rledf$x = 1:nrow(rledf)
  rledf$RestoolsMtchID = rownames(rledf)

  return(rledf)
}

orderSamples <- function(sdata, ordannots) {
  #add sample IDs
  sdata$SampleOrderID = rownames(sdata)

  #order samples based on provided annotations
  sdata = sdata |>
    dplyr::group_by(dplyr::across(!!ordannots)) |>
    dplyr::arrange(.by_group = TRUE)

  return(sdata$SampleOrderID)
}

plotRLE_intl <- function(plotdf, sdata, rl, isSCE = FALSE, ...) {
  #constant - sample size at which standard plot becomes dense
  dense_thresh = 50

  #extract aes
  aesmap = rlang::enquos(...)

  #annotate samples
  plotdf = addSampleAnnot(plotdf, sdata)

  #compute plot
  aesmap = aesmap[!names(aesmap) %in% c('x', 'ymin', 'ymin', 'upper', 'middle', 'lower')] #remove fixed mappings if present

  #split aes params into those that are not aes i.e. static parametrisation
  if (length(aesmap) > 0) {
    is_aes = sapply(aesmap, rlang::quo_is_symbolic)
    defaultmap = lapply(aesmap[!is_aes], rlang::eval_tidy)
    aesmap = aesmap[is_aes]
  } else {
    defaultmap = list()
  }

  #build plot
  if (!isSCE){
    p1 = ggplot2::ggplot(plotdf, aes(x = x, y = middle, group = x, !!!aesmap)) +
      ggplot2::geom_boxplot(
        aes(ymin = ymin, ymax = ymax, upper = upper, middle = middle, lower = lower),
        stat = 'identity') +
      ggplot2::geom_hline(yintercept = 0, colour = 2, lty = 2) +
      ggplot2::ylab('Relative log expression') +
      ggplot2::update_geom_defaults('boxplot', defaultmap) +
      bhuvad_theme(rl) +
      ggplot2::theme(axis.text.x = element_blank())

    #update plot if too many samples are plot
    if (nrow(plotdf) > dense_thresh) {
      ## geom_point will inherit relevant aesthetics from top `aes`, include y=middle
      p1 = p1 + ggplot2::geom_point()
    }
  } else {
    #place-holder for SCE object code
  }

  return(p1)
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
