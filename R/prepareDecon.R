
#' Preparing the inputs for SpatialDecon for doing deconvolution on spatial data
#'
#' @param spe SpatialExperiment object.
#' @param assay2use The name of the assay to use. By default is logcounts.
#' @param negProbeName The name of the negative probe gene. By default is NegProbe-WTX.
#' @param pool A vector indicates the pools of the genes. This is required when there are more than one Negative Probes.
#'
#' @return A list of two dataframes. The first data.frame is the normalised count, the second data.frame is the background for the data.
#' @export
#'
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' 
#' query(eh, "standR")
#' countFile <- eh[["EH7364"]]
#' sampleAnnoFile <- eh[["EH7365"]]
#' 
#' spe <- readGeoMx(countFile, sampleAnnoFile, rmNegProbe = FALSE)
#' 
#' out <- prepareSpatialDecon(spe)
#' 
prepareSpatialDecon <- function(spe, assay2use = "logcounts", negProbeName = "NegProbe-WTX", pool = NA){

  norm <- assay(spe, assay2use)

  neg_norm <- .get_neg_probe_values(spe, norm, assay2use, negProbeName)

  if (is.null(neg_norm)){
    stop(
      paste0(
        paste(negProbeName, collapse = ", "),
        " must be included in the assay rows or stored in metadata(spe)$NegProbes. ",
        "If using readGeoMx(), set rmNegProbe=FALSE to keep negative probes in assays."
      )
    )
  }

  if (all(is.na(pool))){
    mean_bg <- colMeans(neg_norm, na.rm = TRUE)
    bg <- sweep(norm*0, 2, mean_bg, "+")
  } else {
    if (nrow(norm) != length(pool)){
      stop("length(pool) should be equal to nrow(spe).")
    }

    bg <- norm*0

    for (p in unique(pool)){
      genes_p <- rownames(norm)[pool == p]
      negp_p <- base::intersect(rownames(neg_norm), genes_p)

      stopifnot("Negative probes not found in one of the pool." = length(negp_p) > 0)

      mean_bg_p <- colMeans(neg_norm[negp_p, , drop = FALSE])

      bg[pool == p, ] <- sweep(bg[pool == p, ], 2, mean_bg_p, "+")
    }
  }

  out <- list(norm, bg)

  names(out) <- c("normCount","backGround")

  return(out)
}

.get_neg_probe_values <- function(spe, norm, assay2use, negProbeName) {
  neg_rows <- .match_neg_probe_rows(rownames(norm), negProbeName)
  if (length(neg_rows) > 0) {
    return(norm[neg_rows, , drop = FALSE])
  }

  neg_counts <- S4Vectors::metadata(spe)$NegProbes
  if (is.null(neg_counts)) {
    return(NULL)
  }

  neg_rows <- .match_neg_probe_rows(rownames(neg_counts), negProbeName)
  if (length(neg_rows) == 0) {
    return(NULL)
  }
  neg_counts <- as.matrix(neg_counts[neg_rows, , drop = FALSE])
  neg_counts <- neg_counts[, colnames(norm), drop = FALSE]

  if (identical(assay2use, "counts")) {
    return(neg_counts)
  }

  if (identical(assay2use, "logcounts")) {
    counts <- assay(spe, "counts")
    combined_counts <- rbind(counts, neg_counts)
    return(edgeR::cpm(combined_counts, log = TRUE)[rownames(neg_counts), , drop = FALSE])
  }

  stop(
    "Negative probes stored in metadata(spe)$NegProbes are raw counts; ",
    "use assay2use='counts' or assay2use='logcounts', or keep negative probes in the selected assay."
  )
}

.match_neg_probe_rows <- function(row_names, negProbeName) {
  exact <- which(row_names %in% negProbeName)
  if (length(exact) > 0) {
    return(exact)
  }

  if (length(negProbeName) != 1) {
    return(integer())
  }

  which(startsWith(row_names, paste0(negProbeName, "_")))
}
