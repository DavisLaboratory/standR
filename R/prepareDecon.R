
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
  
  if (!all(negProbeName %in% rownames(norm))){
    stop(paste0(negProbeName, " must be included in the dataset. Perhaps specify rmNegProbe=TRUE when using readGeoMx"))
  }
  
  if (length(negProbeName) == 1){
    bg <- sweep(norm*0, 2, norm[negProbeName,], "+")
  } else {
    if (nrow(norm) != length(pool)){
      stop("length(pool) should be equal to nrow(spe).")
    }
    
    bg <- norm*0
    
    for (p in unique(pool)){
      genes_p <- rownames(norm)[pool == p]
      negp_p <- base::intersect(negProbeName, genes_p)
      
      stopifnot("Negative probes not found in one of the pool." = length(negp_p) > 0)
      
      mean_bg_p <- colMeans(norm[negp_p, , drop = FALSE])
      
      bg[pool == p, ] <- sweep(bg[pool == p, ], 2, mean_bg_p, "+")
    }
  }
  
  out <- list(norm, bg)
  
  names(out) <- c("normCount","backGround")
  
  return(out)
}
