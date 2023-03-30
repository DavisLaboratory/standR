
#' Preparing the inputs for SpatialDecon for doing deconvolution on spatial data
#'
#' @param spe SpatialExperiment object.
#' @param assay2use The name of the assay to use. By default is logcounts.
#' @param negProbeName The name of the negative probe gene. By default is NegProbe-WTX.
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
prepareSpatialDecon <- function(spe, assay2use = "logcounts", negProbeName = "NegProbe-WTX"){
  
  norm <- assay(spe, assay2use)
  
  if (!(negProbeName %in% rownames(norm))){
    stop(paste0(negProbeName, " must be included in the dataset. Perhaps specify rmNegProbe=TRUE when using readGeoMx"))
  }
  
  bg <- sweep(norm*0, 2, norm[negProbeName,], "+")
  
  out <- list(norm, bg)
  
  names(out) <- c("normCount","backGround")
  
  return(out)
}
