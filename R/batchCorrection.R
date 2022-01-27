#' Get negative control genes from the data
#'
#' @param spe_object A Spatial Experiment object.
#' @param n_assay Integer to indicate the nth count table in the assay(spe) to be used.
#' @param batch_name Column name indicating batches.
#' @param top_n Integer indicate how many genes to be included as negative control genes.
#'
#' @return A Spatial Experiment object.
#' @export
chooseNegCTRLgenes <- function(spe_object, n_assay = 2, batch_name = "SlideName", top_n = 200){

  . = samples = count = sd = m = cv = rowname = mean_zscore = NULL

  stopifnot(is.numeric(n_assay))
  stopifnot(n_assay <= length(spe_object@assays))
  stopifnot(batch_name %in% colnames(SummarizedExperiment::colData(spe_object)))

  spe <- spe_object

  gene_with_mzscore <- SummarizedExperiment::assay(spe, 2) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tidyr::gather(samples, count, -rowname) %>%
    left_join(SummarizedExperiment::colData(spe) %>%
                as.data.frame() %>%
                dplyr::select(c(batch_name)) %>%
                rownames_to_column(),
              by = c("samples"="rowname")) %>%
    split( ., f = .[,all_of(batch_name)]) %>%
    lapply(., function(x)
      x %>%
        tidyr::spread(samples, count) %>%
        dplyr::select(-batch_name) %>%
        column_to_rownames("rowname") %>%
        mutate(sd = apply(., 1, stats::sd),
               m = rowMeans(.),
               cv = log(100*sqrt(exp(sd^2)-1))) %>%
        dplyr::select(cv)) %>%
    bind_cols() %>%
    magrittr::set_colnames(paste0("cv",seq(ncol(.)))) %>%
    scale() %>%
    as.data.frame() %>%
    mutate(mean_zscore = rowMeans(.)) %>%
    dplyr::select(mean_zscore)

  SummarizedExperiment::rowData(spe)$mean_zscore <- gene_with_mzscore[rownames(spe),]
  SummarizedExperiment::rowData(spe)$mean_expr <- SummarizedExperiment::assay(spe, 2) %>%
    as.data.frame() %>%
    mutate(m = rowMeans(.)) %>%
    .$m

  negative.ctrl.genes <- gene_with_mzscore %>%
    arrange(mean_zscore) %>%
    rownames() %>%
    .[1:top_n]

  S4Vectors::metadata(spe)$negGenes <- negative.ctrl.genes

  return(spe)
}

#' RUV4 normalisation
#'
#' @param spe_object A Spatial Experiment object.
#' @param k The number of unwanted factors to use. Can be 0.
#' @param factors Column name(s) to indicate the factors of interest.
#' @param negctrlGenes Negative control genes.
#'
#' @return A Spatial Experiment object.
#' @export
#'
run_ruv4 <- function(spe_object, k, factors, negctrlGenes){

  spe <- spe_object

  stopifnot(is.numeric(k) & k>0)

  tmat <- spe@assays@data$logcounts %>% as.matrix() %>% t
  factorOfInterest <- SummarizedExperiment::colData(spe) %>% as.data.frame() %>% dplyr::select(all_of(factors))

  test <- ruv::design.matrix(factorOfInterest)

  ruv.out <- ruv::RUV4(tmat,
                  test,
                  ctl = rownames(spe) %in% negctrlGenes,
                  k = k, Z = NULL)

  ruv_w <- ruv.out$W %>%
    as.data.frame() %>%
    magrittr::set_colnames(paste0("ruv_W",seq(k)))

  for(i in seq(ncol(ruv_w))){
    n <- colnames(ruv_w)[i]
    #print(n)
    SummarizedExperiment::colData(spe)[,n] <- ruv_w[,i]
  }

  summary <- ruv::ruv_summary(tmat, ruv.out)

  ruv_norm_count <- ruv::ruv_residuals(summary, type = "adjusted.Y") %>% t

  spe@assays@data$logcounts <- ruv_norm_count

  return(spe)

}


#' Combat normalisation
#'
#' @param spe_object A Spatial Experiment object.
#' @param n_assay Integer, choose the assay from spe_object.
#' @param batch A vector indicating batches.
#' @param bio_factor A vector indicating biology.
#'
#' @return A Spatial Experiment object.
#' @export
run_combat <- function(spe_object, n_assay, batch, bio_factor){

  stopifnot(is.numeric(n_assay))
  stopifnot(n_assay <= length(spe_object@assays))
  stopifnot(length(batch)==ncol(spe_object))
  stopifnot(length(bio_factor)==ncol(spe_object))

  corrected_data <- sva::ComBat_seq(as.matrix(SummarizedExperiment::assay(spe_object, n_assay)),
                                    batch = batch, group = bio_factor)

  spe_combat <- spe_object

  spe_combat@assays@data$logcounts <- corrected_data

  return(spe_combat)
}

