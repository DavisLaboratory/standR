#' Get negative control genes from each batch of the data
#'
#' @param spe A Spatial Experiment object.
#' @param n_assay Integer to indicate the nth count table in the assay(spe) to be used.
#' @param batch_name Column name indicating batches.
#' @param top_n Integer indicate how many genes to be included as negative control genes.
#'
#' @return A Spatial Experiment object, conatining negative control genes in the metadata.
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#'
#' spe <- findNCGs(dkd_spe_subset, top_n = 100)
#' S4Vectors::metadata(spe)$NCGs
#'
findNCGs <- function(spe, n_assay = 2, batch_name = "SlideName", top_n = 200) {
  stopifnot(is.numeric(n_assay))
  stopifnot(n_assay <= length(SummarizedExperiment::assayNames(spe)))
  stopifnot(batch_name %in% colnames(colData(spe)))
  stopifnot(top_n <= nrow(spe))

  # compute coefficient of variance for each batch
  gene_with_mzscore_list <- assay(spe, "logcounts") |>
    as.data.frame() |>
    rownames_to_column() |>
    tidyr::gather(samples, count, -rowname) |>
    left_join(colData(spe) |>
      as.data.frame(optional = TRUE) |>
      dplyr::select(c(batch_name)) |>
      rownames_to_column(),
    by = c("samples" = "rowname")
    ) |>
    (\(.) split(., f = .[, all_of(batch_name)]))() # split data into list of batches
    
  gene_with_mzscore <- lapply(gene_with_mzscore_list, function(x) {
      y <- x |>
        tidyr::spread(samples, count) |>
        dplyr::select(-batch_name) |>
        column_to_rownames("rowname")
      sd <- apply(y, 1, stats::sd)
      m <- rowMeans(y)
      cv <- log(100 * sqrt(exp(sd^2) - 1))
      return(data.frame(cv))
    }) |>
    bind_cols()
  colnames(gene_with_mzscore) <- paste0("cv", seq(ncol(gene_with_mzscore)))
  gene_with_mzscore <- scale(gene_with_mzscore) |>  # compute z-score
    as.data.frame()
  gene_with_mzscore$mean_zscore <- rowMeans(gene_with_mzscore)
  gene_with_mzscore <- gene_with_mzscore |> dplyr::select(mean_zscore)

  rowData(spe)$mean_zscore <- gene_with_mzscore[rownames(spe), ]
  mean_expr <- rowMeans(assay(spe, "logcounts") |> # get mean expression
                          as.data.frame())
  
  rowData(spe)$mean_expr <- mean_expr

  negative.ctrl.genes <- gene_with_mzscore |> # arrange by z-score, top N genes as negative control genes
    arrange(mean_zscore) |>
    rownames() |>
    (\(.) .[seq(top_n)])()

  S4Vectors::metadata(spe)$NCGs <- negative.ctrl.genes

  return(spe)
}


#' Batch correction for GeoMX data
#'
#' @param spe A Spatial Experiment object.
#' @param k The number of unwanted factors to use. Can be 0. This is required for the RUV4 method.
#' @param factors Column name(s) to indicate the factors of interest. This is required for the RUV4 method.
#' @param NCGs Negative control genes. This is required for the RUV4 method.
#' @param n_assay Integer to indicate the nth count table in the assay(spe) to be used.
#' @param batch A vector indicating batches. This is required for the Limma method.
#' @param batch2 A vector indicating the second series of batches. This is specific for the Limma method.
#' @param covariates A matrix or vector of numeric covariates to be adjusted for.
#' @param design A design matrix relating to treatment conditions to be preserved, can be generated using `stats::model.matrix` function with all biological factors included.
#' @param method Can be either RUV4 or Limma or RUVg, by default is RUV4.
#' @param isLog Logical vector, indicating if the count table is log or not.
#'
#' @return A Spatial Experiment object, containing the ruv4-normalized count and normalization factor.
#' @export
#'
#' @references Gagnon-Bartsch, J. A., Jacob, L., & Speed, T. P. (2013). Removing unwanted variation from high dimensional data with negative controls. Berkeley: Tech Reports from Dep Stat Univ California, 1-112.
#' @references Ritchie, M. E., Phipson, B., Wu, D. I., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic acids research, 43(7), e47-e47.

#'
#' @examples
#' data("dkd_spe_subset")
#' spe <- findNCGs(dkd_spe_subset, top_n = 100)
#' spe_ruv <- geomxBatchCorrection(spe,
#'   k = 3,
#'   factors = c("disease_status", "region"),
#'   NCGs = S4Vectors::metadata(spe)$NCGs
#' )
#'
geomxBatchCorrection <- function(spe, k, factors, NCGs, n_assay = 2,
                                 batch = NULL, batch2 = NULL, covariates = NULL, design = matrix(1, ncol(spe), 1),
                                 method = c("RUV4", "Limma", "RUVg"), isLog = TRUE) {
  if (length(method) == 2) {
    method <- "RUV4"
  } else {
    stopifnot(method %in% c("RUV4", "Limma", "RUVg"))
    stopifnot(length(method) == 1)
  }

  if (method == "RUV4") {
    k <- as.integer(k)

    stopifnot(k > 0)

    # get count matrix, and transpose
    tmat <- assay(spe, n_assay) |>
      as.matrix() |>
      t()

    stopifnot(all(factors %in% colnames(colData(spe))))

    # get factor of interest matrix
    factorOfInterest <- colData(spe) |>
      as.data.frame(optional = TRUE) |>
      dplyr::select(all_of(factors))

    test <- ruv::design.matrix(factorOfInterest)

    # run ruv4
    ruv.out <- ruv::RUV4(tmat,
      test,
      ctl = rownames(spe) %in% NCGs,
      k = k, Z = NULL
    )

    # store results
    ruv_w <- ruv.out$W |>
      as.data.frame()
    colnames(ruv_w) <- paste0("ruv_W", seq(k))

    for (i in seq(ncol(ruv_w))) {
      n <- colnames(ruv_w)[i]
      colData(spe)[, n] <- ruv_w[, i]
    }

    summary <- ruv::ruv_summary(tmat, ruv.out)

    ruv_norm_count <- ruv::ruv_residuals(summary, type = "adjusted.Y") |> t()

    assay(spe, "logcounts") <- ruv_norm_count[rownames(spe), colnames(spe)]
  } else if (method == "Limma") {
    assay(spe, "logcounts") <- limma::removeBatchEffect(assay(spe, n_assay),
      batch = batch, batch2 = batch2,
      covariates = covariates,
      design = design
    ) |>
      (\(.) .[rownames(spe), colnames(spe)])()
  } else if (method == "RUVg") {
    k <- as.integer(k)
    
    stopifnot(k > 0)
    
    # get count matrix, and transpose
    tmat <- assay(spe, n_assay) |>
      as.matrix()
    
    ruvg_out <- RUVSeq::RUVg(tmat, k, isLog = isLog, cIdx = NCGs)
    
    # store results
    ruv_w <- ruvg_out$W |>
      as.data.frame()
    colnames(ruv_w) <- paste0("ruv_W", seq(k))
    
    for (i in seq(ncol(ruv_w))) {
      n <- colnames(ruv_w)[i]
      colData(spe)[, n] <- ruv_w[, i]
    }
    
    assay(spe, "logcounts") <- ruvg_out$normalizedCounts[rownames(spe), colnames(spe)]
    
  }

  return(spe)
}


# calculate silhouette score
getSilhouette <- function(pca_object, spe, foiColumn) {
  # compute euclidean distance
  distm <- stats::dist(pca_object)

  # grouping
  label_by_factors <- colData(spe) |>
    as.data.frame(optional = TRUE) |>
    dplyr::select(all_of(foiColumn)) |>
    rownames_to_column() |>
    mutate(grp_id = as.numeric(factor(!!sym(foiColumn))))

  c <- label_by_factors$grp_id
  names(c) <- label_by_factors$rowname

  # calculate silhouette score
  si <- cluster::silhouette(c, distm)

  ss <- mean(si[, 3])

  return(list(ss, c))
}

#' Testing multiple K for RUV4 batch correction to find the best K.
#'
#' @param spe A Spatial Experiment object.
#' @param maxK Integer. The max k to test, will test k from 1 to maxK, by default is 10.
#' @param factor_of_int Column name(s) to indicate the factors of interest. This is required for the RUV4 method.
#' @param factor_batch Column name to indicate the batch.
#' @param NCGs Negative control genes. This is required for the RUV4 method.
#' @param point_size Numeric. Plotting parameter.
#' @param line_col Character. Plotting parameter.
#' @param point_col Character. Plotting parameter.
#' @param text_size Numeric. Plotting parameter.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data("dkd_spe_subset")
#' spe <- findNCGs(dkd_spe_subset, top_n = 100)
#' findBestK(spe,
#'   factor_of_int = c("disease_status"),
#'   factor_batch = "SlideName", NCGs = S4Vectors::metadata(spe)$NCGs
#' )
#'
findBestK <- function(spe, maxK = 10, factor_of_int, factor_batch, NCGs, point_size = 3, line_col = "black", point_col = "black", text_size = 13) {
  kdata <- data.frame()
  for (i in seq(maxK)) {
    spe_ruv <- geomxBatchCorrection(spe, k = i, factors = factor_of_int, NCGs = NCGs)
    pca_object <- calcPCA(assay(spe_ruv, "logcounts"), dims = c(1, 2))
    ss <- getSilhouette(pca_object, spe_ruv, foiColumn = factor_batch)[[1]]
    kdata <- rbind(kdata, data.frame("k" = i, "silhouette" = ss))
  }

  kdata |>
    ggplot(aes(k, silhouette)) +
    geom_path(col = line_col) +
    geom_point(size = point_size, col = point_col) +
    scale_x_continuous(breaks = function(x) int_breaks(x, n = 10)) +
    theme_bw() +
    xlab("K") +
    ylab("Silhouette Score (Batch cluster)") +
    theme(text = element_text(size = text_size))
}


# in findNCGs function
utils::globalVariables(c(".", "samples", "count", "sd", "m", "cv", "rowname", "mean_zscore"))

# in findBestK
utils::globalVariables(c("k", "silhouette"))
