#' Calculate statistics for evaluating batch correction
#'
#' @param spe_object A Spatial Experiment object.
#' @param foiColumn A column name indicating the factor of interest to be tested, can be biological factor or batch factor.
#' @param n_dimension The top n dimensions to be plotted
#' @param precomputed a dimensional reduction results from `stats::prcomp`.
#'   result in `reducedDims(object)` to plot. Default is NULL,
#'   we will compute for you.
#' @param assay a numeric or character, specifying the assay to use (for
#'   `SummarizedExperiment` and its derivative classes).
#'
#' @return A dataframe object containing the clustering evaluating statistics.
#' @export
#'
#' @examples
#' library(scater)
#' data("dkd_spe_subset")
#' spe <- scater::runPCA(dkd_spe_subset)
#' computeClusterEvalStats(spe, "SlideName")
#'
computeClusterEvalStats <- function(spe_object, foiColumn, precomputed = NULL,
                                    n_dimension = c(1,2), assay = 2){

  stopifnot(foiColumn %in% colnames(SummarizedExperiment::colData(spe_object)))

  # get PCA results

  #compute PCA
  if (is.null(precomputed)) {
    pcdata = calcPCA(SummarizedExperiment::assay(spe_object, assay), n_dimension)
  } else {
    pcdata = checkPrecomputedPCA(spe_object, precomputed)
  }

  pca_object <- pcdata

  # compute euclidean distance
  distm <- stats::dist(pca_object)

  # grouping
  label_by_factors <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame(optional = TRUE) %>%
    dplyr::select(all_of(foiColumn)) %>%
    rownames_to_column() %>%
    mutate(grp_id = as.numeric(factor(!!sym(foiColumn)))) %>%
    mutate(sample_id = row_number())

  c <- label_by_factors$grp_id
  names(c) <- label_by_factors$sample_id

  # calculate silhouette score
  si <- cluster::silhouette(c, distm)

  ss <- mean(si[,3])

  # prepare for other stats to be computed
  k <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame(optional = TRUE) %>%
    dplyr::select(all_of(foiColumn))
  k <- k[,1] %>%
    unique() %>%
    length()

  # clustering
  set.seed(119)

  kc <- stats::kmeans(pca_object[,1:2], k)

  km_clusters <- kc$cluster

  types = NULL

  # compute adjrand and jaccard
  df_out <- mclustcomp::mclustcomp(km_clusters, c) %>%
    filter(types %in% c("adjrand","jaccard"))

  df_out[3,] <- c("Silhouette Coefficient",abs(ss))

  df_out <- df_out %>%
    mutate(types = ifelse(types == "adjrand", "Adjusted Rand Index",
                          ifelse(types == "jaccard", "Jaccard Similarity Coefficient", types)))

  return(df_out)
}


#' Compare and evaluate different batch corrected data with plotting.
#'
#' @param spe_list A list of Spatial Experiment object.
#' @param bio_feature_name The common biological variation name.
#' @param batch_feature_name The common batch variation name.
#' @param data_names Data names.
#' @param colors Color values of filing the bars.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' library(scater)
#' data("dkd_spe_subset")
#' spe <- scater::runPCA(dkd_spe_subset)
#' spe2 <- spe
#' spe3 <- spe
#' plotClusterEvalStats(list(spe, spe2, spe3), bio_feature_name = "region",
#'                            batch_feature_name = "SlideName", c("test1","test2","test3"))
plotClusterEvalStats <- function(spe_list, bio_feature_name, batch_feature_name,
                                 data_names, colors = NA){

  from = scores = types = NULL


  # get stat for bio factor
  stat_bio <- lapply(spe_list ,function(x){
    computeClusterEvalStats(x, bio_feature_name)
  }) %>%
    bind_rows() %>%
    mutate(from = rep(data_names, each = 3))

  # get stat for batch factor
  stat_batch <- lapply(spe_list ,function(x){
    computeClusterEvalStats(x, batch_feature_name)
  }) %>%
    bind_rows() %>%
    mutate(from = rep(data_names, each = 3))

  p_bio <- stat_bio %>%
    mutate(from = factor(from, levels = data_names)) %>%
    mutate(scores = as.numeric(scores)) %>%
    ggplot(aes(from, scores, fill = types)) +
    geom_bar(stat = "identity", col = "black", width = .7) +
    facet_wrap(~types, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("Count data") +
    ylab("Scores") +
    ggtitle("Biology")

  p_batch <- stat_batch %>%
    mutate(from = factor(from, levels = data_names)) %>%
    mutate(scores = as.numeric(scores)) %>%
    ggplot(aes(from, scores, fill = types)) +
    geom_bar(stat = "identity", col = "black", width = .7) +
    facet_wrap(~types, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("Count data") +
    ylab("Scores") +
    ggtitle("Batch")

  if(!is.na(colors)){
    p_bio <- p_bio + scale_fill_manual(values = colors)
    p_batch <- p_batch + scale_fill_manual(values = colors)
  }

  print(p_bio + p_batch + patchwork::plot_layout(1, 2))
}
