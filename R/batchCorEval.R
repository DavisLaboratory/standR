#' Calculate statistics for evaluating batch correction
#'
#' @param spe_object A Spatial Experiment object.
#' @param foiColumn Factors of interest.
#'
#' @return A Spatial Experiment object.
#' @export
#'
computeClusterEvalStats <- function(spe_object, foiColumn){

  pca_object <- SingleCellExperiment::reducedDim(spe_object, "PCA")

  distm <- stats::dist(pca_object)

  label_by_factors <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame() %>%
    dplyr::select(foiColumn) %>%
    rownames_to_column() %>%
    mutate(grp_id = as.numeric(factor(!!sym(foiColumn)))) %>%
    mutate(sample_id = row_number())

  c <- label_by_factors$grp_id
  names(c) <- label_by_factors$sample_id

  si <- cluster::silhouette(c, distm)

  ss <- mean(si[,3])

  k <- SummarizedExperiment::colData(spe_object) %>%
    as.data.frame() %>%
    dplyr::select(foiColumn)
  k <- k[,1] %>%
    unique() %>%
    length()

  set.seed(119)

  kc <- stats::kmeans(pca_object[,1:2], k)

  km_clusters <- kc$cluster

  types = NULL

  df_out <- mclustcomp::mclustcomp(km_clusters, c) %>%
    filter(types %in% c("adjrand","jaccard"))

  df_out[3,] <- c("Silhouette Coefficient",ss)

  df_out <- df_out %>%
    mutate(types = ifelse(types == "adjrand", "Adjusted Rand Index",
                          ifelse(types == "jaccard", "Jaccard Similarity Coefficient", types)))

  return(df_out)
}

#' Plot the evaluation statistics
#'
#' @param spe_list A Spatial Experiment object.
#' @param bio_feature_name Biological variation name.
#' @param batch_feature_name Batch variation name.
#' @param data_names Data names.
#'
#' @return A ggplot object.
#' @export
#'
plotClusterEvalStats <- function(spe_list, bio_feature_name, batch_feature_name,
                                 data_names){

  from = scores = types = NULL

  stat_bio <- lapply(spe_list ,function(x){
    computeClusterEvalStats(x, bio_feature_name)
  }) %>%
    bind_rows() %>%
    mutate(from = rep(data_names, each = 3))

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

  print(p_bio + p_batch + patchwork::plot_layout(1, 2))
}
