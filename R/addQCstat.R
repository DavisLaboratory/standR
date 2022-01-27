
#' Add QC statistics to the Spatial Experiment object
#'
#' @param spe_object A spatial experiment object
#' @param sample_fraction Double. Genes with low count in more than this threshold of the samples will be removed. Default is 0.9
#' @param rm_genes Logical. Decide whether genes with low count in more than sample_fraction of the samples are removed from the dataset. Default is TRUE.
#' @param min_count Interger. Minimum read count to calculate count threshold. Default is 5.
#' @param addGeneinfo Logical. Decide if extra gene information are added to the object. Default is TRUE.
#' @param gene_coords Default is from_bioMart, meaning human gene information including genomic position and gene length will be added based on the information in bioMart. Note: genes that failed getting info from biomart will be removed. Alternatively, this can be a dataframe with 5 columns with colnames hgnc_symbol, ensembl_gene_id, start_position, end_position, geneLength.
#'
#' @return A spatial experiment object
#' @export
#'
addQCstat <- function(spe_object, sample_fraction = 0.9, rm_genes = TRUE, min_count = 5,
                      addGeneinfo = FALSE, gene_coords = "from_bioMart"){

  . = NULL
  spe <- spe_object

  stopifnot(nrow(spe)>0)
  stopifnot(ncol(spe)>0)
  stopifnot(sample_fraction<=1 & sample_fraction>=0)
  stopifnot(min_count>=0)

  ## Library size
  SummarizedExperiment::colData(spe)$lib_size <- SummarizedExperiment::assay(spe,1) %>% colSums()


  ## calculate lcpm threshold for a gene to be expressed
  L <- mean(SummarizedExperiment::colData(spe)$lib_size) * 1e-6
  M <- stats:: median(SummarizedExperiment::colData(spe)$lib_size) * 1e-6
  lcpm_threshold <- log2(min_count/M + 2/L)
  spe@metadata$lcpm_threshold <- lcpm_threshold

  ## Feature-wise QC & Filter
  ###### Genes with low count in more than the threshold of the samples will be removed
  genes_lowCount_overNsamples <- rowSums(SummarizedExperiment::assay(spe, 2) <= lcpm_threshold) >= round(ncol(SummarizedExperiment::assay(spe, 2))*sample_fraction)
  SummarizedExperiment::rowData(spe)$genes_lowCount_overNsamples <- genes_lowCount_overNsamples

  if(rm_genes == TRUE){
    spe@metadata$genes_rm_rawCount <- SummarizedExperiment::assay(spe,1)[SummarizedExperiment::rowData(spe)$genes_lowCount_overNsamples,]
    spe@metadata$genes_rm_logCPM <- SummarizedExperiment::assay(spe,2)[SummarizedExperiment::rowData(spe)$genes_lowCount_overNsamples,]
    spe <- spe[!SummarizedExperiment::rowData(spe)$genes_lowCount_overNsamples,]
  }

  ## Sample-wise QC & Filter
  SummarizedExperiment::colData(spe)$countOfLowEprGene <- colSums(SummarizedExperiment::assay(spe,2) <= lcpm_threshold)
  SummarizedExperiment::colData(spe)$percentOfLowEprGene <- colSums(SummarizedExperiment::assay(spe,2) <= lcpm_threshold)/nrow(SummarizedExperiment::assay(spe,2))*100


  ## add gene info

  if(isTRUE(addGeneinfo)){
    if(gene_coords == "from_bioMart"){
      message("Assuming it's human.. otherwise gene coordinates infomation must be provided as a dataframe to the gene_coords parameter, with 4 columns and colnames are hgnc_symbol, ensembl_gene_id, start_position, end_position, geneLength.")
      #requireNamespace(biomaRt)
      gene_list <- rownames(SummarizedExperiment::rowData(spe))
      human <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl", host="uswest.ensembl.org")
      gene_coords <- biomaRt::getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="hgnc_symbol", values=gene_list, mart=human) %>%
        dplyr::mutate(geneLength = .$end_position - .$start_position)

      g1 <- gene_coords$hgnc_symbol %>% table() %>% tibble::as_tibble() %>% dplyr::filter(.data$n==1) %>% .$.
      g2 <- gene_coords$hgnc_symbol %>% table() %>% tibble::as_tibble() %>% dplyr::filter(.data$n!=1) %>% .$.

      gene_coords1 <- gene_coords %>%
        dplyr::filter(.data$hgnc_symbol %in% g1)

      gene_coords2 <- gene_coords %>%
        dplyr::filter(.$hgnc_symbol %in% g2) %>%
        dplyr::group_split(.$hgnc_symbol) %>%
        lapply(function(x) x %>%
                 dplyr::arrange(-.$geneLength) %>%
                 .[1,]) %>%
        dplyr::bind_rows() %>%
        .[,1:5]

      gene_coords <- rbind(gene_coords1, gene_coords2) %>%
        dplyr::arrange(.$hgnc_symbol)

      gene_pass_annotation <- rownames(SummarizedExperiment::rowData(spe)) %>% toupper() %in% gene_coords$hgnc_symbol

      spe@metadata$genes_fail_annotation_rawCount <- spe[!gene_pass_annotation,]@assays@data$rawCounts
      spe@metadata$genes_fail_annotation_logCPM <- spe[!gene_pass_annotation,]@assays@data$logcounts

      spe <- spe[gene_pass_annotation,]

      SummarizedExperiment::rowData(spe) <- SummarizedExperiment::rowData(spe) %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        dplyr::left_join(gene_coords, by = c("rowname"="hgnc_symbol")) %>%
        tibble::column_to_rownames("rowname")
    } else if(is.data.frame(gene_coords)){
      stopifnot("hgnc_symbol" %in% colnames(gene_coords))
      SummarizedExperiment::rowData(spe) <- SummarizedExperiment::rowData(spe) %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        dplyr::left_join(gene_coords, by = c("rowname"="hgnc_symbol")) %>%
        tibble::column_to_rownames("rowname")
    }
  }

  return(spe)
}
