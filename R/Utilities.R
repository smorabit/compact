#' Calculate Log2 Fold Change of In-Silico Perturbation
#'
#' @description
#' Computes the log2 fold change (Log2FC), mean delta, and average baseline 
#' expression for genes following an in-silico perturbation. It compares the 
#' perturbed assay against the baseline observed expression.
#'
#' @param seurat_obj A Seurat object containing both the baseline and perturbed data.
#' @param perturbation_name Character. The name of the assay containing the perturbed expression data.
#' @param assay_obs Character. The name of the baseline assay. Default is \code{'RNA'}.
#' @param module Character. Optional. The specific hdWGCNA module name to analyze. If provided, the output will include a column indicating whether each gene is a top hub.
#' @param n_hubs Integer. The number of top hub genes to flag if \code{module} is provided. Default is \code{10}.
#' @param pseudocount Numeric. Pseudocount added to average expression before log2 transformation to prevent infinite values. Default is \code{1}.
#' @param layer_norm Character. The layer/layer containing normalized data. Default is \code{'data'}.
#' @param layer_counts Character. The layer/layer containing raw counts. Default is \code{'counts'}.
#'
#' @return A data.frame containing gene names, mean delta (raw counts difference), average baseline expression, Log2FC, and optionally module hub status.
#' 
#' @import Seurat
#' @importFrom Matrix rowMeans rowSums
#' @importFrom dplyr left_join rename mutate
#' @export
PerturbationLog2FC <- function(
    seurat_obj,
    perturbation_name,
    assay_obs = 'RNA',
    module = NULL,
    n_hubs = 10,
    pseudocount = 1,
    layer_norm = 'data',
    layer_counts = 'counts'
) {
    
    # 1. Determine features to analyze
    if (!is.null(module)) {

        # Assuming GetHubGenes is an exported/internal function in your package
        hub_df <- GetHubGenes(seurat_obj, n_hubs = Inf)
        hub_df <- subset(hub_df, module == module)
        
        top_hubs <- GetHubGenes(seurat_obj, n_hubs = n_hubs)
        top_hubs <- subset(top_hubs, module == module)$gene_name
        
        hub_df$hub <- ifelse(hub_df$gene_name %in% top_hubs, 'hub', 'other')
        features <- hub_df$gene_name
    } else {
        # If no module provided, use all genes present in the perturbed assay
        features <- rownames(seurat_obj[[perturbation_name]])
    }
    
    # 2. Extract matrices
    # Note: Using Seurat's GetAssayData (works for v4 and mostly v5 backwards compatibility)
    X_counts <- GetAssayData(seurat_obj, assay = assay_obs, layer = layer_counts)[features, ]
    X_norm <- GetAssayData(seurat_obj, assay = assay_obs, layer = layer_norm)[features, ]
    
    X_per_counts <- GetAssayData(seurat_obj, assay = perturbation_name, layer = layer_counts)[features, ]
    X_per_norm <- GetAssayData(seurat_obj, assay = perturbation_name, layer = layer_norm)[features, ]
    
    # 3. Fast Matrix calculations
    # Compute difference in raw counts
    delta_matrix <- X_per_counts - X_counts
    
    # Use Matrix::rowMeans for lightning-fast sparse matrix calculations
    mean_delta <- Matrix::rowMeans(delta_matrix)
    avg_exp <- Matrix::rowMeans(X_norm)
    avg_per_exp <- Matrix::rowMeans(X_per_norm)
    
    # Calculate Log2FC
    l2fc <- log2((avg_per_exp + pseudocount) / (avg_exp + pseudocount))
    
    # 4. Construct final dataframe
    plot_df <- data.frame(
        gene_name = features,
        mean_delta = mean_delta,
        avg_exp = avg_exp,
        log2fc = l2fc,
        stringsAsFactors = FALSE
    )
    
    # Merge hub info if a module was specified
    if (!is.null(module)) {
        plot_df <- dplyr::left_join(plot_df, hub_df, by = "gene_name")
    }
    
    return(plot_df)
}