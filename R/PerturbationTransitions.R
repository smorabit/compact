
#' PerturbationTransitions
#'
#' This function computes cell-cell transition probabilities based on an in-silico perturbation experiment.
#' 
#' @return A Seurat object with the cell-cell transition probabilities stored in the Graphs slot of the Seurat object.
#' 
#' @param seurat_obj A Seurat object
#' @param perturbation_name A name for the in-silico perturbation that will be stored in the Seurat obejct
#' @param features Selected features to use for the transition probability calculation
#' @param graph Name of the cell-cell graph in the Graphs(seurat_obj)
#' @param corr_sigma A numeric scaling factor for the correlation matrix.
#' @param n_threads Number of threads for the correlation calculation
#' @param slot Slot to extract data for aggregation. Default = 'data'
#' @param assay Assay in seurat_obj containing expression information.
#'
#' @import Matrix
#' @import Seurat
#' @export
PerturbationTransitions <- function(
    seurat_obj,
    perturbation_name,
    features,
    graph, 
    use_velocyto = TRUE,
    use_graph_tp = FALSE,
    corr_sigma=0.05,
    n_threads=4,
    layer='data',
    slot='data',
    assay="RNA"
){

    # check graph exists in the Seurat object
    if(!(graph %in% names(seurat_obj@graphs))){
        stop(paste0("Graph '", graph, "' not found in seurat_obj@graphs. Available graphs: ",
                    paste(names(seurat_obj@graphs), collapse = ', ')))
    }

    # check observed assay exists
    if(!(assay %in% names(seurat_obj@assays))){
        stop(paste0("Assay '", assay, "' not found in seurat_obj@assays. Available assays: ",
                    paste(names(seurat_obj@assays), collapse = ', ')))
    }

    # check perturbation assay exists
    if(!(perturbation_name %in% names(seurat_obj@assays))){
        stop(paste0("Perturbation assay '", perturbation_name, "' not found in seurat_obj@assays. ",
                    "Available assays: ", paste(names(seurat_obj@assays), collapse = ', ')))
    }

    # check velocyto.R is available when use_velocyto = TRUE
    if(use_velocyto && !requireNamespace("velocyto.R", quietly = TRUE)){
        stop("use_velocyto = TRUE requires the 'velocyto.R' package. ",
             "Install it or set use_velocyto = FALSE to use the built-in SparseColDeltaCor.")
    }

    cell_graph <- Graphs(seurat_obj, slot=graph)
    cell_graph <- Matrix::Matrix(cell_graph)
    diag(cell_graph) <- 1

    # get the observed and the perturbed expression matrices:
    if(hdWGCNA::CheckSeurat5()){
        exp_obs <- GetAssayData(seurat_obj, assay=assay, layer=layer)
        exp_per <- GetAssayData(seurat_obj, assay=perturbation_name, layer=layer)
    }
    else{
        exp_obs <- GetAssayData(seurat_obj, assay=assay, slot=slot)
        exp_per <- GetAssayData(seurat_obj, assay=perturbation_name, slot=slot)
    }

    # check that requested features exist in both assays
    missing_obs <- setdiff(features, rownames(exp_obs))
    missing_per <- setdiff(features, rownames(exp_per))
    if(length(missing_obs) > 0){
        stop(paste0("The following features are not found in assay '", assay, "': ",
                    paste(missing_obs, collapse = ', ')))
    }
    if(length(missing_per) > 0){
        stop(paste0("The following features are not found in perturbation assay '",
                    perturbation_name, "': ", paste(missing_per, collapse = ', ')))
    }

    # subset by selected features
    exp_obs <- exp_obs[features,]
    exp_per <- exp_per[features,]

    # run the Velocyto colDeltaCor function
    if(use_velocyto){
        
        # convert to dense (only needed for velocyto)
        delta <- as.matrix(exp_per) - as.matrix(exp_obs)
        cc <- colDeltaCor_velocyto(exp_obs, delta, nthreads=n_threads)
    } else{
        delta <- exp_per - exp_obs
        cc <- SparseColDeltaCor(
          methods::as(exp_obs, "CsparseMatrix"),
          methods::as(delta, "CsparseMatrix"),
          methods::as(cell_graph, "CsparseMatrix")
        )
    }

    # fill the diagonal with zeros (cells won't transition to self)
    diag(cc) <- 0

    # compute transition probs between cells
    tp <- exp(cc / corr_sigma) * cell_graph
    col_sums <- Matrix::colSums(tp)

    # cells with zero or NA column sum (NA arises when SparseColDeltaCor
    # returns NA for cells whose delta is all-zero) cannot be normalized;
    # set those entries to 1 so the division leaves the all-zero column as-is
    zero_cols <- is.na(col_sums) | col_sums == 0
    if(any(zero_cols)){
        warning(sum(zero_cols), " cell(s) have zero or NA column sum in the transition matrix ",
                "(no neighbors after masking, or undefined correlation). ",
                "Their transition probabilities will be set to zero.")
        col_sums[zero_cols] <- 1  # avoid 0/0 or NA/0; column stays all-zero
    }

    tp <- t(t(tp) / col_sums)  # tp shows transition from a given column cell to different row cells
    tp <- as(tp,'dgCMatrix') # cast to sparse matrix

    # add rownames / colnames 
    rownames(tp) <- colnames(seurat_obj)
    colnames(tp) <- colnames(seurat_obj)

    # convert to Seurat Graph object, and add it to the Seurat object with the perturbation name
    tp <- Seurat::as.Graph(tp)
    graph_name <- paste0(perturbation_name, '_tp')
    seurat_obj@graphs[[graph_name]] <- tp

    # return the Seurat object 
    seurat_obj

}



################################################################################
# Velocyo C++ helper functions
# Is there a way that we could remove the dependency on Velocyto?
################################################################################

colDeltaCor_velocyto <- function(e, d, nthreads = 1L) {
    .Call('_velocyto_R_colDeltaCor', PACKAGE = 'velocyto.R', e, d, nthreads)
}

