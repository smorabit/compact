
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

    # check for selected Graph 
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

    # subset by selected features
    # and convert to dense matrix (TODO: figure out a way to avoid this?)
    exp_obs <- as.matrix(exp_obs[features,])
    exp_per <- as.matrix(exp_per[features,])

    delta <- exp_per - exp_obs

    # run the Velocyto colDeltaCor function
    if(use_velocyto){
        cc <- colDeltaCor_velocyto(exp_obs, delta, nthreads=n_threads)
    } else{
        cc <- SparseColDeltaCor(exp_obs, delta, cell_graph)
    }

    # fill the diagnoal with zeros (cells won't transition to self)
    diag(cc) <- 0

    # compute transition probs between cells
    tp <- exp(cc / corr_sigma) * cell_graph
    tp <- t(t(tp)/Matrix::colSums(tp)); # tp shows transition from a given column cell to different row cells
    tp <- as(tp,'dgCMatrix') # cast to sparse matrix

    # add rownames / colnames 
    rownames(tp) <- colnames(seurat_obj)
    colnames(tp) <- colnames(seurat_obj)

    # convert to Seurat Graph object, and add it to the Seurat object with the perturbation name
    tp <- Seurat::as.Graph(tp)
    graph_name <- paste0(perturbation_name, '_tp')
    seurat_obj@graphs[graph_name] <- tp

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

