
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
        if(use_graph_tp){
            cc <- colDeltaCor_knn(exp_obs, delta, cell_graph, n_threads)
        } else{
            cc <- colDeltaCor_knn(exp_obs, delta, NULL, n_threads)
        }
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



#' colDeltaCor_all
#'
#' This function calculates cell-to-cell correlations based on observed and 
#' perturbed gene expression matrices, with an option to limit computations 
#' to pairs of cells connected in a KNN graph. Correlations are computed 
#' in parallel to optimize performance.
#'
#' @return A matrix of cell-to-cell correlation values.
#'
#' @param observed_matrix A matrix of observed gene expression values (genes by cells).
#' @param delta_matrix A matrix of perturbed gene expression values (genes by cells).
#' @param knn_graph Optional; an adjacency matrix representing the k-nearest-neighbor (KNN) graph. 
#' If provided, only correlations for pairs of cells connected in the KNN graph will be computed.
#' @param n_cores The number of threads to use for parallel processing. Defaults to `detectCores() - 1`.
#'
#' @details
#' `colDeltaCor_all` calculates cell-to-cell correlations based on variance-stabilized 
#' transformed differences between observed and perturbed gene expression matrices.
#' The function computes correlations either for all cell pairs or, when a KNN graph is 
#' provided, only for cells connected in the graph.
#'
#' The primary steps of this analysis are:
#'
#' 1. Variance-stabilizing transformation: For each cell, observed and perturbed expression 
#'    values are transformed to reduce variance effects.
#'
#' 2. Cell-to-cell correlation computation: For each cell, the function calculates correlation 
#'    values with either all other cells or only the nearest neighbors, as defined by the KNN graph.
#'
#' 3. Parallel processing: Correlations are computed in parallel for optimized performance.
#'
#' @import parallel pbapply
#' @export
colDeltaCor_knn <- function(observed_matrix, delta_matrix, knn_graph = NULL, n_cores = detectCores() - 1) {
  
  n_cells <- ncol(observed_matrix)
  n_genes <- nrow(observed_matrix)
  
  # Variance-stabilizing transformation applied across all cells at once
  q_delta <- sign(delta_matrix) * sqrt(abs(delta_matrix))
  q_observed <- sign(observed_matrix) * sqrt(abs(observed_matrix))
  
  # Define function to compute correlations for one cell
  compute_correlations_for_cell <- function(i, q_delta, q_observed, n_cells, knn_neighbors = NULL) {

    # Initialize a vector to store correlations for cell i
    correlations <- numeric(n_cells)
    
    # Determine which cells to calculate correlations for
    if (is.null(knn_neighbors)) {
      neighbors <- 1:n_cells  # All cells
    } else {
      neighbors <- knn_neighbors[[i]]  # Only cells linked by KNN
    }
    
    for (j in neighbors) {
      if (j == i) next  # Skip self
      
      # Difference between observed expressions for variance-stabilizing
      q_obs_diff <- q_observed[, j] - q_observed[, i]
      
      # Dot product and magnitudes for correlation
      dot_product <- sum(q_delta[, i] * q_obs_diff)
      magnitude_delta <- sqrt(sum(q_delta[, i]^2))
      magnitude_obs_diff <- sqrt(sum(q_obs_diff^2))
      
      # Calculate correlation, handling cases where magnitudes are zero
      if (magnitude_delta == 0 || magnitude_obs_diff == 0) {
        correlations[j] <- 0
      } else {
        correlations[j] <- dot_product / (magnitude_delta * magnitude_obs_diff)
      }
    }
    
    return(correlations)
  }
  
  # Preprocess the KNN neighbors list, if a graph is provided
  knn_neighbors <- NULL
  if (!is.null(knn_graph)) {
    knn_neighbors <- lapply(1:n_cells, function(i) which(knn_graph[i, ] != 0))
  }
  
  # Parallel processing for each cell with a progress bar
  cl <- parallel::makeCluster(n_cores)
  delta_correlation_list <- pbapply::pblapply(
    1:n_cells, 
    compute_correlations_for_cell, 
    cl = cl, 
    q_delta = q_delta, 
    q_observed = q_observed, 
    n_cells = n_cells,
    knn_neighbors = knn_neighbors
  )
  parallel::stopCluster(cl)
  
  # Combine the list of correlations into a matrix
  delta_correlation <- do.call(rbind, delta_correlation_list)
  delta_correlation <- t(delta_correlation)
  
  return(delta_correlation)
}


# this one doesn't include the KNN option
colDeltaCor_all <- function(observed_matrix, delta_matrix, n_cores = detectCores() - 1) {
  
  n_cells <- ncol(observed_matrix)
  n_genes <- nrow(observed_matrix)
  
  # Variance-stabilizing transformation applied across all cells at once
  q_delta <- sign(delta_matrix) * sqrt(abs(delta_matrix))
  q_observed <- sign(observed_matrix) * sqrt(abs(observed_matrix))
  
  # Define function to compute correlations for one cell
  compute_correlations_for_cell <- function(i, q_delta, q_observed, n_cells) {
    # Initialize a vector to store correlations for cell i
    correlations <- numeric(n_cells)
    
    for (j in 1:n_cells) {
      if (j == i) next  # Skip self
      
      # Difference between observed expressions for variance-stabilizing
      q_obs_diff <- q_observed[, j] - q_observed[, i]
      
      # Dot product and magnitudes for correlation
      dot_product <- sum(q_delta[, i] * q_obs_diff)
      magnitude_delta <- sqrt(sum(q_delta[, i]^2))
      magnitude_obs_diff <- sqrt(sum(q_obs_diff^2))
      
      # Calculate correlation, handling cases where magnitudes are zero
      if (magnitude_delta == 0 || magnitude_obs_diff == 0) {
        correlations[j] <- 0
      } else {
        correlations[j] <- dot_product / (magnitude_delta * magnitude_obs_diff)
      }
    }
    
    return(correlations)
  }
  
  # Parallel processing for each cell with a progress bar
  # NOTE: 
  cl <- makeCluster(n_cores)
  delta_correlation_list <- pblapply(
    1:n_cells, 
    compute_correlations_for_cell, 
    cl = cl, 
    q_delta = q_delta, 
    q_observed = q_observed, 
    n_cells = n_cells
  )
  stopCluster(cl)
  
  # Combine the list of correlations into a matrix
  delta_correlation <- do.call(rbind, delta_correlation_list)
  delta_correlation <- t(delta_correlation)

  return(delta_correlation)
}


################################################################################
# Velocyo C++ helper functions
# Is there a way that we could remove the dependency on Velocyto?
################################################################################

colDeltaCor_velocyto <- function(e, d, nthreads = 1L) {
    .Call('_velocyto_R_colDeltaCor', PACKAGE = 'velocyto.R', e, d, nthreads)
}

