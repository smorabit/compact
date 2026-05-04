
#' Predict Perturbation Pseudotime (Absorbing Markov Chain Hitting Time)
#'
#' @description
#' Calculates the expected hitting time for cells to reach a defined "sink" or target state 
#' under a simulated perturbation. Mathematically, this function models the perturbation 
#' transition probability matrix as an Absorbing Markov Chain. It calculates the expected 
#' number of random walk steps required for any transient cell to reach the absorbing 
#' sink states, effectively quantifying the transcriptomic "distance" or resistance to 
#' the perturbation-driven trajectory.
#'
#' @param seurat_obj A Seurat object containing the required graphs.
#' @param perturbation_name Character. The base name of the simulated perturbation. 
#'   The function expects to find a transition probability graph named 
#'   \code{paste0(perturbation_name, '_tp')} in \code{seurat_obj@graphs}.
#' @param graph Character. The name of the K-Nearest Neighbors (KNN) graph stored in 
#'   \code{seurat_obj@graphs} (e.g., "RNA_nn"). This is used to mask the transition 
#'   probabilities and force the random walk to strictly obey the biological manifold.
#' @param sink_cells Character vector. An explicit list of cell barcodes to be defined 
#'   as the target/absorbing states. If \code{NULL}, \code{group.by} and \code{group_name} 
#'   must be provided.
#' @param group.by Character. The name of a column in \code{seurat_obj@meta.data} containing 
#'   cluster or group identities. Used to automatically define sink cells.
#' @param group_name Character vector. The specific identity class(es) within the 
#'   \code{group.by} column to define as the sink/target states.
#' @param output_name Character. The column name used to store the resulting pseudotime 
#'   values in \code{seurat_obj@meta.data}. Default is \code{"perturbation_pseudotime"}.
#' @param max_iter Numeric. The maximum number of iterations for the fixed-point 
#'   iterative sparse solver. Default is \code{10000}.
#' @param tolerance Numeric. The convergence tolerance for the iterative solver. 
#'   Default is \code{1e-4}.
#' @param return_seurat Logical. Whether to return a Seurat object (default), or to return
#'   the calculated perturbation pseudotime as a numeric vector.
#'
#' @return A Seurat object with the calculated perturbation times appended to the 
#'   metadata under the column name specified by \code{output_name}. Sink cells 
#'   are strictly assigned a time of 0.
#'
#' @details 
#' This function uses a highly efficient, matrix-free iterative solver to calculate 
#' hitting times, ensuring it scales comfortably to hundreds of thousands of cells 
#' without dense matrix inversion. Cells that hit the \code{max_iter} ceiling are likely 
#' on disconnected manifold "islands" or face severe topological bottlenecks preventing 
#' them from reaching the specified sink states under the current perturbation.
#'
#' @export
PredictPerturbationTime <- function(
    seurat_obj,
    perturbation_name,
    graph,
    sink_cells = NULL,
    group.by = NULL,
    group_name = NULL,
    output_name = "perturbation_pseudotime",
    max_iter = 10000,
    tolerance = 1e-4,
    return_seurat = TRUE
){

    # -------------------------------------------------------------------------
    # input validation
    # -------------------------------------------------------------------------
    
    # check if the requested graphs exist in the Seurat object
    tp_graph_name <- paste0(perturbation_name, '_tp')
    
    if (!(graph %in% names(seurat_obj@graphs))) {
        stop(paste0("Error: KNN graph '", graph, "' not found in seurat_obj@graphs."))
    }
    if (!(tp_graph_name %in% names(seurat_obj@graphs))) {
        stop(paste0("Error: Transition graph '", tp_graph_name, "' not found. Did you run PerturbationTransitions?"))
    }

    # determine sink cells based on user input
    if (is.null(sink_cells)) {
        if (is.null(group.by) || is.null(group_name)) {
            stop("Error: You must provide either 'sink_cells' directly, OR both 'group.by' and 'group_name'.")
        }
        
        # Check if group.by exists in metadata
        if (!(group.by %in% colnames(seurat_obj@meta.data))) {
            stop(paste0("Error: Column '", group.by, "' not found in seurat_obj@meta.data."))
        }
        
        # Extract cell barcodes matching the group_name
        sink_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] %in% group_name]
        
        if (length(sink_cells) == 0) {
            stop(paste0("Error: No cells found matching '", paste(group_name, collapse=", "), "' in column '", group.by, "'."))
        }
        message(paste("Identified", length(sink_cells), "sink cells using", group.by, "=", paste(group_name, collapse=", ")))
    } else {
        # Validate provided sink_cells
        valid_sinks <- intersect(sink_cells, colnames(seurat_obj))
        if (length(valid_sinks) == 0) {
            stop("Error: None of the provided 'sink_cells' were found in the Seurat object.")
        }
        sink_cells <- valid_sinks
        message(paste("Using", length(sink_cells), "user-provided sink cells."))
    }

    # -------------------------------------------------------------------------
    # prepare matrices
    # -------------------------------------------------------------------------
    
    # get the KNN graph and set diagonal to 1
    cell_graph <- Graphs(seurat_obj, slot = graph)
    cell_graph <- Matrix::Matrix(cell_graph, sparse = TRUE)
    diag(cell_graph) <- 1

    # get the transition probability matrix
    # the stored TP is column-stochastic: tp[i,j] = P(j -> i).
    # transpose so that tp[i,j] = P(i -> j), then row-normalize to get
    # the correct row-stochastic transition matrix for downstream computations.
    tp <- Graphs(seurat_obj, slot = tp_graph_name)
    tp <- Matrix::Matrix(tp, sparse = TRUE)
    tp[is.na(tp)] <- 0
    tp <- t(tp)

    # mask transition probabilities with KNN graph
    tp <- tp * cell_graph

    # row-normalize to create Markov transition matrix P
    row_sums <- rowSums(tp)
    row_sums[row_sums == 0] <- 1
    P <- tp / row_sums

    # -------------------------------------------------------------------------
    # set up absorbing markov chain
    # -------------------------------------------------------------------------
    
    all_cells <- colnames(seurat_obj)
    transient_cells <- setdiff(all_cells, sink_cells)

    if (length(transient_cells) == 0) {
        stop("Error: All cells are defined as sink cells. No transient cells remaining to calculate time for.")
    }

    # extract the sub-matrix of transitions between transient states (Q matrix)
    Q <- P[transient_cells, transient_cells]

    # initialize starting times at 0
    t_current <- rep(0, length(transient_cells))
    ones <- rep(1, length(transient_cells))
    converged <- FALSE

    # -------------------------------------------------------------------------
    # iterative solver
    # -------------------------------------------------------------------------
    
    message("Starting iterative solver for hitting times...")

    for (i in 1:max_iter) {
        # core Markov update step: t_(n+1) = 1 + Q * t_n
        t_new <- ones + as.vector(Q %*% t_current)
        
        # check maximum difference for convergence
        max_diff <- max(abs(t_new - t_current))
        
        if (max_diff < tolerance) {
            message(paste("Converged successfully in", i, "iterations!"))
            t_current <- t_new
            converged <- TRUE
            break
        }
        
        t_current <- t_new
        
        if (i %% 100 == 0) {
            message(paste("Iteration:", i, "| Max change:", round(max_diff, 4)))
        }
    }

    if (!converged) {
        message("Warning: Reached max iterations. Some cells may be on disconnected 'islands' with no path to the sink.")
    }

    # -------------------------------------------------------------------------
    # format the output
    # -------------------------------------------------------------------------
    
    # map calculated times back to transient cells, and set sink cells strictly to 0
    perturbation_pseudotime <- data.frame(
        bc = c(transient_cells, sink_cells),
        time = c(t_current, rep(0, length(sink_cells)))
    )
    rownames(perturbation_pseudotime) <- perturbation_pseudotime$bc

    # add to Seurat metadata under the user-specified column name
    final_times <- perturbation_pseudotime[colnames(seurat_obj), 'time']
    names(final_times) <- colnames(seurat_obj)
    seurat_obj <- AddMetaData(
        object = seurat_obj,
        metadata = final_times,
        col.name = output_name
    )
    
    if(return_seurat){
        return(seurat_obj)
    } else{
        return(final_times)
    }
}

#' Predict Perturbation Attractor States
#'
#' @description
#' Performs unsupervised discovery of terminal attractor states (sinks) within a simulated 
#' perturbation vector field. Mathematically, this function calculates the stationary 
#' distribution of the directed Markov transition matrix by finding the dominant left 
#' eigenvector. Cells with the highest stationary probabilities represent deep phenotypic 
#' basins where the perturbation trajectory naturally pools .
#'
#' @param seurat_obj A Seurat object containing the required graphs.
#' @param perturbation_name Character. The base name of the simulated perturbation. 
#'   The function expects to find a transition probability graph named 
#'   \code{paste0(perturbation_name, '_tp')} in \code{seurat_obj@graphs}.
#' @param graph Character. The name of the K-Nearest Neighbors (KNN) graph stored in 
#'   \code{seurat_obj@graphs} (e.g., "RNA_nn") used to mask the transition probabilities.
#' @param output_name Character. The column name used to store the resulting attractor 
#'   scores in \code{seurat_obj@meta.data}. Default is \code{"attractor_score"}.
#' @param quantile_threshold Numeric. The percentile threshold used to explicitly define 
#'   which cells belong to the terminal sink state. Default is \code{0.98} (top 2%).
#' @param return_seurat Logical. If \code{TRUE}, returns the updated Seurat object. 
#'   If \code{FALSE}, returns a list containing the numeric attractor scores and a 
#'   character vector of the identified sink cells. Default is \code{TRUE}.
#'
#' @return If \code{return_seurat = TRUE}, a Seurat object with the calculated attractor 
#'   scores appended to the metadata. If \code{return_seurat = FALSE}, a list containing 
#'   \code{attractor_score} and \code{sink_cells}.
#'
#' @export
PredictAttractors <- function(
    seurat_obj,
    perturbation_name,
    graph,
    output_name = "attractor_score",
    quantile_threshold = 0.98,
    return_seurat = TRUE
){

    # -------------------------------------------------------------------------
    # input validation
    # -------------------------------------------------------------------------
    
    # check if the requested graphs exist in the Seurat object
    tp_graph_name <- paste0(perturbation_name, '_tp')
    
    if (!(graph %in% names(seurat_obj@graphs))) {
        stop(paste0("Error: KNN graph '", graph, "' not found in seurat_obj@graphs."))
    }
    if (!(tp_graph_name %in% names(seurat_obj@graphs))) {
        stop(paste0("Error: Transition graph '", tp_graph_name, "' not found. Did you run PerturbationTransitions?"))
    }

    # -------------------------------------------------------------------------
    # prepare matrices
    # -------------------------------------------------------------------------
    
    # get the KNN graph and set diagonal to 1
    cell_graph <- Graphs(seurat_obj, slot = graph)
    cell_graph <- Matrix::Matrix(cell_graph, sparse = TRUE)
    diag(cell_graph) <- 1

    # get the transition probability matrix
    # the stored TP is column-stochastic: tp[i,j] = P(j -> i).
    # transpose so that tp[i,j] = P(i -> j), then row-normalize to get
    # the correct row-stochastic transition matrix for downstream computations.
    tp <- Graphs(seurat_obj, slot = tp_graph_name)
    tp <- Matrix::Matrix(tp, sparse = TRUE)
    tp[is.na(tp)] <- 0
    tp <- t(tp)

    # mask transition probabilities with KNN graph
    tp <- tp * cell_graph

    # row-normalize to create Markov transition matrix P
    row_sums <- rowSums(tp)
    row_sums[row_sums == 0] <- 1
    P <- tp / row_sums

    # -------------------------------------------------------------------------
    # identify attractors
    # -------------------------------------------------------------------------

    # calculate the dominant left eigenvector of P (= dominant right eigenvector
    # of t(P)), which gives the stationary distribution of the random walk
    eig_res <- RSpectra::eigs(t(P), k = 1, which = "LR")

    # extract the real part; take absolute value before normalizing because
    # ARPACK's Lanczos solver can return sign-ambiguous eigenvectors and a
    # reducible P (disconnected graph components) can produce near-zero negative
    # entries that would otherwise yield negative attractor scores
    stat_dist <- abs(Re(eig_res$vectors[, 1]))
    stat_dist <- stat_dist / sum(stat_dist)
    names(stat_dist) <- rownames(P)

    # define the sinks based on the stationary distribution
    threshold <- quantile(stat_dist, quantile_threshold) 
    auto_sink_cells <- names(which(stat_dist >= threshold))

    # -------------------------------------------------------------------------
    # return results
    # -------------------------------------------------------------------------

    if(return_seurat){
        # add to the seurat object metadata
        seurat_obj <- AddMetaData(
            object = seurat_obj,
            metadata = stat_dist,
            col.name = output_name
        )
        return(seurat_obj)
        
    } else {
        # return a list containing both the continuous scores and the discrete sink cells
        return(list(
            attractor_score = stat_dist,
            sink_cells = auto_sink_cells
        ))
    }
}

#' Predict Perturbation Fates (Forward Diffusion)
#'
#' @description
#' Simulates the forward diffusion of a specific cell population under a perturbed 
#' vector field. This function models how a starting group of cells (sources) will 
#' transition over time, revealing their terminal fate bias. It computes the probability 
#' distribution of the population across the manifold after a specified number of 
#' Markov steps .
#'
#' @param seurat_obj A Seurat object containing the required graphs.
#' @param perturbation_name Character. The base name of the simulated perturbation. 
#'   The function expects to find a transition probability graph named 
#'   \code{paste0(perturbation_name, '_tp')} in \code{seurat_obj@graphs}.
#' @param graph Character. The name of the K-Nearest Neighbors (KNN) graph stored in 
#'   \code{seurat_obj@graphs} (e.g., "RNA_nn") used to mask the transition probabilities.
#' @param source_cells Character vector. An explicit list of cell barcodes to be defined 
#'   as the starting population. If \code{NULL}, \code{group.by} and \code{group_name} 
#'   must be provided.
#' @param group.by Character. The name of a column in \code{seurat_obj@meta.data} containing 
#'   cluster or group identities. Used to automatically define source cells.
#' @param group_name Character vector. The specific identity class(es) within the 
#'   \code{group.by} column to define as the starting population.
#' @param output_name Character. The column name used to store the resulting fate 
#'   probabilities in \code{seurat_obj@meta.data}. Default is \code{"forward_fate"}.
#' @param t_steps Numeric. The maximum number of forward diffusion steps to simulate. 
#'   Default is \code{100}.
#' @param tolerance Numeric. The convergence tolerance. If the probability distribution 
#'   stops changing between steps (reaches a stationary state early), the simulation 
#'   will stop. Default is \code{1e-5}.
#' @param rank_transform Logical. If \code{TRUE}, transforms the raw probabilities into 
#'   percentile ranks (ECDF) for better dynamic range during visualization. Default is \code{TRUE}.
#' @param return_seurat Logical. If \code{TRUE}, returns the updated Seurat object. 
#'   If \code{FALSE}, returns a list containing the probabilities and source cells. 
#'   Default is \code{TRUE}.
#'
#' @return If \code{return_seurat = TRUE}, a Seurat object with the calculated fate 
#'   scores appended to the metadata. If \code{return_seurat = FALSE}, a list containing 
#'   the fate probabilities and the initial source cells.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate where the "Pre-Exhausted" cells go after down-regulating the exhaustion module
#' seurat_obj <- PredictFates(
#'   seurat_obj = seurat_obj,
#'   perturbation_name = "Exhaustion_down",
#'   graph = "RNA_nn",
#'   group.by = "seurat_clusters",
#'   group_name = "7", # Assuming cluster 7 is the pre-exhausted state
#'   output_name = "fate_from_cluster7",
#'   t_steps = 100,
#'   rank_transform = TRUE
#' )
#' FeaturePlot(seurat_obj, features = "fate_from_cluster7")
#' }
PredictFates <- function(
    seurat_obj,
    perturbation_name,
    graph,
    source_cells = NULL,
    group.by = NULL,
    group_name = NULL,
    output_name = "forward_fate",
    t_steps = 100,
    tolerance = 1e-5,
    rank_transform = TRUE,
    return_seurat = TRUE
){

    # -------------------------------------------------------------------------
    # input validation
    # -------------------------------------------------------------------------
    
    # check if the requested graphs exist in the Seurat object
    tp_graph_name <- paste0(perturbation_name, '_tp')
    
    if (!(graph %in% names(seurat_obj@graphs))) {
        stop(paste0("Error: KNN graph '", graph, "' not found in seurat_obj@graphs."))
    }
    if (!(tp_graph_name %in% names(seurat_obj@graphs))) {
        stop(paste0("Error: Transition graph '", tp_graph_name, "' not found. Did you run PerturbationTransitions?"))
    }

    # determine source cells based on user input
    if (is.null(source_cells)) {
        if (is.null(group.by) || is.null(group_name)) {
            stop("Error: You must provide either 'source_cells' directly, OR both 'group.by' and 'group_name'.")
        }
        
        # check if group.by exists in metadata
        if (!(group.by %in% colnames(seurat_obj@meta.data))) {
            stop(paste0("Error: Column '", group.by, "' not found in seurat_obj@meta.data."))
        }
        
        # extract cell barcodes matching the group_name
        source_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] %in% group_name]
        
        if (length(source_cells) == 0) {
            stop(paste0("Error: No cells found matching '", paste(group_name, collapse=", "), "' in column '", group.by, "'."))
        }
        message(paste("Identified", length(source_cells), "source cells using", group.by, "=", paste(group_name, collapse=", ")))
    } else {
        # validate provided source_cells
        valid_sources <- intersect(source_cells, colnames(seurat_obj))
        if (length(valid_sources) == 0) {
            stop("Error: None of the provided 'source_cells' were found in the Seurat object.")
        }
        source_cells <- valid_sources
        message(paste("Using", length(source_cells), "user-provided source cells."))
    }

    # -------------------------------------------------------------------------
    # prepare matrices
    # -------------------------------------------------------------------------
    
    # get the KNN graph and set diagonal to 1
    cell_graph <- Graphs(seurat_obj, slot = graph)
    cell_graph <- Matrix::Matrix(cell_graph, sparse = TRUE)
    diag(cell_graph) <- 1

    # get the transition probability matrix
    # the stored TP is column-stochastic: tp[i,j] = P(j -> i).
    # transpose so that tp[i,j] = P(i -> j), then row-normalize to get
    # the correct row-stochastic transition matrix for downstream computations.
    tp <- Graphs(seurat_obj, slot = tp_graph_name)
    tp <- Matrix::Matrix(tp, sparse = TRUE)
    tp[is.na(tp)] <- 0
    tp <- t(tp)

    # mask transition probabilities with KNN graph
    tp <- tp * cell_graph

    # row-normalize to create Markov transition matrix P
    row_sums <- rowSums(tp)
    row_sums[row_sums == 0] <- 1
    P <- tp / row_sums

    # -------------------------------------------------------------------------
    # perform forward simulation
    # -------------------------------------------------------------------------

    # initialize the starting probability distribution (p_0)
    p_current <- rep(0, nrow(P))
    names(p_current) <- rownames(P)
    p_current[source_cells] <- 1 / length(source_cells)

    message(paste("Simulating forward diffusion for up to", t_steps, "steps..."))

    # diffuse the probability forward for t_steps with early stopping
    for (i in 1:t_steps) {
        p_new <- as.vector(p_current %*% P)
        
        # check for convergence (early stopping if stationary)
        if (max(abs(p_new - p_current)) < tolerance) {
            message(paste("Converged to a stationary distribution early at step", i))
            p_current <- p_new
            break
        }
        p_current <- p_new
    }

    # apply rank transformation if requested
    if (rank_transform) {
        final_scores <- ecdf(p_current)(p_current)
    } else {
        final_scores <- p_current
    }
    names(final_scores) <- rownames(P)

    # -------------------------------------------------------------------------
    # return results
    # -------------------------------------------------------------------------

    if (return_seurat) {
        # add to the seurat object metadata
        seurat_obj <- AddMetaData(
            object = seurat_obj,
            metadata = final_scores,
            col.name = output_name
        )
        return(seurat_obj)
        
    } else {
        # return a list containing the scores and the source cells
        return(list(
            fate_score = final_scores,
            source_cells = source_cells
        ))
    }
}



#' Predict Fate Commitment (Committor Probabilities)
#'
#' @description
#' Calculates the committor probability for every cell under a simulated perturbation. 
#' This function models the vector field as a Dirichlet boundary value problem to answer 
#' a specific question: "What is the exact probability that a cell will successfully 
#' reach the target state (Sink) before it falls back into the origin state (Source)?" 
#' This is highly effective for identifying epigenetic "points of no return" and 
#' transcriptomic bottlenecks .
#'
#' @param seurat_obj A Seurat object containing the required graphs.
#' @param perturbation_name Character. The base name of the simulated perturbation. 
#'   The function expects to find a transition probability graph named 
#'   \code{paste0(perturbation_name, '_tp')} in \code{seurat_obj@graphs}.
#' @param graph Character. The name of the K-Nearest Neighbors (KNN) graph stored in 
#'   \code{seurat_obj@graphs} (e.g., "RNA_nn") used to mask the transition probabilities.
#' @param source_cells Character vector. Cell barcodes to be defined as the starting 
#'   population (Boundary condition = 0).
#' @param sink_cells Character vector. Cell barcodes to be defined as the target 
#'   population (Boundary condition = 1).
#' @param group.by Character. The name of a column in \code{seurat_obj@meta.data}. 
#'   Used to automatically define source and sink cells if explicit lists are not provided.
#' @param source_group Character vector. The identity class(es) in \code{group.by} 
#'   representing the source state.
#' @param sink_group Character vector. The identity class(es) in \code{group.by} 
#'   representing the sink state.
#' @param output_name Character. The column name used to store the resulting probabilities 
#'   in \code{seurat_obj@meta.data}. Default is \code{"commitment_score"}.
#' @param max_iter Numeric. The maximum number of iterations for the solver. Default is \code{2000}.
#' @param tolerance Numeric. The convergence tolerance. Default is \code{1e-5}.
#' @param return_seurat Logical. If \code{TRUE}, returns the updated Seurat object. 
#'   If \code{FALSE}, returns a list containing the probabilities. Default is \code{TRUE}.
#'
#' @return If \code{return_seurat = TRUE}, a Seurat object with the calculated committor 
#'   probabilities appended to the metadata. If \code{return_seurat = FALSE}, a list 
#'   containing the probability scores.
#'
#' @export
PredictCommitment <- function(
    seurat_obj,
    perturbation_name,
    graph,
    source_cells = NULL,
    sink_cells = NULL,
    group.by = NULL,
    source_group = NULL,
    sink_group = NULL,
    output_name = "commitment_score",
    max_iter = 2000,
    tolerance = 1e-5,
    return_seurat = TRUE
){

    # -------------------------------------------------------------------------
    # input validation
    # -------------------------------------------------------------------------
    
    # check if the requested graphs exist in the Seurat object
    tp_graph_name <- paste0(perturbation_name, '_tp')
    
    if (!(graph %in% names(seurat_obj@graphs))) {
        stop(paste0("Error: KNN graph '", graph, "' not found in seurat_obj@graphs."))
    }
    if (!(tp_graph_name %in% names(seurat_obj@graphs))) {
        stop(paste0("Error: Transition graph '", tp_graph_name, "' not found. Did you run PerturbationTransitions?"))
    }

    # determine source cells
    if (is.null(source_cells)) {
        if (is.null(group.by) || is.null(source_group)) {
            stop("Error: You must provide either 'source_cells' directly, OR both 'group.by' and 'source_group'.")
        }
        source_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] %in% source_group]
        if (length(source_cells) == 0) stop("Error: No source cells found matching the provided group.")
    }

    # determine sink cells
    if (is.null(sink_cells)) {
        if (is.null(group.by) || is.null(sink_group)) {
            stop("Error: You must provide either 'sink_cells' directly, OR both 'group.by' and 'sink_group'.")
        }
        sink_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[group.by]] %in% sink_group]
        if (length(sink_cells) == 0) stop("Error: No sink cells found matching the provided group.")
    }

    # validate intersection
    if (length(intersect(source_cells, sink_cells)) > 0) {
        stop("Error: 'source_cells' and 'sink_cells' cannot overlap. A cell cannot be both the origin and the target.")
    }

    message(paste("Identified", length(source_cells), "source cells and", length(sink_cells), "sink cells."))

    # -------------------------------------------------------------------------
    # prepare matrices
    # -------------------------------------------------------------------------
    
    # get the KNN graph and set diagonal to 1
    cell_graph <- Graphs(seurat_obj, slot = graph)
    cell_graph <- Matrix::Matrix(cell_graph, sparse = TRUE)
    diag(cell_graph) <- 1

    # get the transition probability matrix
    # the stored TP is column-stochastic: tp[i,j] = P(j -> i).
    # transpose so that tp[i,j] = P(i -> j), then row-normalize to get
    # the correct row-stochastic transition matrix for downstream computations.
    tp <- Graphs(seurat_obj, slot = tp_graph_name)
    tp <- Matrix::Matrix(tp, sparse = TRUE)
    tp[is.na(tp)] <- 0
    tp <- t(tp)

    # mask transition probabilities with KNN graph
    tp <- tp * cell_graph

    # row-normalize to create Markov transition matrix P
    row_sums <- rowSums(tp)
    row_sums[row_sums == 0] <- 1
    P <- tp / row_sums

    # -------------------------------------------------------------------------
    # setup markov chain boundary value problem
    # -------------------------------------------------------------------------

    all_cells <- rownames(P) 
    transient_cells <- setdiff(all_cells, c(source_cells, sink_cells))
    
    message(paste("Calculating commitment for", length(transient_cells), "transient cells..."))

    # extract sub-matrices
    P_UU <- P[transient_cells, transient_cells]
    P_UB <- P[transient_cells, sink_cells]

    # calculate 'b' (probability of jumping directly to the sink in one step)
    b <- rowSums(P_UB)

    # -------------------------------------------------------------------------
    # solve for committor probabilities
    # -------------------------------------------------------------------------

    q_current <- rep(0, length(transient_cells))
    converged <- FALSE

    for (i in 1:max_iter) {
        # the core update step: q_new = P_UU * q_current + b
        q_new <- as.vector(P_UU %*% q_current) + b
        
        max_diff <- max(abs(q_new - q_current))
        
        if (max_diff < tolerance) {
            message(paste("Converged successfully in", i, "iterations!"))
            q_current <- q_new
            converged <- TRUE
            break
        }

        # optional: print progress every 100 iterations
        if (i %% 100 == 0) {
            message(paste("Iteration:", i, "| Max change:", round(max_diff, 4)))
        }
        
        q_current <- q_new
    }

    if (!converged) {
        message("Warning: Reached max iterations. Some cells may be mathematically disconnected.")
    }

    # -------------------------------------------------------------------------
    # return results
    # -------------------------------------------------------------------------

    # format the results back into a full-length vector
    committor_results <- data.frame(
        cell = all_cells,
        score = NA 
    )
    rownames(committor_results) <- committor_results$cell

    # assign the boundary conditions and the solved probabilities
    committor_results[sink_cells, "score"] <- 1
    committor_results[source_cells, "score"] <- 0
    committor_results[transient_cells, "score"] <- q_current

    final_scores <- committor_results[colnames(seurat_obj), "score"]
    names(final_scores) <- colnames(seurat_obj)

    if (return_seurat) {
        seurat_obj <- AddMetaData(
            object = seurat_obj,
            metadata = final_scores,
            col.name = output_name
        )
        return(seurat_obj)
        
    } else {
        return(final_scores)
    }
}