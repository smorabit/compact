#' FindConnectedComponents
#'
#' Inspects an SNN (or any named) graph from a Seurat object, or a raw
#' adjacency/transition matrix, and labels each cell by its connected
#' component. When run on a Seurat object, the result is written to a new
#' column in \code{seurat_obj@meta.data}. This is a useful diagnostic before
#' running \code{PerturbationTransitions}, because cells that belong to
#' different connected components cannot exchange transition probability mass
#' through the graph. It is also used internally by \code{PredictAttractors}
#' to check the post-masking transition matrix for disconnected components,
#' which can otherwise make the dominant eigenvalue non-unique and cause the
#' stationary-distribution eigensolver to fail to converge.
#'
#' @param seurat_obj A Seurat object containing at least one graph in
#'   \code{seurat_obj@graphs}. Ignored if \code{matrix} is provided.
#' @param graph Character. Name of the graph in \code{Graphs(seurat_obj)} to
#'   analyse. If \code{NULL} (default), the function auto-detects a graph whose
#'   name ends in \code{"_snn"}. When multiple SNN graphs exist the first one is
#'   used and a message is emitted. If no SNN graph is present the first
#'   available graph is used instead (with a warning). Ignored if
#'   \code{matrix} is provided.
#' @param matrix A square adjacency or transition matrix (sparse or dense) to
#'   analyse directly, bypassing the Seurat object entirely. Row/column names,
#'   if present, are used as cell identifiers. When supplied, the function
#'   returns a plain list of component info (see Value) instead of a Seurat
#'   object, and \code{seurat_obj}/\code{graph}/\code{meta_data_name} are
#'   ignored.
#' @param meta_data_name Character. Name of the new column written to
#'   \code{seurat_obj@meta.data}. Default \code{"connected_component"}.
#'   Ignored if \code{matrix} is provided.
#' @param verbose Logical. Print a summary of the component structure. Default
#'   \code{TRUE}.
#'
#' @return If \code{matrix} is \code{NULL} (default), the Seurat object with a
#'   new integer-factor column in \code{seurat_obj@meta.data}. Levels are
#'   ordered by component size (largest component = level 1) so the dominant
#'   component always has the lowest label. If \code{matrix} is supplied
#'   instead, a list with elements \code{membership} (a named integer vector,
#'   one entry per row/column of \code{matrix}, ranked the same way),
#'   \code{n_components}, and \code{component_sizes} (sizes in decreasing
#'   order).
#'
#' @details
#' The graph/matrix is treated as undirected for component finding — any
#' non-zero edge weight in either direction is interpreted as a connection,
#' and the diagonal is ignored. The underlying computation uses
#' \code{igraph::components()}, which implements a fast depth-first-search.
#'
#' @import Seurat
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom Matrix Matrix drop0
#' @export
FindConnectedComponents <- function(
    seurat_obj = NULL,
    graph = NULL,
    matrix = NULL,
    meta_data_name = "connected_component",
    verbose = TRUE
) {

    # -------------------------------------------------------------------------
    # matrix mode: analyse a raw adjacency/transition matrix directly
    # -------------------------------------------------------------------------

    if (!is.null(matrix)) {
        adj <- Matrix::Matrix(matrix, sparse = TRUE)
        comp <- .RankedComponents(adj)

        if (verbose) {
            .ReportComponents(comp, label = "provided matrix")
        }

        return(list(
            membership = comp$membership,
            n_components = comp$n_components,
            component_sizes = comp$component_sizes
        ))
    }

    # -------------------------------------------------------------------------
    # seurat object mode (original behavior, unchanged)
    # -------------------------------------------------------------------------

    # --- validate seurat object has graphs ---
    available_graphs <- names(seurat_obj@graphs)
    if (length(available_graphs) == 0) {
        stop("No graphs found in seurat_obj@graphs. Run FindNeighbors() first.")
    }

    # --- auto-detect or validate graph name ---
    if (is.null(graph)) {
        snn_graphs <- available_graphs[grepl("_snn$", available_graphs, ignore.case = TRUE)]
        if (length(snn_graphs) >= 1) {
            graph <- snn_graphs[1]
            if (length(snn_graphs) > 1 && verbose) {
                message("Multiple SNN graphs found: ",
                        paste(snn_graphs, collapse = ", "),
                        ". Using '", graph, "'.")
            }
        } else {
            graph <- available_graphs[1]
            warning("No graph ending in '_snn' found. Falling back to '", graph,
                    "'. Pass graph = '<name>' to specify explicitly.")
        }
        if (verbose) {
            message("Using graph: '", graph, "'")
        }
    } else {
        if (!(graph %in% available_graphs)) {
            stop("Graph '", graph, "' not found in seurat_obj@graphs. ",
                 "Available graphs: ", paste(available_graphs, collapse = ", "))
        }
    }

    # --- extract adjacency matrix ---
    adj <- seurat_obj@graphs[[graph]]
    adj <- Matrix::Matrix(adj, sparse = TRUE)

    comp <- .RankedComponents(adj)

    # --- write to meta.data (preserving cell order) ---
    seurat_obj@meta.data[[meta_data_name]] <- factor(
        comp$membership,
        levels = seq_len(comp$n_components)
    )

    if (verbose) {
        .ReportComponents(comp, label = paste0("graph '", graph, "'"))
        message("Component labels written to seurat_obj@meta.data[['",
                meta_data_name, "']].")
    }

    seurat_obj
}

# computes connected components of a square adjacency/weight matrix, treated
# as undirected with the diagonal ignored, and ranks component labels by size
# (largest = 1). returns a named-membership vector aligned to rownames(adj).
.RankedComponents <- function(adj) {

    # remove self-loops before component analysis
    diag(adj) <- 0

    g <- igraph::graph_from_adjacency_matrix(
        adj,
        mode     = "undirected",
        weighted = TRUE,
        diag     = FALSE
    )
    comp <- igraph::components(g)

    # order levels by component size (largest = 1)
    size_order   <- order(comp$csize, decreasing = TRUE)
    rank_map     <- integer(comp$no)
    rank_map[size_order] <- seq_len(comp$no)
    ranked_membership <- rank_map[comp$membership]
    names(ranked_membership) <- rownames(adj)

    list(
        membership = ranked_membership,
        n_components = comp$no,
        component_sizes = sort(comp$csize, decreasing = TRUE)
    )
}

# prints the same verbose component summary used by both the seurat-object
# and matrix-input modes of FindConnectedComponents
.ReportComponents <- function(comp, label) {

    n_cells   <- length(comp$membership)
    n_comp    <- comp$n_components
    main_size <- comp$component_sizes[1]
    main_pct  <- round(100 * main_size / n_cells, 1)

    message(sprintf(
        "Found %d connected component%s across %d cells.",
        n_comp, if (n_comp == 1) "" else "s", n_cells
    ))
    message(sprintf(
        "Largest component: %d cells (%.1f%% of total).",
        main_size, main_pct
    ))

    if (n_comp > 1) {
        small_sizes <- comp$component_sizes[-1]
        message(sprintf(
            "Remaining %d component%s: size%s %s.",
            n_comp - 1,
            if (n_comp - 1 == 1) "" else "s",
            if (n_comp - 1 == 1) "" else "s",
            paste(small_sizes, collapse = ", ")
        ))
        warning(
            n_comp, " connected components detected in ", label, ". ",
            "Cells in different components cannot exchange transition probability ",
            "mass via the graph. Consider increasing FindNeighbors() k.param or ",
            "prune.SNN to reduce the number of isolated components.",
            call. = FALSE
        )
    }
}
