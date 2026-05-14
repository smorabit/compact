#' FindConnectedComponents
#'
#' Inspects an SNN (or any named) graph from a Seurat object and labels each
#' cell by its connected component. The result is written to a new column in
#' \code{seurat_obj@meta.data}. This is a useful diagnostic before running
#' \code{PerturbationTransitions}, because cells that belong to different
#' connected components cannot exchange transition probability mass through the
#' graph.
#'
#' @param seurat_obj A Seurat object containing at least one graph in
#'   \code{seurat_obj@graphs}.
#' @param graph Character. Name of the graph in \code{Graphs(seurat_obj)} to
#'   analyse. If \code{NULL} (default), the function auto-detects a graph whose
#'   name ends in \code{"_snn"}. When multiple SNN graphs exist the first one is
#'   used and a message is emitted. If no SNN graph is present the first
#'   available graph is used instead (with a warning).
#' @param meta_data_name Character. Name of the new column written to
#'   \code{seurat_obj@meta.data}. Default \code{"connected_component"}.
#' @param verbose Logical. Print a summary of the component structure. Default
#'   \code{TRUE}.
#'
#' @return The Seurat object with a new integer-factor column in
#'   \code{seurat_obj@meta.data}. Levels are ordered by component size
#'   (largest component = level 1) so the dominant component always has the
#'   lowest label.
#'
#' @details
#' The graph is treated as undirected for component finding — any non-zero edge
#' weight is interpreted as a connection. The underlying computation uses
#' \code{igraph::components()}, which implements a fast depth-first-search.
#'
#' @import Seurat
#' @importFrom igraph graph_from_adjacency_matrix components
#' @importFrom Matrix Matrix drop0
#' @export
FindConnectedComponents <- function(
    seurat_obj,
    graph = NULL,
    meta_data_name = "connected_component",
    verbose = TRUE
) {

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

    # remove self-loops before component analysis
    diag(adj) <- 0

    # --- build undirected igraph and find components ---
    g <- igraph::graph_from_adjacency_matrix(
        adj,
        mode     = "undirected",
        weighted = TRUE,
        diag     = FALSE
    )
    comp <- igraph::components(g)

    # --- order levels by component size (largest = 1) ---
    size_order   <- order(comp$csize, decreasing = TRUE)
    rank_map     <- integer(comp$no)
    rank_map[size_order] <- seq_len(comp$no)
    ranked_membership <- rank_map[comp$membership]

    # store as an ordered factor
    component_labels <- factor(
        ranked_membership,
        levels = seq_len(comp$no)
    )

    # --- write to meta.data (preserving cell order) ---
    seurat_obj@meta.data[[meta_data_name]] <- component_labels

    # --- verbose summary ---
    if (verbose) {
        n_cells <- length(comp$membership)
        n_comp  <- comp$no
        main_size <- comp$csize[size_order[1]]
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
            small_sizes <- sort(comp$csize[size_order[-1]], decreasing = TRUE)
            message(sprintf(
                "Remaining %d component%s: size%s %s.",
                n_comp - 1,
                if (n_comp - 1 == 1) "" else "s",
                if (n_comp - 1 == 1) "" else "s",
                paste(small_sizes, collapse = ", ")
            ))
            warning(
                n_comp, " connected components detected in graph '", graph, "'. ",
                "Cells in different components cannot exchange transition probability ",
                "mass via the SNN. Consider increasing FindNeighbors() k.param or ",
                "prune.SNN to reduce the number of isolated components.",
                call. = FALSE
            )
        }

        message("Component labels written to seurat_obj@meta.data[['",
                meta_data_name, "']].")
    }

    seurat_obj
}
