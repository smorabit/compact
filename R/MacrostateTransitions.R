#' MacrostateTransitions
#'
#' @description
#' Coarse-grains the cell-level perturbation transition probability (TP) matrix
#' into a group-level summary matrix (K × K), analogous to the macrostate
#' transition probability heatmap in CellRank (Lange et al. 2022, Nature Methods,
#' Fig. 2c). For each pair of user-defined cell groups, computes the average
#' transition probability flowing from cells in the source group to cells in the
#' destination group under the specified perturbation.
#'
#' The diagonal of the resulting matrix gives a \emph{stability index} per group:
#' the mean probability that a cell in that group transitions to another cell
#' \emph{within the same group}. High stability (close to 1) indicates
#' attractor-like behavior — the perturbation does not move cells out of this
#' group. Low stability indicates a transient source state from which cells are
#' being redirected toward other groups.
#'
#' Unlike \code{PlotTransitionVectors}, this analysis is entirely independent of
#' the low-dimensional embedding. All computations operate on the full-dimensional
#' transition matrix, avoiding the distortions introduced by UMAP/t-SNE
#' projections.
#'
#' @param seurat_obj A Seurat object with a computed perturbation transition
#'   probability matrix stored in \code{seurat_obj@graphs} (produced by
#'   \code{ModulePerturbation}, \code{TFPerturbation}, or
#'   \code{PerturbationTransitions}).
#' @param perturbation_name Character. The base name of the perturbation. The
#'   function looks up the TP graph as \code{paste0(perturbation_name, '_tp')}.
#' @param graph Character. Name of the KNN graph stored in
#'   \code{seurat_obj@graphs} (e.g., \code{"RNA_nn"}). Used to mask transition
#'   probabilities to the biological manifold, identical to the masking applied
#'   in \code{PredictAttractors} and \code{PredictPerturbationTime}.
#' @param group.by Character. Metadata column containing group labels used for
#'   coarse-graining (e.g., \code{"functional_anno"}, \code{"seurat_clusters"}).
#' @param normalize_rows Logical. If \code{TRUE} (default), each row of the
#'   coarse-grained matrix Q is normalized to sum to 1, making it row-stochastic.
#'   If \code{FALSE}, rows sum to the fraction of total outflow captured within
#'   the non-dropped groups (useful for diagnosing how much flow escapes excluded
#'   small groups).
#' @param min_group_size Integer. Groups with fewer than this many cells are
#'   excluded before computing Q, with a warning. Default: \code{5}.
#' @param store_result Logical. If \code{TRUE} (default), the result list is
#'   stored in \code{seurat_obj@misc$MacrostateTransitions[[result_name]]}.
#' @param result_name Character. Key used for \code{@misc} storage. Defaults to
#'   \code{perturbation_name}.
#' @param verbose Logical. Print progress messages. Default: \code{TRUE}.
#'
#' @return The Seurat object (if \code{store_result = TRUE}) with the following
#'   list stored in \code{seurat_obj@misc$MacrostateTransitions[[result_name]]}:
#'   \describe{
#'     \item{\code{Q}}{K × K numeric matrix. Coarse-grained (row-stochastic)
#'       transition probability matrix. \code{Q[k1, k2]} is the average
#'       probability of a cell in group \code{k1} transitioning to a cell in
#'       group \code{k2}.}
#'     \item{\code{stability}}{Named numeric vector of length K. Diagonal of Q;
#'       the self-transition probability (stability index) for each group.}
#'     \item{\code{group_sizes}}{Named integer vector. Number of cells per group
#'       after small-group filtering.}
#'     \item{\code{group.by}}{Character. The metadata column used.}
#'     \item{\code{perturbation_name}}{Character. The perturbation name.}
#'   }
#'   If \code{store_result = FALSE}, returns the result list directly.
#'
#' @details
#' \strong{Mathematical summary:}
#'
#' Let P be the n × n row-stochastic cell-level transition matrix (after
#' transposing the stored column-stochastic TP, applying KNN masking, and
#' row-normalizing). Let M be the n × K sparse indicator matrix where
#' \code{M[i, k] = 1} if cell i belongs to group k. Then:
#'
#' \deqn{Q = \mathrm{diag}(1/n_k) \cdot M^\top \cdot P \cdot M}
#'
#' where \eqn{n_k} is the number of cells in group k. Q is row-stochastic
#' (each row sums to 1) because P is row-stochastic and M sums each row over
#' all K groups.
#'
#' \strong{Complexity:} O(n · k_nn · K). The bottleneck is \code{P \%*\% M}
#' (sparse n × n times sparse n × K). For n = 150,000, k_nn = 20, K = 20 this
#' requires ~3M floating-point operations and ~24 MB of sparse storage. No dense
#' n × n matrix is ever materialized.
#'
#' @seealso \code{\link{PlotMacrostateTransitions}},
#'   \code{\link{PredictAttractors}}, \code{\link{PredictFates}},
#'   \code{\link{VectorFieldCoherence}}
#'
#' @references
#' Lange, M. et al. CellRank for directed single-cell fate mapping.
#' \emph{Nature Methods} \strong{19}, 159–170 (2022).
#' \doi{10.1038/s41592-021-01346-6}
#'
#' @export
MacrostateTransitions <- function(
    seurat_obj,
    perturbation_name,
    graph,
    group.by,
    normalize_rows  = TRUE,
    min_group_size  = 5,
    store_result    = TRUE,
    result_name     = NULL,
    verbose         = TRUE
) {

    # -------------------------------------------------------------------------
    # input validation
    # -------------------------------------------------------------------------

    tp_graph_name <- paste0(perturbation_name, '_tp')

    if (!(graph %in% names(seurat_obj@graphs))) {
        stop(paste0(
            "Error: KNN graph '", graph,
            "' not found in seurat_obj@graphs. ",
            "Available graphs: ", paste(names(seurat_obj@graphs), collapse = ", ")
        ))
    }

    if (!(tp_graph_name %in% names(seurat_obj@graphs))) {
        stop(paste0(
            "Error: Transition probability graph '", tp_graph_name,
            "' not found in seurat_obj@graphs. ",
            "Did you run ModulePerturbation() or PerturbationTransitions()?",
            " Available graphs: ", paste(names(seurat_obj@graphs), collapse = ", ")
        ))
    }

    if (!(group.by %in% colnames(seurat_obj@meta.data))) {
        stop(paste0(
            "Error: Metadata column '", group.by,
            "' not found in seurat_obj@meta.data. ",
            "Available columns: ",
            paste(colnames(seurat_obj@meta.data), collapse = ", ")
        ))
    }

    if (is.null(result_name)) {
        result_name <- perturbation_name
    }

    # -------------------------------------------------------------------------
    # prepare group labels and filter small groups
    # -------------------------------------------------------------------------

    group_labels <- seurat_obj@meta.data[[group.by]]

    # convert to character to avoid factor level issues
    group_labels <- as.character(group_labels)

    # compute group sizes
    group_tab <- table(group_labels)
    all_groups <- names(group_tab)

    # identify and warn about small groups
    small_groups <- names(group_tab[group_tab < min_group_size])
    if (length(small_groups) > 0) {
        warning(paste0(
            "The following groups have fewer than ", min_group_size,
            " cells and will be excluded: ",
            paste(small_groups, collapse = ", ")
        ))
    }

    keep_groups <- names(group_tab[group_tab >= min_group_size])
    if (length(keep_groups) == 0) {
        stop(paste0(
            "Error: No groups with >= ", min_group_size,
            " cells remain after filtering. ",
            "Reduce min_group_size or check your group.by column."
        ))
    }

    # filter cells to kept groups
    keep_cells <- which(group_labels %in% keep_groups)
    group_labels_keep <- group_labels[keep_cells]
    group_levels <- sort(unique(group_labels_keep))
    K <- length(group_levels)

    if (verbose) {
        message(paste0(
            "MacrostateTransitions: ", K, " groups retained from '", group.by,
            "' (", length(keep_cells), " of ", ncol(seurat_obj), " cells)."
        ))
    }

    # -------------------------------------------------------------------------
    # prepare the row-stochastic transition matrix P
    # -------------------------------------------------------------------------

    if (verbose) message("Loading and preparing transition probability matrix...")

    # load KNN graph; set self-loops to 1 (consistent with PredictAttractors)
    cell_graph <- SeuratObject::Graphs(seurat_obj, slot = graph)
    cell_graph <- Matrix::Matrix(cell_graph, sparse = TRUE)
    diag(cell_graph) <- 1

    # load the stored TP (column-stochastic: tp[i,j] = P(j -> i))
    tp <- SeuratObject::Graphs(seurat_obj, slot = tp_graph_name)
    tp <- Matrix::Matrix(tp, sparse = TRUE)
    tp[is.na(tp)] <- 0

    # transpose to row-stochastic: P[i,j] = P(i -> j)
    tp <- Matrix::t(tp)

    # mask to KNN neighborhood
    tp <- tp * cell_graph

    # row-normalize
    row_sums <- Matrix::rowSums(tp)
    row_sums[row_sums == 0] <- 1
    P <- tp / row_sums

    # subset to cells in kept groups (rows and columns)
    P <- P[keep_cells, keep_cells]

    # re-normalize rows after subsetting (some outflow to excluded cells is lost)
    row_sums2 <- Matrix::rowSums(P)
    row_sums2[row_sums2 == 0] <- 1
    P <- P / row_sums2

    # -------------------------------------------------------------------------
    # build sparse membership matrix M  (n_keep × K)
    # -------------------------------------------------------------------------

    if (verbose) message("Building membership matrix and computing coarse-grained Q...")

    group_idx <- as.integer(factor(group_labels_keep, levels = group_levels))
    n_keep    <- length(keep_cells)

    M <- Matrix::sparseMatrix(
        i    = seq_len(n_keep),
        j    = group_idx,
        x    = 1,
        dims = c(n_keep, K),
        dimnames = list(NULL, group_levels)
    )

    # -------------------------------------------------------------------------
    # compute Q = diag(1/n_k) * M^T * P * M
    # -------------------------------------------------------------------------

    group_sizes <- Matrix::colSums(M)          # K-vector of group sizes
    M_norm <- Matrix::Diagonal(x = 1 / group_sizes) %*% Matrix::t(M)  # K × n

    # two sparse multiplications: P %*% M is the expensive step (n × K)
    # then M_norm %*% (P %*% M) is K × n times n × K = K × K
    PM <- P %*% M                              # n × K sparse
    Q  <- as.matrix(M_norm %*% PM)            # K × K dense (tiny)

    rownames(Q) <- group_levels
    colnames(Q) <- group_levels

    # optional row re-normalization after subsetting may have left rows not
    # summing exactly to 1; apply normalize_rows if requested
    if (normalize_rows) {
        Q_row_sums <- rowSums(Q)
        Q_row_sums[Q_row_sums == 0] <- 1
        Q <- Q / Q_row_sums
    }

    # -------------------------------------------------------------------------
    # compute stability index and package result
    # -------------------------------------------------------------------------

    stability    <- diag(Q)
    names(stability) <- group_levels

    group_sizes_int          <- as.integer(group_sizes)
    names(group_sizes_int)   <- group_levels

    result <- list(
        Q                = Q,
        stability        = stability,
        group_sizes      = group_sizes_int,
        group.by         = group.by,
        perturbation_name = perturbation_name
    )

    if (verbose) {
        message("Done. Stability index summary:")
        message(paste(
            sprintf("  %-30s %.3f  (n = %d)", group_levels, stability, group_sizes_int),
            collapse = "\n"
        ))
    }

    # -------------------------------------------------------------------------
    # store result and return
    # -------------------------------------------------------------------------

    if (store_result) {
        if (is.null(seurat_obj@misc$MacrostateTransitions)) {
            seurat_obj@misc$MacrostateTransitions <- list()
        }
        seurat_obj@misc$MacrostateTransitions[[result_name]] <- result
        return(seurat_obj)
    } else {
        return(result)
    }
}


#' PlotMacrostateTransitions
#'
#' @description
#' Visualizes the coarse-grained group-level transition probability matrix
#' computed by \code{\link{MacrostateTransitions}} as a heatmap, analogous to
#' the macrostate transition probability heatmap in CellRank (Lange et al. 2022,
#' Fig. 2c).
#'
#' Each tile \code{[source, dest]} is colored by the average transition
#' probability flowing from the source group to the destination group under the
#' specified perturbation. The diagonal tiles represent self-transitions
#' (stability index) and are optionally annotated with the numeric stability
#' value and a heavier border for visual emphasis.
#'
#' @param seurat_obj A Seurat object with results stored in
#'   \code{seurat_obj@misc$MacrostateTransitions} (produced by
#'   \code{\link{MacrostateTransitions}}).
#' @param perturbation_name Character. Name of the perturbation to visualize.
#'   Used to look up \code{result_name} in \code{@misc$MacrostateTransitions}.
#' @param result_name Character. Key to retrieve from
#'   \code{seurat_obj@misc$MacrostateTransitions}. Defaults to
#'   \code{perturbation_name}.
#' @param group_order Character vector. Order of groups on both axes. Defaults
#'   to alphabetical order (the order stored in the result). Both the x-axis
#'   (destination) and the y-axis (source, reversed) follow this order.
#' @param color_scale Character vector. Two or more colors defining the
#'   low-to-high color gradient, passed to
#'   \code{ggplot2::scale_fill_gradientn}. Default:
#'   \code{c("white", "#2166AC")} (white to dark blue). For a diverging scale
#'   (e.g., when plotting ΔQ between perturbations), supply three colors with
#'   the midpoint as the second element.
#' @param show_stability Logical. If \code{TRUE} (default), overlay the
#'   numeric stability index on each diagonal tile.
#' @param stability_size Numeric. Font size for the stability text annotation.
#'   Default: \code{3}.
#' @param diagonal_border Logical. If \code{TRUE} (default), draw a heavier
#'   black border around the diagonal tiles to visually highlight
#'   self-transitions.
#' @param text_size Numeric. Base text size for axis labels and theme. Default:
#'   \code{10}.
#' @param legend_title Character. Title for the fill legend. Default:
#'   \code{"Transition\nProbability"}.
#' @param title Character. Plot title. Defaults to \code{perturbation_name}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{MacrostateTransitions}}, \code{\link{PredictAttractors}},
#'   \code{\link{VectorFieldCoherence}}
#'
#' @references
#' Lange, M. et al. CellRank for directed single-cell fate mapping.
#' \emph{Nature Methods} \strong{19}, 159–170 (2022).
#' \doi{10.1038/s41592-021-01346-6}
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradientn
#'   scale_x_discrete scale_y_discrete labs theme_bw theme element_text
#'   element_blank element_line element_rect
#'
#' @export
PlotMacrostateTransitions <- function(
    seurat_obj,
    perturbation_name,
    result_name      = NULL,
    group_order      = NULL,
    color_scale      = c("white", "#2166AC"),
    show_stability   = TRUE,
    stability_size   = 3,
    diagonal_border  = TRUE,
    text_size        = 10,
    legend_title     = "Transition\nProbability",
    title            = NULL
) {

    # -------------------------------------------------------------------------
    # retrieve stored result
    # -------------------------------------------------------------------------

    if (is.null(result_name)) {
        result_name <- perturbation_name
    }

    if (is.null(seurat_obj@misc$MacrostateTransitions)) {
        stop(paste0(
            "Error: No MacrostateTransitions results found in seurat_obj@misc. ",
            "Run MacrostateTransitions() first."
        ))
    }

    if (!(result_name %in% names(seurat_obj@misc$MacrostateTransitions))) {
        stop(paste0(
            "Error: Result '", result_name, "' not found in ",
            "seurat_obj@misc$MacrostateTransitions. ",
            "Available results: ",
            paste(names(seurat_obj@misc$MacrostateTransitions), collapse = ", ")
        ))
    }

    res       <- seurat_obj@misc$MacrostateTransitions[[result_name]]
    Q         <- res$Q
    stability <- res$stability

    # -------------------------------------------------------------------------
    # group ordering
    # -------------------------------------------------------------------------

    all_groups <- rownames(Q)

    if (is.null(group_order)) {
        group_order <- all_groups
    } else {
        missing_groups <- setdiff(group_order, all_groups)
        if (length(missing_groups) > 0) {
            warning(paste0(
                "The following groups in 'group_order' are not in the result and will be ignored: ",
                paste(missing_groups, collapse = ", ")
            ))
        }
        group_order <- intersect(group_order, all_groups)
        if (length(group_order) == 0) {
            stop("Error: No valid groups remain after filtering 'group_order'. Check group names.")
        }
        # subset Q and stability to the requested groups (in case group_order is a subset)
        Q         <- Q[group_order, group_order, drop = FALSE]
        stability <- stability[group_order]
    }

    # -------------------------------------------------------------------------
    # reshape Q to long format for ggplot2
    # -------------------------------------------------------------------------

    df_long <- expand.grid(
        source = group_order,
        dest   = group_order,
        stringsAsFactors = FALSE
    )
    df_long$probability <- Q[cbind(df_long$source, df_long$dest)]

    # enforce factor ordering
    df_long$source <- factor(df_long$source, levels = group_order)
    df_long$dest   <- factor(df_long$dest,   levels = group_order)

    # diagonal data frame for annotation
    diag_df <- data.frame(
        source = factor(group_order, levels = group_order),
        dest   = factor(group_order, levels = group_order),
        label  = sprintf("%.2f", stability[group_order]),
        stringsAsFactors = FALSE
    )

    # -------------------------------------------------------------------------
    # plot title
    # -------------------------------------------------------------------------

    if (is.null(title)) {
        title <- perturbation_name
    }

    # -------------------------------------------------------------------------
    # build heatmap
    # -------------------------------------------------------------------------

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = dest, y = source, fill = probability)) +
        ggplot2::geom_tile(color = "grey85", linewidth = 0.3) +
        ggplot2::scale_fill_gradientn(
            colors = color_scale,
            limits = c(0, max(Q, na.rm = TRUE)),
            name   = legend_title
        ) +
        ggplot2::scale_y_discrete(limits = rev(group_order)) +
        ggplot2::scale_x_discrete(limits = group_order) +
        ggplot2::labs(
            x     = "Destination group",
            y     = "Source group",
            title = title
        ) +
        ggplot2::theme_bw(base_size = text_size) +
        ggplot2::theme(
            axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid    = ggplot2::element_blank(),
            plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold")
        )

    # stability text on diagonal tiles
    if (show_stability) {
        p <- p + ggplot2::geom_text(
            data         = diag_df,
            mapping      = ggplot2::aes(x = dest, y = source, label = label),
            size         = stability_size,
            color        = "white",
            fontface     = "bold",
            inherit.aes  = FALSE
        )
    }

    # heavier border around diagonal tiles
    if (diagonal_border) {
        p <- p + ggplot2::geom_tile(
            data        = diag_df,
            mapping     = ggplot2::aes(x = dest, y = source),
            fill        = NA,
            color       = "black",
            linewidth   = 0.8,
            inherit.aes = FALSE
        )
    }

    return(p)
}
