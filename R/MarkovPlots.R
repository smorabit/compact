
#' Plot a Markov Score on a Dimensionality Reduction Embedding
#'
#' Visualizes a Markov analysis score (from \code{\link{PredictAttractors}},
#' \code{\link{PredictPerturbationTime}}, \code{\link{PredictFates}}, or
#' \code{\link{PredictCommitment}}) as a colored scatter plot on any 2D
#' reduction stored in the Seurat object. Optionally marks source and/or sink
#' cell populations with a centroid diamond and greys out unconverged cells so
#' they do not contaminate the color scale.
#'
#' @param seurat_obj A Seurat object containing the dimensionality reduction
#'   and the score in \code{@meta.data}.
#' @param feature Character. Metadata column name to color cells by (e.g.,
#'   \code{"attractor_score"}, \code{"perturbation_pseudotime"},
#'   \code{"commitment_score"}).
#' @param reduction Character. Name of the reduction stored in
#'   \code{seurat_obj@reductions}. Default: \code{"umap"}.
#' @param source_cells Character vector or NULL. Cell barcodes belonging to the
#'   source population. Their mean embedding coordinates are computed and drawn
#'   as a single centroid diamond (shape 18). NULL = no source centroid.
#' @param sink_cells Character vector or NULL. Cell barcodes belonging to the
#'   sink population. Their mean embedding coordinates are computed and drawn as
#'   a single centroid diamond (shape 18). NULL = no sink centroid.
#' @param source_label Character. Text label shown next to the source centroid
#'   via \code{ggrepel::geom_label_repel}. Only used when
#'   \code{label_centroids = TRUE} and \code{source_cells} is not NULL.
#'   Default: \code{"Source"}.
#' @param sink_label Character. Text label shown next to the sink centroid via
#'   \code{ggrepel::geom_label_repel}. Only used when
#'   \code{label_centroids = TRUE} and \code{sink_cells} is not NULL.
#'   Default: \code{"Sink"}.
#' @param centroid_size Numeric. Size of the centroid diamond marker.
#'   Default: \code{5}.
#' @param centroid_color Character. Color of the centroid diamond marker and its
#'   label. Default: \code{"black"}.
#' @param label_centroids Logical. If \code{TRUE}, add a text label to each
#'   centroid using \code{ggrepel::geom_label_repel}. Default: \code{TRUE}.
#' @param label_size Numeric. Font size for centroid labels. Default: \code{3}.
#' @param unconverged_cells Character vector or NULL. Cell barcodes that hit
#'   \code{max_iter} in \code{\link{PredictPerturbationTime}}. These cells are
#'   drawn in \code{unconverged_color} below all other layers and excluded from
#'   the color scale limits. NULL = no special handling.
#' @param unconverged_color Character. Color for unconverged cells.
#'   Default: \code{"grey70"}.
#' @param color_scale Character vector. Two or more colors for the sequential
#'   color gradient passed to \code{scale_color_gradientn}.
#'   Default: \code{c("lightgrey", "#2166AC")}.
#' @param pt_size Numeric. Point size for the main cell layer. Default: \code{0.8}.
#' @param alpha Numeric. Transparency for the main cell layer. Default: \code{0.8}.
#' @param legend_title Character or NULL. Legend title. If NULL, defaults to
#'   \code{feature}.
#' @param title Character or NULL. Plot title. If NULL, defaults to
#'   \code{feature}.
#' @param ... Additional arguments (currently unused; reserved for future use).
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' Layers are drawn bottom to top: unconverged cells (grey) → converged cells
#' (colored gradient) → source centroid diamond → sink centroid diamond →
#' centroid labels (ggrepel). A single centroid per group is computed from the
#' mean of the source/sink cells' embedding coordinates, so the marker position
#' represents the group center rather than cluttering the plot with one marker
#' per cell.
#'
#' The color scale limits are computed from converged cells only, so
#' unconverged cells (which often carry a meaningless \code{max_iter} value)
#' do not distort the range.
#'
#' @seealso \code{\link{PredictAttractors}}, \code{\link{PredictPerturbationTime}},
#'   \code{\link{PredictFates}}, \code{\link{PredictCommitment}},
#'   \code{\link{PlotMarkovScatter}}
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn labs
#' @importFrom ggrepel geom_label_repel
#'
#' @export
PlotMarkovEmbedding <- function(
    seurat_obj,
    feature,
    reduction         = "umap",
    source_cells      = NULL,
    sink_cells        = NULL,
    source_label      = "Source",
    sink_label        = "Sink",
    centroid_size     = 5,
    centroid_color    = "black",
    label_centroids   = TRUE,
    label_size        = 3,
    unconverged_cells = NULL,
    unconverged_color = "grey70",
    color_scale       = c("lightgrey", "#2166AC"),
    pt_size           = 0.8,
    alpha             = 0.8,
    legend_title      = NULL,
    title             = NULL,
    ...
) {
    # validate feature
    if (!feature %in% colnames(seurat_obj@meta.data)) {
        stop(sprintf(
            "feature '%s' not found in seurat_obj@meta.data.",
            feature
        ))
    }

    # validate reduction
    if (!reduction %in% names(seurat_obj@reductions)) {
        stop(sprintf(
            "reduction '%s' not found in seurat_obj@reductions. Available: %s",
            reduction, paste(names(seurat_obj@reductions), collapse = ", ")
        ))
    }

    legend_title <- if (is.null(legend_title)) feature else legend_title
    title        <- if (is.null(title)) feature else title

    # extract embedding coordinates
    emb <- as.data.frame(Seurat::Embeddings(seurat_obj, reduction = reduction))
    colnames(emb)[1:2] <- c("dim1", "dim2")
    emb$cell  <- rownames(emb)
    emb$score <- seurat_obj@meta.data[rownames(emb), feature]

    # split unconverged from converged cells
    if (!is.null(unconverged_cells)) {
        emb_unconverged <- emb[emb$cell %in% unconverged_cells, ]
        emb_converged   <- emb[!emb$cell %in% unconverged_cells, ]
    } else {
        emb_unconverged <- emb[0, ]
        emb_converged   <- emb
    }

    # color scale limits exclude unconverged cells
    score_limits <- range(emb_converged$score, na.rm = TRUE)

    # base plot — no data yet so we can add layers in the correct order
    p <- ggplot2::ggplot()

    # layer 1: unconverged cells (grey, bottom-most)
    if (nrow(emb_unconverged) > 0) {
        p <- p + ggplot2::geom_point(
            data  = emb_unconverged,
            ggplot2::aes(x = dim1, y = dim2),
            color = unconverged_color,
            size  = pt_size,
            alpha = alpha
        )
    }

    # layer 2: converged cells colored by score — sorted so high-value cells render on top
    emb_converged <- emb_converged[order(emb_converged$score, na.last = FALSE), ]

    p <- p +
        ggplot2::geom_point(
            data  = emb_converged,
            ggplot2::aes(x = dim1, y = dim2, color = score),
            size  = pt_size,
            alpha = alpha
        ) +
        ggplot2::scale_color_gradientn(
            colors = color_scale,
            limits = score_limits,
            name   = legend_title
        )

    # layers 3–4: source and sink centroid diamonds
    centroid_rows <- list()
    if (!is.null(source_cells)) {
        src_in_obj <- source_cells[source_cells %in% rownames(emb)]
        centroid_rows[["source"]] <- data.frame(
            dim1  = mean(emb[src_in_obj, "dim1"]),
            dim2  = mean(emb[src_in_obj, "dim2"]),
            label = source_label,
            stringsAsFactors = FALSE
        )
    }
    if (!is.null(sink_cells)) {
        snk_in_obj <- sink_cells[sink_cells %in% rownames(emb)]
        centroid_rows[["sink"]] <- data.frame(
            dim1  = mean(emb[snk_in_obj, "dim1"]),
            dim2  = mean(emb[snk_in_obj, "dim2"]),
            label = sink_label,
            stringsAsFactors = FALSE
        )
    }

    if (length(centroid_rows) > 0) {
        centroid_df <- do.call(rbind, centroid_rows)

        p <- p + ggplot2::geom_point(
            data        = centroid_df,
            ggplot2::aes(x = dim1, y = dim2),
            shape       = 18,
            size        = centroid_size,
            color       = centroid_color,
            inherit.aes = FALSE
        )

        if (label_centroids) {
            p <- p + ggrepel::geom_label_repel(
                data        = centroid_df,
                ggplot2::aes(x = dim1, y = dim2, label = label),
                size        = label_size,
                color       = centroid_color,
                fill        = "white",
                box.padding = 0.5,
                inherit.aes = FALSE
            )
        }
    }

    p <- p +
        ggplot2::labs(title = title) +
        hdWGCNA::umap_theme()

    p
}


#' Scatter Plot of Two Markov or Biological Features
#'
#' Plots any two metadata columns as a scatter plot, colored by an optional
#' grouping variable. Designed for visualizing relationships between Markov
#' analysis scores (e.g., attractor score vs. hitting time, commitment
#' probability vs. module eigengene) but works with any numeric metadata pair.
#'
#' @param seurat_obj A Seurat object with all requested features in
#'   \code{@meta.data}.
#' @param x_feature Character. Metadata column for the x-axis.
#' @param y_feature Character. Metadata column for the y-axis.
#' @param color.by Character or NULL. Metadata column for point color. If
#'   NULL, all points are drawn in a single color (\code{"#2166AC"}).
#' @param color_palette Named character vector or NULL. Named color vector for
#'   categorical \code{color.by} groups passed to \code{scale_color_manual}.
#'   If NULL, \code{ggplot2} default colors are used. Ignored when
#'   \code{color.by} is numeric.
#' @param pt_size Numeric. Point size. Default: \code{0.8}.
#' @param alpha Numeric. Point transparency. Default: \code{0.6}.
#' @param add_smooth Logical. Add a smoothing line via \code{geom_smooth}.
#'   Default: \code{FALSE}.
#' @param smooth_method Character. Smoothing method passed to
#'   \code{geom_smooth} (e.g., \code{"loess"}, \code{"lm"}).
#'   Default: \code{"loess"}.
#' @param smooth_se Logical. Show confidence interval band around the smooth
#'   line. Default: \code{TRUE}.
#' @param hline Numeric or NULL. Y-value for a horizontal reference line.
#'   NULL = no line. Drawn under the point layer. Default: \code{NULL}.
#' @param hline_color Character. Color for the horizontal reference line.
#'   Default: \code{"grey40"}.
#' @param hline_linetype Character. Linetype for the horizontal reference line.
#'   Default: \code{"dashed"}.
#' @param vline Numeric or NULL. X-value for a vertical reference line.
#'   NULL = no line. Drawn under the point layer. Default: \code{NULL}.
#' @param vline_color Character. Color for the vertical reference line.
#'   Default: \code{"grey40"}.
#' @param vline_linetype Character. Linetype for the vertical reference line.
#'   Default: \code{"dashed"}.
#' @param x_label Character or NULL. X-axis label. If NULL, defaults to
#'   \code{x_feature}.
#' @param y_label Character or NULL. Y-axis label. If NULL, defaults to
#'   \code{y_feature}.
#' @param title Character or NULL. Plot title.
#' @param legend_title Character or NULL. Legend title. If NULL, defaults to
#'   \code{color.by}.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' Color handling adapts to the type of \code{color.by}:
#' \itemize{
#'   \item NULL: fixed single color.
#'   \item Numeric column: sequential gradient via \code{scale_color_gradientn}.
#'   \item Factor/character column: discrete colors via \code{scale_color_manual}
#'     (if \code{color_palette} provided) or \code{scale_color_discrete}.
#' }
#'
#' Reference lines (\code{hline}, \code{vline}) are drawn before the point
#' layer so they sit under the points.
#'
#' @seealso \code{\link{PredictAttractors}}, \code{\link{PredictPerturbationTime}},
#'   \code{\link{PredictCommitment}}, \code{\link{PlotMarkovEmbedding}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline geom_smooth
#'   scale_color_gradientn scale_color_manual scale_color_discrete labs theme_bw
#'
#' @export
PlotMarkovScatter <- function(
    seurat_obj,
    x_feature,
    y_feature,
    color.by          = NULL,
    color_palette     = NULL,
    pt_size           = 0.8,
    alpha             = 0.6,
    add_smooth        = FALSE,
    smooth_method     = "loess",
    smooth_se         = TRUE,
    hline             = NULL,
    hline_color       = "grey40",
    hline_linetype    = "dashed",
    vline             = NULL,
    vline_color       = "grey40",
    vline_linetype    = "dashed",
    x_label           = NULL,
    y_label           = NULL,
    title             = NULL,
    legend_title      = NULL
) {
    # validate x_feature
    if (!x_feature %in% colnames(seurat_obj@meta.data)) {
        stop(sprintf(
            "x_feature '%s' not found in seurat_obj@meta.data.",
            x_feature
        ))
    }

    # validate y_feature
    if (!y_feature %in% colnames(seurat_obj@meta.data)) {
        stop(sprintf(
            "y_feature '%s' not found in seurat_obj@meta.data.",
            y_feature
        ))
    }

    x_label      <- if (is.null(x_label)) x_feature else x_label
    y_label      <- if (is.null(y_label)) y_feature else y_label
    legend_title <- if (is.null(legend_title) && !is.null(color.by)) color.by else legend_title

    # build plotting data frame
    df <- data.frame(
        x    = seurat_obj@meta.data[[x_feature]],
        y    = seurat_obj@meta.data[[y_feature]],
        cell = colnames(seurat_obj),
        stringsAsFactors = FALSE
    )
    if (!is.null(color.by)) {
        df$color_var <- seurat_obj@meta.data[[color.by]]
    }
    df <- df[!is.na(df$x) & !is.na(df$y), ]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y))

    # reference lines drawn under points
    if (!is.null(hline)) {
        p <- p + ggplot2::geom_hline(
            yintercept = hline,
            color      = hline_color,
            linetype   = hline_linetype
        )
    }
    if (!is.null(vline)) {
        p <- p + ggplot2::geom_vline(
            xintercept = vline,
            color      = vline_color,
            linetype   = vline_linetype
        )
    }

    # point layer and color scale
    if (is.null(color.by)) {
        p <- p + ggplot2::geom_point(
            color = "#2166AC",
            size  = pt_size,
            alpha = alpha
        )
    } else if (is.numeric(df$color_var)) {
        p <- p +
            ggplot2::geom_point(
                ggplot2::aes(color = color_var),
                size  = pt_size,
                alpha = alpha
            ) +
            ggplot2::scale_color_gradientn(
                colors = c("lightgrey", "#2166AC"),
                name   = legend_title
            )
    } else {
        p <- p + ggplot2::geom_point(
            ggplot2::aes(color = color_var),
            size  = pt_size,
            alpha = alpha
        )
        if (!is.null(color_palette)) {
            p <- p + ggplot2::scale_color_manual(
                values = color_palette,
                name   = legend_title
            )
        } else {
            p <- p + ggplot2::scale_color_discrete(name = legend_title)
        }
    }

    if (add_smooth) {
        p <- p + ggplot2::geom_smooth(method = smooth_method, se = smooth_se)
    }

    p <- p +
        ggplot2::labs(x = x_label, y = y_label, title = title) +
        ggplot2::theme_bw()

    p
}
