# ============================================================
# COMPACT Alert System heatmap — annotated standalone preview
# ============================================================
#
# This version is compatible with AlertSystemScore_UPDATED.R. It preserves the
# existing selection and plotting behavior while making the response
# transformation explicit. Major additions are marked "UPDATED".
#
# Suggested package location after review: R/AlertHeatmap.R


#' Plot paired Alert System pathway responses
#'
#' Pathways are selected independently on the positive and negative sides using
#' an evidence-aware priority metric. Tile color displays either the raw paired
#' pathway effect or its null-standardized value. Text annotations display
#' response prevalence (`r`) and directional consistency (`c`).
#'
#' UPDATED: UCell difference and UCell log2-ratio results are labeled
#' explicitly, and mixed transformations are rejected rather than silently
#' plotted together.
#'
#' @param seurat_obj Seurat object processed by `AlertSystemScore()`.
#' @param misc_name Name of the Alert System state in `seurat_obj@misc`.
#' @param score_method Score channel: `"expression"` or `"ucell"`.
#' @param top_n Number of pathways selected separately from each requested
#'   effect side.
#' @param side One of `"both"`, `"positive"`, or `"negative"`.
#' @param select_by Ranking metric: `"widespread_priority"`,
#'   `"focal_priority"`, `"abs_effect"`, `"null_z"`,
#'   `"robust_fraction"`, or `"response_rate"`.
#' @param fill_by Tile fill: raw `"effect"` or calibrated `"null_z"`.
#' @param reorder_by Order selected pathways by mean effect or selection score.
#' @param scale_mode `"none"` retains the native response scale; `"unit"`
#'   divides by a robust magnitude and clips display values to `[-1, 1]`.
#' @param unit_q Quantile of absolute displayed values used by unit scaling.
#' @param show_annotations Display `r` and `c` in each tile.
#' @param show_fdr Display `*` when FDR is at or below `fdr_cutoff`.
#' @param fdr_cutoff FDR threshold for the symbol.
#' @param fill_limits Optional explicit symmetric or asymmetric fill limits.
#'
#' @return A ggplot object. Selection data and display metadata are attached as
#'   attributes: `selected_pathways`, `selection_stats`, `change_measure`,
#'   `response_threshold`, and `fill_limits`.
#' @export
AlertHeatmap <- function(
  seurat_obj,
  misc_name = "AlertSystem",
  score_method = c("expression", "ucell"),
  top_n = 5L,
  side = c("both", "positive", "negative"),
  select_by = c(
    "widespread_priority", "focal_priority", "abs_effect",
    "null_z", "robust_fraction", "response_rate"
  ),
  fill_by = c("effect", "null_z"),
  reorder_by = c("effect", "select"),
  scale_mode = c("none", "unit"),
  unit_q = 0.9,
  show_annotations = TRUE,
  show_fdr = TRUE,
  fdr_cutoff = 0.05,
  fill_limits = NULL
) {
  score_method <- match.arg(score_method)
  side <- match.arg(side)
  select_by <- match.arg(select_by)
  fill_by <- match.arg(fill_by)
  reorder_by <- match.arg(reorder_by)
  scale_mode <- match.arg(scale_mode)

  # ---------------------------------------------------------------------------
  # 1. Validate the stored state and display arguments.
  # ---------------------------------------------------------------------------
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }
  state <- seurat_obj@misc[[misc_name]]
  if (is.null(state) || is.null(state$metrics)) {
    stop(
      "No Alert System metrics found in @misc[['", misc_name,
      "']]. Run AlertSystemScore() first.",
      call. = FALSE
    )
  }
  if (!is.data.frame(state$metrics)) {
    stop("Stored Alert System metrics must be a data frame.", call. = FALSE)
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1 || top_n != as.integer(top_n)) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(unit_q) || length(unit_q) != 1L || !is.finite(unit_q) ||
      unit_q <= 0 || unit_q > 1) {
    stop("unit_q must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L ||
      !is.finite(fdr_cutoff) || fdr_cutoff < 0 || fdr_cutoff > 1) {
    stop("fdr_cutoff must be in [0, 1].", call. = FALSE)
  }
  logical_args <- list(
    show_annotations = show_annotations,
    show_fdr = show_fdr
  )
  bad_logical <- vapply(
    logical_args,
    function(x) !is.logical(x) || length(x) != 1L || is.na(x),
    logical(1)
  )
  if (any(bad_logical)) {
    stop(
      paste(names(logical_args)[bad_logical], collapse = ", "),
      " must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  if (!is.null(fill_limits)) {
    if (!is.numeric(fill_limits) || length(fill_limits) != 2L ||
        any(!is.finite(fill_limits)) || fill_limits[1] >= fill_limits[2]) {
      stop("fill_limits must be NULL or two increasing finite numbers.",
           call. = FALSE)
    }
  }

  metrics <- state$metrics
  required <- c(
    "group", "pathway", "score_method", "effect", "response_threshold",
    "response_rate", "directional_consistency", "robust_fraction",
    "widespread_priority", "focal_priority"
  )
  missing_columns <- setdiff(required, colnames(metrics))
  if (length(missing_columns) > 0L) {
    stop(
      "Alert System metrics are missing column(s): ",
      paste(missing_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }

  # Keep only the requested scoring channel before resolving its transformation.
  df <- metrics[metrics$score_method == score_method, , drop = FALSE]
  if (nrow(df) == 0L) {
    available <- unique(metrics$score_method)
    stop(
      "No metrics for score_method = '", score_method,
      "'. Available: ", paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 2. UPDATED: determine and validate the response transformation.
  # ---------------------------------------------------------------------------
  # Updated results carry change_measure in every row. For the immediately
  # preceding scorer, expression and UCell were both post-minus-pre differences,
  # so the fallback remains backward compatible with that version.
  if ("change_measure" %in% colnames(df)) {
    change_measures <- unique(as.character(df$change_measure))
    change_measures <- change_measures[
      !is.na(change_measures) & nzchar(change_measures)
    ]
  } else if (score_method == "expression") {
    change_measures <- "difference"
  } else if (!is.null(state$params$ucell_change)) {
    change_measures <- as.character(state$params$ucell_change)
  } else {
    change_measures <- "difference"
  }

  if (length(change_measures) != 1L) {
    stop(
      "The selected metrics contain multiple or undefined change measures: ",
      paste(change_measures, collapse = ", "),
      ". Plot each transformation separately.",
      call. = FALSE
    )
  }
  change_measure <- change_measures[[1L]]
  allowed_measures <- if (score_method == "expression") {
    "difference"
  } else {
    c("difference", "log2_ratio")
  }
  if (!change_measure %in% allowed_measures) {
    stop(
      "Unsupported change_measure '", change_measure,
      "' for score_method = '", score_method, "'.",
      call. = FALSE
    )
  }

  thresholds <- unique(df$response_threshold[is.finite(df$response_threshold)])
  if (length(thresholds) != 1L) {
    stop(
      "The selected metrics contain multiple or undefined response thresholds. ",
      "Plot analyses with different thresholds separately.",
      call. = FALSE
    )
  }
  response_threshold <- thresholds[[1L]]

  if (fill_by == "null_z" &&
      (!"null_z" %in% colnames(df) || all(!is.finite(df$null_z)))) {
    stop(
      "fill_by = 'null_z' requires calibrated null_z values. Supply sample.by ",
      "and enable sample-level calibration in AlertSystemScore().",
      call. = FALSE
    )
  }
  if (select_by == "null_z" &&
      (!"null_z" %in% colnames(df) || all(!is.finite(df$null_z)))) {
    stop("select_by = 'null_z' requires calibrated null_z values.",
         call. = FALSE)
  }

  # Small helpers prevent missing/non-finite values from breaking pathway-level
  # selection summaries.
  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  }
  safe_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
  }

  # ---------------------------------------------------------------------------
  # 3. Select pathways independently from positive and negative effect sides.
  # ---------------------------------------------------------------------------
  # The pathway sign uses mean effect across displayed groups. Ranking uses the
  # strongest group so a localized response is not averaged out of selection.
  split_pathways <- split(seq_len(nrow(df)), df$pathway)
  selection_stats <- do.call(rbind, lapply(names(split_pathways), function(pathway) {
    x <- df[split_pathways[[pathway]], , drop = FALSE]
    data.frame(
      pathway = pathway,
      mean_effect = safe_mean(x$effect),
      abs_effect = safe_max(abs(x$effect)),
      widespread_priority = safe_max(x$widespread_priority),
      focal_priority = safe_max(x$focal_priority),
      null_z = if ("null_z" %in% colnames(x)) {
        safe_max(abs(x$null_z))
      } else {
        NA_real_
      },
      robust_fraction = safe_max(x$robust_fraction),
      response_rate = safe_max(x$response_rate),
      stringsAsFactors = FALSE
    )
  }))
  rownames(selection_stats) <- NULL
  selection_stats$rank_value <- selection_stats[[select_by]]
  selection_stats$rank_value[!is.finite(selection_stats$rank_value)] <- -Inf

  choose_side <- function(direction) {
    keep <- if (direction == "positive") {
      selection_stats$mean_effect > 0
    } else {
      selection_stats$mean_effect < 0
    }
    candidates <- selection_stats[keep, , drop = FALSE]
    candidates <- candidates[
      order(candidates$rank_value, decreasing = TRUE), , drop = FALSE
    ]
    head(candidates$pathway, as.integer(top_n))
  }

  selected <- character(0)
  if (side %in% c("both", "positive")) {
    selected <- c(selected, choose_side("positive"))
  }
  if (side %in% c("both", "negative")) {
    selected <- c(selected, choose_side("negative"))
  }
  selected <- unique(selected)
  if (length(selected) == 0L) {
    stop("No pathways were selected for side = '", side, "'.",
         call. = FALSE)
  }

  plot_df <- df[df$pathway %in% selected, , drop = FALSE]
  selected_stats <- selection_stats[
    match(selected, selection_stats$pathway), , drop = FALSE
  ]

  ordering <- if (reorder_by == "select") {
    selected_stats$pathway[
      order(selected_stats$rank_value, decreasing = TRUE)
    ]
  } else {
    selected_stats$pathway[order(selected_stats$mean_effect)]
  }
  plot_df$pathway <- factor(plot_df$pathway, levels = ordering)
  plot_df$group <- factor(plot_df$group, levels = unique(df$group))

  # ---------------------------------------------------------------------------
  # 4. Resolve raw values, transformation-aware labels, and color scaling.
  # ---------------------------------------------------------------------------
  plot_df$plot_value <- plot_df[[fill_by]]

  # UPDATED: labels now distinguish absolute UCell differences from relative
  # log2-ratio responses. Expression remains an average log-expression delta.
  measure_label <- if (score_method == "expression") {
    "mean expression difference (post - pre)"
  } else if (change_measure == "log2_ratio") {
    pseudocount <- state$params$ucell_pseudocount
    if (!is.null(pseudocount) && length(pseudocount) == 1L &&
        is.finite(pseudocount)) {
      paste0(
        "UCell log2 ratio (post/pre; pseudocount=",
        signif(pseudocount, 3), ")"
      )
    } else {
      "UCell log2 ratio (post/pre)"
    }
  } else {
    "UCell difference (post - pre)"
  }

  legend_title <- if (fill_by == "null_z") {
    "Null-standardized effect"
  } else if (score_method == "expression") {
    "Mean expression difference\n(post - pre)"
  } else if (change_measure == "log2_ratio") {
    "UCell log2 ratio\n(post/pre)"
  } else {
    "UCell difference\n(post - pre)"
  }

  if (scale_mode == "unit") {
    scale_value <- as.numeric(stats::quantile(
      abs(plot_df$plot_value), unit_q, na.rm = TRUE, names = FALSE
    ))
    if (is.finite(scale_value) && scale_value > 0) {
      plot_df$plot_value <- pmax(
        -1, pmin(1, plot_df$plot_value / scale_value)
      )
      legend_title <- paste0(legend_title, "\n(scaled)")
    }
  }

  if (is.null(fill_limits)) {
    max_abs <- max(abs(plot_df$plot_value), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
    fill_limits <- c(-max_abs, max_abs)
  }

  # ---------------------------------------------------------------------------
  # 5. Build tile annotations.
  # ---------------------------------------------------------------------------
  plot_df$label <- NA_character_
  if (isTRUE(show_annotations)) {
    consistency <- ifelse(
      is.finite(plot_df$directional_consistency),
      sprintf("%.2f", plot_df$directional_consistency),
      "NA"
    )
    plot_df$label <- paste0(
      "r=", sprintf("%.2f", plot_df$response_rate),
      "\nc=", consistency
    )
  }
  if (isTRUE(show_fdr) && "fdr" %in% colnames(plot_df)) {
    significant <- is.finite(plot_df$fdr) & plot_df$fdr <= fdr_cutoff
    if (isTRUE(show_annotations)) {
      plot_df$label[significant] <- paste0(
        plot_df$label[significant], "\n*"
      )
    } else {
      plot_df$label[significant] <- "*"
    }
  }

  group_label <- "group"
  if (!is.null(state$params$group.by)) group_label <- state$params$group.by
  title <- switch(
    side,
    both = "Pathway Responses: Positive and Negative",
    positive = "Pathway Responses: Positive",
    negative = "Pathway Responses: Negative"
  )

  # UPDATED: the subtitle documents the transformation and the threshold used
  # to calculate r, c, robust_fraction, and the priority metrics.
  subtitle <- paste0(
    "Measure: ", measure_label,
    "; response threshold=", signif(response_threshold, 3),
    "; selected by ", select_by,
    "; r=response rate, c=directional consistency",
    if (isTRUE(show_fdr)) ", *=FDR threshold" else ""
  )

  # ---------------------------------------------------------------------------
  # 6. Construct the heatmap; selection calculations above remain unchanged.
  # ---------------------------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = group, y = pathway, fill = plot_value)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = "#1682C0",
      mid = "#FFF5CD",
      high = "#A02248",
      midpoint = 0,
      limits = fill_limits,
      oob = scales::squish,
      name = legend_title
    ) +
    ggplot2::labs(
      x = group_label,
      y = "Pathway",
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (isTRUE(show_annotations) ||
      (isTRUE(show_fdr) && any(!is.na(plot_df$label)))) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 2.5,
      lineheight = 0.85,
      na.rm = TRUE
    )
  }

  # Attributes make the automatic selection and display scale auditable.
  attr(p, "selected_pathways") <- selected
  attr(p, "selection_stats") <- selected_stats
  attr(p, "change_measure") <- change_measure
  attr(p, "response_threshold") <- response_threshold
  attr(p, "fill_limits") <- fill_limits
  p
}
