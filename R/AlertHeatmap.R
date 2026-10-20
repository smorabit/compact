# ============================================================
# COMPACT AlertSystem heatmap
# ============================================================

#' Plot sensitivity-oriented AlertSystem pathway results
#'
#' Selects pathways independently on the positive and negative sides using an
#' evidence-aware priority metric. Tile color shows either the paired pathway
#' effect or its null-standardized value. Text annotations keep response
#' prevalence and consistency visible so small but reproducible changes are not
#' confused with numerical noise.
#'
#' @param seurat_obj Seurat object processed by \code{AlertSystemScore()}.
#' @param misc_name Misc slot containing  results.
#' @param score_method Score channel: \code{"expression"} or \code{"ucell"}.
#' @param top_n Number of pathways selected per effect side.
#' @param side One of \code{"both"}, \code{"positive"}, or \code{"negative"}.
#' @param select_by Selection metric. Defaults to
#'   \code{"widespread_priority"}; alternatives are \code{"focal_priority"},
#'   \code{"abs_effect"}, \code{"null_z"}, \code{"robust_fraction"}, and
#'   \code{"response_rate"}.
#' @param fill_by Tile fill metric: raw paired \code{"effect"} or
#'   \code{"null_z"}.
#' @param reorder_by Order selected pathways by \code{"effect"} or selection
#'   priority.
#' @param scale_mode Use raw values or divide by a robust global magnitude.
#' @param unit_q Quantile used for unit scaling.
#' @param show_annotations Show response and consistency text.
#' @param show_fdr Add an FDR symbol when calibrated FDR is available.
#' @param fdr_cutoff FDR threshold for the symbol.
#' @param fill_limits Optional explicit fill limits.
#'
#' @return A ggplot object. Selected pathway names and selection statistics are
#'   attached as \code{"selected_pathways"} and \code{"selection_stats"}.
#'
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

  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  state <- seurat_obj@misc[[misc_name]]
  if (is.null(state) || is.null(state$metrics)) {
    stop(
      "No AlertSystem  metrics found in @misc[['", misc_name,
      "']]. Run AlertSystemScore() first."
    )
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1) {
    stop("top_n must be a positive integer.")
  }
  if (!is.numeric(unit_q) || length(unit_q) != 1L || !is.finite(unit_q) ||
      unit_q <= 0 || unit_q > 1) {
    stop("unit_q must be in (0, 1].")
  }

  metrics <- state$metrics
  required <- c(
    "group", "pathway", "score_method", "effect", "response_rate",
    "directional_consistency", "robust_fraction", "widespread_priority",
    "focal_priority"
  )
  missing <- setdiff(required, colnames(metrics))
  if (length(missing) > 0L) {
    stop(" metrics missing column(s): ", paste(missing, collapse = ", "), ".")
  }

  df <- metrics[metrics$score_method == score_method, , drop = FALSE]
  if (nrow(df) == 0L) {
    available <- unique(metrics$score_method)
    stop(
      "No metrics for score_method = '", score_method,
      "'. Available: ", paste(available, collapse = ", "), "."
    )
  }
  if (fill_by == "null_z" && (!"null_z" %in% colnames(df) ||
                              all(!is.finite(df$null_z)))) {
    stop(
      "fill_by = 'null_z' requires calibrated null_z values. Supply sample.by ",
      "and enable sample-level calibration in AlertSystemScore()."
    )
  }
  if (select_by == "null_z" && (!"null_z" %in% colnames(df) ||
                                all(!is.finite(df$null_z)))) {
    stop(
      "select_by = 'null_z' requires calibrated null_z values."
    )
  }

  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
  }
  safe_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
  }

  # Pathway sign is based on its mean effect across displayed groups. Selection
  # priorities use the strongest group so a localized cell-type response is not
  # diluted out of the plot.
  split_path <- split(seq_len(nrow(df)), df$pathway)
  selection_stats <- do.call(rbind, lapply(names(split_path), function(pathway) {
    x <- df[split_path[[pathway]], , drop = FALSE]
    data.frame(
      pathway = pathway,
      mean_effect = safe_mean(x$effect),
      abs_effect = safe_max(abs(x$effect)),
      widespread_priority = safe_max(x$widespread_priority),
      focal_priority = safe_max(x$focal_priority),
      null_z = if ("null_z" %in% colnames(x)) safe_max(abs(x$null_z)) else NA_real_,
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
    stop("No pathways were selected for side = '", side, "'.")
  }

  plot_df <- df[df$pathway %in% selected, , drop = FALSE]
  selected_stats <- selection_stats[
    match(selected, selection_stats$pathway), , drop = FALSE
  ]

  if (reorder_by == "select") {
    ordering <- selected_stats$pathway[
      order(selected_stats$rank_value, decreasing = TRUE)
    ]
  } else {
    ordering <- selected_stats$pathway[
      order(selected_stats$mean_effect)
    ]
  }
  plot_df$pathway <- factor(plot_df$pathway, levels = ordering)
  plot_df$group <- factor(plot_df$group, levels = unique(df$group))

  plot_df$plot_value <- plot_df[[fill_by]]
  legend_title <- if (fill_by == "effect") {
    paste0(score_method, " paired effect")
  } else {
    "null-standardized effect"
  }

  if (scale_mode == "unit") {
    scale_value <- as.numeric(stats::quantile(
      abs(plot_df$plot_value), unit_q, na.rm = TRUE, names = FALSE
    ))
    if (is.finite(scale_value) && scale_value > 0) {
      plot_df$plot_value <- pmax(-1, pmin(1, plot_df$plot_value / scale_value))
      legend_title <- paste0(legend_title, " (scaled)")
    }
  }

  if (is.null(fill_limits)) {
    max_abs <- max(abs(plot_df$plot_value), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
    fill_limits <- c(-max_abs, max_abs)
  }

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
      plot_df$label[significant] <- paste0(plot_df$label[significant], "\n*")
    } else {
      plot_df$label[significant] <- "*"
    }
  }

  group_label <- "group"
  if (!is.null(state$params$group.by)) group_label <- state$params$group.by
  title <- switch(
    side,
    both = "Subtle Pathway Responses: Positive and Negative",
    positive = "Subtle Pathway Responses: Positive",
    negative = "Subtle Pathway Responses: Negative"
  )

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
      subtitle = paste0(
        "Selected by ", select_by,
        "; r=response rate, c=directional consistency",
        if (isTRUE(show_fdr)) ", *=FDR threshold" else ""
      )
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

  attr(p, "selected_pathways") <- selected
  attr(p, "selection_stats") <- selected_stats
  p
}
