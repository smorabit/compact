# ============================================================
# COMPACT AlertSystem
# Paired, sensitivity-oriented pathway scoring
# ============================================================

#' Score subtle pathway responses after an in-silico perturbation
#'
#' \code{AlertSystemScore()} is an experimental, paired redesign of
#' \code{AlertSystemScore()}. It compares every cell with itself and supports
#' two complementary score channels:
#' \itemize{
#'   \item \code{"expression"}: the mean post-minus-pre log-expression across
#'         genes in a pathway. This is sensitive to coordinated changes that do
#'         not alter within-cell gene ranks.
#'   \item \code{"ucell"}: the post-minus-pre UCell score. This is more robust
#'         to expression scale but less sensitive to small rank-preserving
#'         changes.
#' }
#'
#' The primary \code{effect} is a paired difference, not a ratio. Response
#' prevalence and directional consistency are reported separately. When
#' \code{sample.by} is supplied, sample-level sign-flip calibration produces
#' an empirical P value, null-standardized effect, and FDR without treating
#' every cell as an independent biological replicate.
#'
#' @param seurat_obj A Seurat object containing baseline and perturbation
#'   assays.
#' @param gene_sets Named list of pathway gene vectors.
#' @param perturbation_name Perturbation assay name.
#' @param baseline_assay Baseline assay name. Default \code{"RNA"}.
#' @param group.by Optional metadata column for cell-type/state summaries.
#' @param sample.by Optional metadata column identifying biological samples.
#'   Strongly recommended when inferential calibration is desired.
#' @param layer Log-normalized layer used by the expression score.
#' @param min_genes Minimum detected genes required for a pathway.
#' @param score_methods Any of \code{c("expression", "ucell")}. If UCell is
#'   unavailable, the expression channel is retained and a warning is issued.
#' @param summary_fun Primary within-group summary: \code{"median"},
#'   \code{"trimmed_mean"}, or \code{"mean"}.
#' @param trim Fraction trimmed from each tail for \code{"trimmed_mean"}.
#' @param expression_response_threshold Minimum absolute expression delta used
#'   to define a responding cell.
#' @param ucell_response_threshold Minimum absolute UCell delta used to define
#'   a responding cell.
#' @param null_unit One of \code{"auto"}, \code{"sample"}, \code{"cell"}, or
#'   \code{"none"}. \code{"auto"} uses samples when \code{sample.by} is set
#'   and otherwise reports descriptive metrics without inferential P values.
#' @param n_permutations Number of sign-flip null draws.
#' @param seed Random seed for null calibration.
#' @param leading_n Number of genes reported as leading contributors.
#' @param gene_change_tolerance Tolerance for inferring directly changed counts
#'   and propagated data-layer genes.
#' @param misc_name Misc slot used for storing results.
#' @param return_details Return a result list instead of only the Seurat object.
#' @param verbose Emit progress messages.
#'
#' @return A Seurat object, or a list containing \code{seurat},
#'   \code{per_cell}, \code{metrics}, \code{gene_sets}, and \code{params}.
#'
#' @export
AlertSystemScore <- function(
  seurat_obj,
  gene_sets,
  perturbation_name,
  baseline_assay = "RNA",
  group.by = NULL,
  sample.by = NULL,
  layer = "data",
  min_genes = 5L,
  score_methods = c("expression", "ucell"),
  summary_fun = c("trimmed_mean", "median", "mean"),
  trim = 0.1,
  expression_response_threshold = 0.01,
  ucell_response_threshold = 0.005,
  null_unit = c("auto", "sample", "cell", "none"),
  n_permutations = 1000L,
  seed = 1L,
  leading_n = 5L,
  gene_change_tolerance = 1e-8,
  misc_name = "AlertSystem",
  return_details = FALSE,
  verbose = TRUE
) {
  summary_fun <- match.arg(summary_fun)
  null_unit <- match.arg(null_unit)

  say <- function(...) {
    if (isTRUE(verbose)) message("[COMPACT:Alert] ", ...)
  }
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets)) ||
      any(!nzchar(names(gene_sets)))) {
    stop("gene_sets must be a named list with non-empty pathway names.")
  }
  if (anyDuplicated(names(gene_sets))) {
    stop("gene_sets contains duplicated pathway names.")
  }
  if (!is.numeric(min_genes) || length(min_genes) != 1L ||
      !is.finite(min_genes) || min_genes < 1) {
    stop("min_genes must be a positive integer.")
  }
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim < 0 || trim >= 0.5) {
    stop("trim must be in [0, 0.5).")
  }
  for (item in list(
    expression_response_threshold = expression_response_threshold,
    ucell_response_threshold = ucell_response_threshold,
    gene_change_tolerance = gene_change_tolerance
  )) {
    if (!is.numeric(item) || length(item) != 1L || !is.finite(item) || item < 0) {
      stop("Response thresholds and gene_change_tolerance must be finite and non-negative.")
    }
  }

  assays <- Seurat::Assays(seurat_obj)
  missing_assays <- setdiff(c(baseline_assay, perturbation_name), assays)
  if (length(missing_assays) > 0L) {
    stop("Missing assay(s): ", paste(missing_assays, collapse = ", "), ".")
  }

  score_methods <- unique(as.character(score_methods))
  bad_methods <- setdiff(score_methods, c("expression", "ucell"))
  if (length(score_methods) == 0L || length(bad_methods) > 0L) {
    stop("score_methods must contain 'expression', 'ucell', or both.")
  }
  if ("ucell" %in% score_methods && !requireNamespace("UCell", quietly = TRUE)) {
    warning(
      "UCell is not installed; continuing with score_methods = 'expression'.",
      call. = FALSE
    )
    score_methods <- setdiff(score_methods, "ucell")
  }
  if (length(score_methods) == 0L) {
    stop("No available score method remains.")
  }

  metadata <- seurat_obj@meta.data
  if (!is.null(group.by) && !group.by %in% colnames(metadata)) {
    stop("group.by column not found: ", group.by)
  }
  if (!is.null(sample.by) && !sample.by %in% colnames(metadata)) {
    stop("sample.by column not found: ", sample.by)
  }

  is_v5 <- inherits(seurat_obj[[baseline_assay]], "Assay5")
  get_layer <- function(assay, layer_name) {
    out <- if (is_v5) {
      Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer_name)
    } else {
      Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer_name)
    }
    if (nrow(out) == 0L || ncol(out) == 0L) {
      stop("Assay '", assay, "' has an empty '", layer_name, "' layer.")
    }
    out
  }

  say("Loading paired baseline and perturbation data.")
  baseline <- get_layer(baseline_assay, layer)
  perturbed <- get_layer(perturbation_name, layer)

  features <- rownames(baseline)
  cells <- colnames(baseline)
  if (!setequal(cells, colnames(perturbed))) {
    stop("Baseline and perturbation assays contain different cells.")
  }
  missing_features <- setdiff(features, rownames(perturbed))
  if (length(missing_features) > 0L) {
    stop(
      "Perturbation assay is missing ", length(missing_features),
      " baseline feature(s)."
    )
  }
  perturbed <- perturbed[features, cells, drop = FALSE]
  metadata <- metadata[cells, , drop = FALSE]

  gene_sets_original_n <- vapply(
    gene_sets, function(x) length(unique(as.character(x))), integer(1L)
  )
  gene_sets_filtered <- lapply(
    gene_sets,
    function(x) intersect(unique(as.character(x)), features)
  )
  gene_sets_filtered <- gene_sets_filtered[
    vapply(gene_sets_filtered, length, integer(1L)) >= as.integer(min_genes)
  ]
  if (length(gene_sets_filtered) == 0L) {
    stop("No gene sets retained after feature filtering.")
  }

  group_vec <- if (is.null(group.by)) {
    rep("all_cells", length(cells))
  } else {
    as.character(metadata[[group.by]])
  }
  if (anyNA(group_vec)) {
    group_vec[is.na(group_vec)] <- "NA"
  }
  sample_vec <- if (is.null(sample.by)) {
    rep(NA_character_, length(cells))
  } else {
    as.character(metadata[[sample.by]])
  }
  if (!is.null(sample.by) && anyNA(sample_vec)) {
    stop("sample.by contains missing values.")
  }

  if (null_unit == "auto") {
    null_unit <- if (is.null(sample.by)) "none" else "sample"
  }
  if (null_unit == "sample" && is.null(sample.by)) {
    stop("null_unit = 'sample' requires sample.by.")
  }
  if (null_unit == "cell") {
    warning(
      "Cell-level null calibration treats cells as exchangeable and is descriptive, ",
      "not biological-replicate inference. Prefer sample.by with null_unit = 'sample'.",
      call. = FALSE
    )
  }
  if (null_unit != "none" &&
      (!is.numeric(n_permutations) || length(n_permutations) != 1L ||
       !is.finite(n_permutations) || n_permutations < 100)) {
    stop("n_permutations must be at least 100 when null calibration is enabled.")
  }

  summarize_values <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    switch(
      summary_fun,
      median = stats::median(x),
      trimmed_mean = mean(x, trim = trim),
      mean = mean(x)
    )
  }

  # Data-layer delta is retained for pathway scoring and leading-gene analysis.
  gene_delta <- perturbed - baseline

  pathway_scores <- list()
  if ("expression" %in% score_methods) {
    say("Computing paired expression-sensitive pathway deltas.")
    pathway_scores$expression <- lapply(gene_sets_filtered, function(genes) {
      pre <- Matrix::colMeans(baseline[genes, , drop = FALSE])
      post <- Matrix::colMeans(perturbed[genes, , drop = FALSE])
      list(pre = as.numeric(pre), post = as.numeric(post), delta = as.numeric(post - pre))
    })
  }

  if ("ucell" %in% score_methods) {
    say("Computing paired UCell pathway deltas.")
    original_meta <- seurat_obj@meta.data
    suffix_pre <- paste0("__AlertPre", abs(as.integer(seed)))
    suffix_post <- paste0("__AlertPost", abs(as.integer(seed)))

    seurat_obj <- UCell::AddModuleScore_UCell(
      seurat_obj,
      features = gene_sets_filtered,
      name = suffix_pre,
      assay = baseline_assay
    )
    seurat_obj <- UCell::AddModuleScore_UCell(
      seurat_obj,
      features = gene_sets_filtered,
      name = suffix_post,
      assay = perturbation_name
    )

    scored_meta <- seurat_obj@meta.data[cells, , drop = FALSE]
    pre_cols <- grep(paste0(suffix_pre, "$"), colnames(scored_meta), value = TRUE)
    post_cols <- grep(paste0(suffix_post, "$"), colnames(scored_meta), value = TRUE)
    pre_names <- sub(paste0(suffix_pre, "$"), "", pre_cols)
    post_names <- sub(paste0(suffix_post, "$"), "", post_cols)
    names(pre_cols) <- pre_names
    names(post_cols) <- post_names

    missing_scores <- setdiff(names(gene_sets_filtered), intersect(pre_names, post_names))
    if (length(missing_scores) > 0L) {
      seurat_obj@meta.data <- original_meta
      stop(
        "UCell did not return matched scores for: ",
        paste(missing_scores, collapse = ", "), "."
      )
    }

    pathway_scores$ucell <- lapply(names(gene_sets_filtered), function(pathway) {
      pre <- as.numeric(scored_meta[[pre_cols[[pathway]]]])
      post <- as.numeric(scored_meta[[post_cols[[pathway]]]])
      list(pre = pre, post = post, delta = post - pre)
    })
    names(pathway_scores$ucell) <- names(gene_sets_filtered)
    seurat_obj@meta.data <- original_meta
  }

  # Infer which pathway genes were directly changed in counts and which changed
  # only in data. These are descriptive provenance fields, not score filters.
  direct_genes <- character(0)
  propagated_genes <- character(0)
  counts_ok <- TRUE
  base_counts <- pert_counts <- NULL
  tryCatch({
    base_counts <- get_layer(baseline_assay, "counts")
    pert_counts <- get_layer(perturbation_name, "counts")
    pert_counts <- pert_counts[rownames(base_counts), colnames(base_counts), drop = FALSE]
  }, error = function(e) {
    counts_ok <<- FALSE
  })
  if (counts_ok) {
    count_delta <- pert_counts - base_counts
    direct_genes <- rownames(base_counts)[
      Matrix::rowSums(abs(count_delta) > gene_change_tolerance) > 0
    ]
    data_changed <- features[
      Matrix::rowSums(abs(gene_delta) > gene_change_tolerance) > 0
    ]
    propagated_genes <- setdiff(data_changed, direct_genes)
  }

  say("Building per-cell paired score table.")
  per_cell_parts <- list()
  cursor <- 0L
  for (method in names(pathway_scores)) {
    for (pathway in names(gene_sets_filtered)) {
      cursor <- cursor + 1L
      score <- pathway_scores[[method]][[pathway]]
      per_cell_parts[[cursor]] <- data.frame(
        cell = cells,
        group = group_vec,
        sample = sample_vec,
        pathway = pathway,
        score_method = method,
        pre_score = score$pre,
        post_score = score$post,
        delta = score$delta,
        stringsAsFactors = FALSE
      )
    }
  }
  per_cell <- do.call(rbind, per_cell_parts)
  rownames(per_cell) <- NULL

  say("Summarizing effects, prevalence, consistency, and null calibration.")
  set.seed(seed)
  metric_parts <- list()
  cursor <- 0L

  for (method in names(pathway_scores)) {
    threshold <- if (method == "expression") {
      expression_response_threshold
    } else {
      ucell_response_threshold
    }

    for (pathway in names(gene_sets_filtered)) {
      values_all <- pathway_scores[[method]][[pathway]]$delta
      for (group_name in unique(group_vec)) {
        idx <- which(group_vec == group_name & is.finite(values_all))
        vals <- values_all[idx]
        if (length(vals) == 0L) next

        effect <- summarize_values(vals)
        effect_direction <- sign(effect)
        responder <- abs(vals) >= threshold
        tie_fraction <- mean(abs(vals) < threshold)
        response_rate <- mean(responder)

        if (effect_direction == 0 || !any(responder)) {
          directional_consistency <- NA_real_
          robust_fraction <- 0
          opposite_fraction <- 0
        } else {
          directional_consistency <- mean(
            sign(vals[responder]) == effect_direction
          )
          robust_fraction <- mean(
            responder & sign(vals) == effect_direction
          )
          opposite_fraction <- mean(
            responder & sign(vals) == -effect_direction
          )
        }
        concordance_all <- if (effect_direction == 0) {
          0
        } else {
          mean(sign(vals) == effect_direction)
        }

        sample_effects <- numeric(0)
        if (!is.null(sample.by)) {
          local_samples <- sample_vec[idx]
          sample_effects <- vapply(
            split(vals, local_samples),
            summarize_values,
            numeric(1L)
          )
          sample_effects <- sample_effects[is.finite(sample_effects)]
        }
        n_samples <- if (is.null(sample.by)) NA_integer_ else length(sample_effects)
        sample_consistency <- if (length(sample_effects) > 0L &&
                                  effect_direction != 0) {
          mean(sign(sample_effects) == effect_direction)
        } else {
          NA_real_
        }

        empirical_p <- null_z <- null_center <- null_scale <- NA_real_
        calibration_effect <- NA_real_
        if (null_unit != "none") {
          units <- if (null_unit == "sample") sample_effects else vals
          units <- units[is.finite(units)]
          if (length(units) >= 2L) {
            calibration_effect <- summarize_values(units)
            null_effects <- replicate(as.integer(n_permutations), {
              summarize_values(units * sample(c(-1, 1), length(units), replace = TRUE))
            })
            null_center <- stats::median(null_effects, na.rm = TRUE)
            null_scale <- stats::mad(null_effects, center = null_center, na.rm = TRUE)
            if (is.finite(null_scale) && null_scale > 0) {
              null_z <- (calibration_effect - null_center) / null_scale
            }
            empirical_p <- (
              1 + sum(abs(null_effects) >= abs(calibration_effect), na.rm = TRUE)
            ) / (1 + sum(is.finite(null_effects)))
          }
        }

        genes <- gene_sets_filtered[[pathway]]
        gene_effect <- Matrix::rowMeans(gene_delta[genes, idx, drop = FALSE])
        gene_effect <- as.numeric(gene_effect)
        names(gene_effect) <- genes
        lead_order <- order(abs(gene_effect), decreasing = TRUE)
        lead_order <- head(lead_order, min(as.integer(leading_n), length(lead_order)))
        leading_genes <- paste(
          paste0(names(gene_effect)[lead_order], "=", signif(gene_effect[lead_order], 3)),
          collapse = ";"
        )
        gene_direction_fraction <- if (effect_direction == 0) {
          NA_real_
        } else {
          mean(sign(gene_effect) == effect_direction)
        }

        standardized_strength <- if (is.finite(null_z)) {
          abs(null_z)
        } else if (threshold > 0) {
          abs(effect) / threshold
        } else {
          abs(effect)
        }
        consistency_signal <- if (is.finite(directional_consistency)) {
          max(0, 2 * directional_consistency - 1)
        } else {
          0
        }
        widespread_priority <- standardized_strength *
          sqrt(response_rate) * consistency_signal
        directional_values <- effect_direction * vals
        focal_tail <- if (effect_direction == 0) {
          0
        } else {
          max(0, as.numeric(stats::quantile(
            directional_values, 0.9, na.rm = TRUE, names = FALSE
          )))
        }
        focal_priority <- if (threshold > 0) {
          (focal_tail / threshold) * sqrt(robust_fraction)
        } else {
          focal_tail * sqrt(robust_fraction)
        }

        cursor <- cursor + 1L
        metric_parts[[cursor]] <- data.frame(
          group = group_name,
          pathway = pathway,
          score_method = method,
          n_cells = length(vals),
          n_samples = n_samples,
          effect = effect,
          mean_delta = mean(vals),
          median_delta = stats::median(vals),
          trimmed_mean_delta = mean(vals, trim = trim),
          q10_delta = as.numeric(stats::quantile(vals, 0.1, names = FALSE)),
          q90_delta = as.numeric(stats::quantile(vals, 0.9, names = FALSE)),
          response_threshold = threshold,
          response_rate = response_rate,
          directional_consistency = directional_consistency,
          robust_fraction = robust_fraction,
          opposite_fraction = opposite_fraction,
          tie_fraction = tie_fraction,
          concordance_all = concordance_all,
          sample_consistency = sample_consistency,
          calibration_effect = calibration_effect,
          null_center = null_center,
          null_scale = null_scale,
          null_z = null_z,
          empirical_p = empirical_p,
          genes_in_set = gene_sets_original_n[[pathway]],
          genes_used = length(genes),
          gene_coverage = length(genes) / gene_sets_original_n[[pathway]],
          n_direct_overlap = sum(genes %in% direct_genes),
          n_propagated_overlap = sum(genes %in% propagated_genes),
          gene_direction_fraction = gene_direction_fraction,
          leading_genes = leading_genes,
          widespread_priority = widespread_priority,
          focal_priority = focal_priority,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  metrics <- do.call(rbind, metric_parts)
  rownames(metrics) <- NULL
  metrics$fdr <- NA_real_
  calibration_groups <- interaction(
    metrics$group, metrics$score_method, drop = TRUE
  )
  for (cal_group in levels(calibration_groups)) {
    idx <- which(calibration_groups == cal_group & is.finite(metrics$empirical_p))
    if (length(idx) > 0L) {
      metrics$fdr[idx] <- stats::p.adjust(metrics$empirical_p[idx], method = "BH")
    }
  }

  params <- list(
    perturbation_name = perturbation_name,
    baseline_assay = baseline_assay,
    group.by = group.by,
    sample.by = sample.by,
    layer = layer,
    min_genes = min_genes,
    score_methods = names(pathway_scores),
    summary_fun = summary_fun,
    trim = trim,
    expression_response_threshold = expression_response_threshold,
    ucell_response_threshold = ucell_response_threshold,
    null_unit = null_unit,
    n_permutations = if (null_unit == "none") 0L else as.integer(n_permutations),
    seed = seed,
    leading_n = leading_n,
    timestamp = Sys.time()
  )

  seurat_obj@misc[[misc_name]] <- list(
    per_cell = per_cell,
    metrics = metrics,
    gene_sets = gene_sets_filtered,
    params = params
  )

  say(
    "Completed: ", length(gene_sets_filtered), " pathways, ",
    length(unique(group_vec)), " group(s), methods = ",
    paste(names(pathway_scores), collapse = ", "), "."
  )

  if (!isTRUE(return_details)) {
    return(seurat_obj)
  }
  list(
    seurat = seurat_obj,
    per_cell = per_cell,
    metrics = metrics,
    gene_sets = gene_sets_filtered,
    params = params
  )
}
