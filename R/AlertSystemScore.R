# ============================================================
# COMPACT Alert System scorer — annotated standalone preview
# ============================================================
#
# Source AlertSystemHelper.R before this file when using out_dir, or when using
# SaveAlertSystem(), LoadAlertSystem(), or LoadAlertSystemPathways(). This scorer
# can run by itself when out_dir = NULL.
#
# Major changes relative to the current repository copy are marked "UPDATED".
# Existing V2 calculations are also annotated by section so the data flow is
# explicit. Suggested package location after review: R/AlertSystemScore.R


#' Score paired pathway responses after an in-silico perturbation
#'
#' `AlertSystemScore()` compares each cell with itself between a baseline assay
#' and a virtually perturbed assay. It supports two complementary channels:
#'
#' * `expression`: mean post-minus-pre log-normalized expression over retained
#'   pathway genes.
#' * `ucell`: either post-minus-pre UCell score (`"difference"`) or the legacy
#'   log2 ratio of post and pre UCell scores (`"log2_ratio"`).
#'
#' The cell-level pathway responses are summarized within cell groups. Response
#' prevalence, directionality, leading genes, and optional sample-level
#' sign-flip calibration are reported separately.
#'
#' @param seurat_obj Seurat object containing baseline and perturbation assays.
#' @param gene_sets Named list of pathway gene vectors.
#' @param perturbation_name Name of the virtually perturbed assay.
#' @param baseline_assay Name of the baseline assay; default `"RNA"`.
#' @param group.by Optional metadata column defining cell groups. When `NULL`,
#'   all cells form the `"all_cells"` group.
#' @param sample.by Optional metadata column identifying biological samples.
#'   Strongly recommended for inferential calibration.
#' @param layer Log-normalized layer used by expression scoring.
#' @param min_genes Minimum retained genes required for a pathway.
#' @param score_methods Any of `c("expression", "ucell")`.
#' @param ucell_change UCell response transformation. `"difference"` computes
#'   `post - pre`; `"log2_ratio"` computes
#'   `log2((post + ucell_pseudocount)/(pre + ucell_pseudocount))`.
#' @param ucell_pseudocount Positive pseudocount used only for `"log2_ratio"`.
#' @param summary_fun Group-level effect summary: `"trimmed_mean"`, `"median"`,
#'   or `"mean"`.
#' @param trim Fraction removed from each tail by the trimmed mean.
#' @param expression_response_threshold Absolute expression-delta threshold.
#' @param ucell_response_threshold Absolute UCell-response threshold. When
#'   `NULL`, defaults to `0.005` for differences and `0.1` for log2 ratios.
#' @param null_unit `"auto"`, `"sample"`, `"cell"`, or `"none"`. `"auto"`
#'   uses samples when `sample.by` is supplied and otherwise disables null
#'   calibration.
#' @param n_permutations Number of sign-flip permutations.
#' @param seed Random seed for null calibration and temporary UCell column names.
#' @param leading_n Number of leading genes recorded per pathway/group.
#' @param gene_change_tolerance Numeric tolerance used to classify directly and
#'   propagation-only changed genes.
#' @param misc_name Name used under `seurat_obj@misc`.
#' @param out_dir Optional output directory. UPDATED: when supplied, the complete
#'   state is saved as RDS and metrics are exported as CSV.
#' @param output_prefix Optional filename prefix. Defaults to a sanitized
#'   perturbation-based prefix.
#' @param overwrite_outputs Whether existing RDS/CSV outputs may be replaced.
#' @param return_details Return a result list instead of only the Seurat object.
#' @param verbose Emit progress messages.
#'
#' @return A Seurat object, or a list containing `seurat`, `per_cell`, `metrics`,
#'   `gene_sets`, `params`, `info_file`, and `metrics_file`.
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
  ucell_change = c("difference", "log2_ratio"),
  ucell_pseudocount = 1e-4,
  summary_fun = c("trimmed_mean", "median", "mean"),
  trim = 0.1,
  expression_response_threshold = 0.01,
  ucell_response_threshold = NULL,
  null_unit = c("auto", "sample", "cell", "none"),
  n_permutations = 1000L,
  seed = 1L,
  leading_n = 5L,
  gene_change_tolerance = 1e-8,
  misc_name = "AlertSystem",
  out_dir = NULL,
  output_prefix = NULL,
  overwrite_outputs = FALSE,
  return_details = FALSE,
  verbose = TRUE
) {
  summary_fun <- match.arg(summary_fun)
  null_unit <- match.arg(null_unit)
  ucell_change <- match.arg(ucell_change)

  # A local logger keeps the scoring function self-contained.
  say <- function(...) {
    if (isTRUE(verbose)) message("[COMPACT:Alert] ", ...)
  }

  # ---------------------------------------------------------------------------
  # 1. Validate inputs before loading large matrices.
  # ---------------------------------------------------------------------------
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets)) ||
      length(gene_sets) == 0L || anyNA(names(gene_sets)) ||
      any(!nzchar(names(gene_sets)))) {
    stop("gene_sets must be a non-empty named list.", call. = FALSE)
  }
  if (anyDuplicated(names(gene_sets))) {
    stop("gene_sets contains duplicated pathway names.", call. = FALSE)
  }
  if (!is.numeric(min_genes) || length(min_genes) != 1L ||
      !is.finite(min_genes) || min_genes < 1) {
    stop("min_genes must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim < 0 || trim >= 0.5) {
    stop("trim must be in [0, 0.5).", call. = FALSE)
  }
  if (!is.numeric(ucell_pseudocount) || length(ucell_pseudocount) != 1L ||
      !is.finite(ucell_pseudocount) || ucell_pseudocount <= 0) {
    stop("ucell_pseudocount must be one finite positive number.", call. = FALSE)
  }
  for (item in list(
    expression_response_threshold = expression_response_threshold,
    gene_change_tolerance = gene_change_tolerance
  )) {
    if (!is.numeric(item) || length(item) != 1L || !is.finite(item) || item < 0) {
      stop("Response thresholds and gene_change_tolerance must be finite and non-negative.",
           call. = FALSE)
    }
  }

  # UPDATED: choose a scale-appropriate UCell default when the user does not
  # explicitly supply one. Difference and log-ratio values are not comparable.
  if (is.null(ucell_response_threshold)) {
    ucell_response_threshold <- if (ucell_change == "difference") 0.005 else 0.1
  }
  if (!is.numeric(ucell_response_threshold) ||
      length(ucell_response_threshold) != 1L ||
      !is.finite(ucell_response_threshold) || ucell_response_threshold < 0) {
    stop("ucell_response_threshold must be NULL or one finite non-negative number.",
         call. = FALSE)
  }

  logical_args <- list(overwrite_outputs = overwrite_outputs,
                       return_details = return_details,
                       verbose = verbose)
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

  assays <- Seurat::Assays(seurat_obj)
  missing_assays <- setdiff(c(baseline_assay, perturbation_name), assays)
  if (length(missing_assays) > 0L) {
    stop("Missing assay(s): ", paste(missing_assays, collapse = ", "), ".",
         call. = FALSE)
  }

  score_methods <- unique(as.character(score_methods))
  bad_methods <- setdiff(score_methods, c("expression", "ucell"))
  if (length(score_methods) == 0L || length(bad_methods) > 0L) {
    stop("score_methods must contain 'expression', 'ucell', or both.",
         call. = FALSE)
  }
  if ("ucell" %in% score_methods && !requireNamespace("UCell", quietly = TRUE)) {
    warning(
      "UCell is not installed; continuing with score_methods = 'expression'.",
      call. = FALSE
    )
    score_methods <- setdiff(score_methods, "ucell")
  }
  if (length(score_methods) == 0L) {
    stop("No available score method remains.", call. = FALSE)
  }

  metadata <- seurat_obj@meta.data
  if (!is.null(group.by) && !group.by %in% colnames(metadata)) {
    stop("group.by column not found: ", group.by, call. = FALSE)
  }
  if (!is.null(sample.by) && !sample.by %in% colnames(metadata)) {
    stop("sample.by column not found: ", sample.by, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 2. Read and align paired assay matrices.
  # ---------------------------------------------------------------------------
  is_v5 <- inherits(seurat_obj[[baseline_assay]], "Assay5")
  get_layer <- function(assay, layer_name) {
    out <- if (is_v5) {
      Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer_name)
    } else {
      Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer_name)
    }
    if (nrow(out) == 0L || ncol(out) == 0L) {
      stop("Assay '", assay, "' has an empty '", layer_name, "' layer.",
           call. = FALSE)
    }
    out
  }

  say("Loading paired baseline and perturbation data.")
  baseline <- get_layer(baseline_assay, layer)
  perturbed <- get_layer(perturbation_name, layer)
  features <- rownames(baseline)
  cells <- colnames(baseline)

  if (!setequal(cells, colnames(perturbed))) {
    stop("Baseline and perturbation assays contain different cells.",
         call. = FALSE)
  }
  missing_features <- setdiff(features, rownames(perturbed))
  if (length(missing_features) > 0L) {
    stop(
      "Perturbation assay is missing ", length(missing_features),
      " baseline feature(s).",
      call. = FALSE
    )
  }
  # Enforce identical row/column order before paired subtraction.
  perturbed <- perturbed[features, cells, drop = FALSE]
  metadata <- metadata[cells, , drop = FALSE]

  # ---------------------------------------------------------------------------
  # 3. Filter gene sets while retaining coverage information.
  # ---------------------------------------------------------------------------
  gene_sets_original_n <- vapply(
    gene_sets,
    function(x) length(unique(as.character(x))),
    integer(1L)
  )
  gene_sets_filtered <- lapply(gene_sets, function(x) {
    intersect(unique(as.character(x)), features)
  })
  gene_sets_filtered <- gene_sets_filtered[
    vapply(gene_sets_filtered, length, integer(1L)) >= as.integer(min_genes)
  ]
  if (length(gene_sets_filtered) == 0L) {
    stop("No gene sets retained after feature filtering.", call. = FALSE)
  }

  # Define the biological grouping and replicate labels used downstream.
  group_vec <- if (is.null(group.by)) {
    rep("all_cells", length(cells))
  } else {
    as.character(metadata[[group.by]])
  }
  if (anyNA(group_vec)) group_vec[is.na(group_vec)] <- "NA"

  sample_vec <- if (is.null(sample.by)) {
    rep(NA_character_, length(cells))
  } else {
    as.character(metadata[[sample.by]])
  }
  if (!is.null(sample.by) && anyNA(sample_vec)) {
    stop("sample.by contains missing values.", call. = FALSE)
  }

  # `auto` deliberately avoids treating cells as biological replicates.
  if (null_unit == "auto") {
    null_unit <- if (is.null(sample.by)) "none" else "sample"
  }
  if (null_unit == "sample" && is.null(sample.by)) {
    stop("null_unit = 'sample' requires sample.by.", call. = FALSE)
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
    stop("n_permutations must be at least 100 when null calibration is enabled.",
         call. = FALSE)
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

  # Gene-level post-minus-pre values support expression scoring and leading
  # gene descriptions. They are not the UCell response themselves.
  gene_delta <- perturbed - baseline

  # ---------------------------------------------------------------------------
  # 4. Produce one paired score per cell × pathway × scoring method.
  # ---------------------------------------------------------------------------
  pathway_scores <- list()
  if ("expression" %in% score_methods) {
    say("Computing paired expression-sensitive pathway deltas.")
    pathway_scores$expression <- lapply(gene_sets_filtered, function(genes) {
      # colMeans averages pathway genes separately in every cell.
      pre <- Matrix::colMeans(baseline[genes, , drop = FALSE])
      post <- Matrix::colMeans(perturbed[genes, , drop = FALSE])
      list(
        pre = as.numeric(pre),
        post = as.numeric(post),
        response = as.numeric(post - pre),
        change_measure = "difference"
      )
    })
  }

  if ("ucell" %in% score_methods) {
    say("Computing paired UCell pathway responses using '", ucell_change, "'.")
    original_meta <- seurat_obj@meta.data
    suffix_pre <- paste0("__AlertPre", abs(as.integer(seed)))
    suffix_post <- paste0("__AlertPost", abs(as.integer(seed)))

    # UCell temporarily appends score columns to meta.data.
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

    missing_scores <- setdiff(
      names(gene_sets_filtered),
      intersect(pre_names, post_names)
    )
    if (length(missing_scores) > 0L) {
      seurat_obj@meta.data <- original_meta
      stop(
        "UCell did not return matched scores for: ",
        paste(missing_scores, collapse = ", "), ".",
        call. = FALSE
      )
    }

    pathway_scores$ucell <- lapply(names(gene_sets_filtered), function(pathway) {
      pre <- as.numeric(scored_meta[[pre_cols[[pathway]]]])
      post <- as.numeric(scored_meta[[post_cols[[pathway]]]])

      # UPDATED: retain the stable paired difference and restore the legacy
      # log2-ratio as an explicit, reproducible option.
      response <- if (ucell_change == "difference") {
        post - pre
      } else {
        log2((post + ucell_pseudocount) / (pre + ucell_pseudocount))
      }
      list(
        pre = pre,
        post = post,
        response = response,
        change_measure = ucell_change
      )
    })
    names(pathway_scores$ucell) <- names(gene_sets_filtered)

    # Do not permanently bloat the Seurat metadata with temporary UCell columns.
    seurat_obj@meta.data <- original_meta
  }

  # ---------------------------------------------------------------------------
  # 5. Describe direct versus propagation-only pathway overlap.
  # ---------------------------------------------------------------------------
  # These fields are annotations only; they do not filter or weight the scores.
  direct_genes <- character(0)
  propagated_genes <- character(0)
  counts_ok <- TRUE
  base_counts <- pert_counts <- NULL
  tryCatch({
    base_counts <- get_layer(baseline_assay, "counts")
    pert_counts <- get_layer(perturbation_name, "counts")
    if (!setequal(colnames(base_counts), colnames(pert_counts)) ||
        !all(rownames(base_counts) %in% rownames(pert_counts))) {
      stop("Counts layers cannot be aligned.")
    }
    pert_counts <- pert_counts[
      rownames(base_counts), colnames(base_counts), drop = FALSE
    ]
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

  # ---------------------------------------------------------------------------
  # 6. Store cell-level results in one consistent long-format table.
  # ---------------------------------------------------------------------------
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
        change_measure = score$change_measure,
        pre_score = score$pre,
        post_score = score$post,
        delta = score$response,
        stringsAsFactors = FALSE
      )
    }
  }
  per_cell <- do.call(rbind, per_cell_parts)
  rownames(per_cell) <- NULL

  # ---------------------------------------------------------------------------
  # 7. Summarize effect, prevalence, direction, and optional null calibration.
  # ---------------------------------------------------------------------------
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
      score_entry <- pathway_scores[[method]][[pathway]]
      values_all <- score_entry$response

      for (group_name in unique(group_vec)) {
        idx <- which(group_vec == group_name & is.finite(values_all))
        vals <- values_all[idx]
        if (length(vals) == 0L) next

        effect <- summarize_values(vals)
        effect_direction <- sign(effect)
        responder <- abs(vals) >= threshold
        response_rate <- mean(responder)
        tie_fraction <- mean(!responder)

        # Directional consistency uses responding cells as its denominator;
        # robust_fraction uses all cells as its denominator.
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

        # Summarize cells within biological samples before inference.
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

        # Two-sided sign-flip calibration. Sample-level units are preferable;
        # cell-level calibration is deliberately labeled descriptive.
        empirical_p <- null_z <- null_center <- null_scale <- NA_real_
        calibration_effect <- NA_real_
        if (null_unit != "none") {
          units <- if (null_unit == "sample") sample_effects else vals
          units <- units[is.finite(units)]
          if (length(units) >= 2L) {
            calibration_effect <- summarize_values(units)
            null_effects <- replicate(as.integer(n_permutations), {
              randomized_sign <- sample(
                c(-1, 1), length(units), replace = TRUE
              )
              summarize_values(units * randomized_sign)
            })
            null_center <- stats::median(null_effects, na.rm = TRUE)
            null_scale <- stats::mad(
              null_effects, center = null_center, na.rm = TRUE
            )
            if (is.finite(null_scale) && null_scale > 0) {
              null_z <- (calibration_effect - null_center) / null_scale
            }
            empirical_p <- (
              1 + sum(abs(null_effects) >= abs(calibration_effect), na.rm = TRUE)
            ) / (1 + sum(is.finite(null_effects)))
          }
        }

        # Leading genes always describe data-layer expression changes, including
        # when the pathway-level score_method is UCell.
        genes <- gene_sets_filtered[[pathway]]
        gene_effect <- as.numeric(
          Matrix::rowMeans(gene_delta[genes, idx, drop = FALSE])
        )
        names(gene_effect) <- genes
        lead_order <- head(
          order(abs(gene_effect), decreasing = TRUE),
          min(as.integer(leading_n), length(gene_effect))
        )
        leading_genes <- paste(
          paste0(
            names(gene_effect)[lead_order], "=",
            signif(gene_effect[lead_order], 3)
          ),
          collapse = ";"
        )
        gene_direction_fraction <- if (effect_direction == 0) {
          NA_real_
        } else {
          mean(sign(gene_effect) == effect_direction)
        }

        # Priorities are ranking aids, not inferential statistics.
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
          change_measure = score_entry$change_measure,
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

  # Benjamini-Hochberg adjustment is performed separately for each group and
  # scoring method, matching the family of pathways compared together.
  metrics$fdr <- NA_real_
  calibration_groups <- interaction(
    metrics$group, metrics$score_method, drop = TRUE
  )
  for (calibration_group in levels(calibration_groups)) {
    idx <- which(
      calibration_groups == calibration_group &
        is.finite(metrics$empirical_p)
    )
    if (length(idx) > 0L) {
      metrics$fdr[idx] <- stats::p.adjust(
        metrics$empirical_p[idx], method = "BH"
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 8. Record a reproducible, versioned state in @misc.
  # ---------------------------------------------------------------------------
  safe_package_version <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    tryCatch(
      as.character(utils::packageVersion(package)),
      error = function(e) NA_character_
    )
  }

  params <- list(
    schema_version = 2L,
    perturbation_name = perturbation_name,
    baseline_assay = baseline_assay,
    group.by = group.by,
    sample.by = sample.by,
    layer = layer,
    min_genes = as.integer(min_genes),
    score_methods = names(pathway_scores),
    ucell_change = ucell_change,
    ucell_pseudocount = if (ucell_change == "log2_ratio") {
      ucell_pseudocount
    } else {
      NA_real_
    },
    summary_fun = summary_fun,
    trim = trim,
    expression_response_threshold = expression_response_threshold,
    ucell_response_threshold = ucell_response_threshold,
    null_unit = null_unit,
    n_permutations = if (null_unit == "none") 0L else as.integer(n_permutations),
    seed = seed,
    leading_n = as.integer(leading_n),
    gene_change_tolerance = gene_change_tolerance,
    compact_version = safe_package_version("compact"),
    seurat_version = safe_package_version("Seurat"),
    ucell_version = if ("ucell" %in% names(pathway_scores)) {
      safe_package_version("UCell")
    } else {
      NA_character_
    },
    timestamp = Sys.time()
  )

  seurat_obj@misc[[misc_name]] <- list(
    schema_version = 2L,
    per_cell = per_cell,
    metrics = metrics,
    gene_sets = gene_sets_filtered,
    params = params
  )

  # ---------------------------------------------------------------------------
  # 9. UPDATED: optional reproducible output using AlertSystemHelper.R.
  # ---------------------------------------------------------------------------
  info_file <- NULL
  metrics_file <- NULL
  if (!is.null(out_dir)) {
    if (!exists("SaveAlertSystem", mode = "function")) {
      stop(
        "out_dir was supplied, but SaveAlertSystem() is unavailable. ",
        "Source AlertSystemHelper.R before calling AlertSystemScore().",
        call. = FALSE
      )
    }
    if (!is.character(out_dir) || length(out_dir) != 1L ||
        is.na(out_dir) || !nzchar(out_dir)) {
      stop("out_dir must be NULL or one non-empty path.", call. = FALSE)
    }
    if (is.null(output_prefix)) {
      safe_perturbation <- gsub("[^A-Za-z0-9]+", "_", perturbation_name)
      output_prefix <- paste0("AlertSystem_", safe_perturbation)
    }
    if (!is.character(output_prefix) || length(output_prefix) != 1L ||
        is.na(output_prefix) || !nzchar(output_prefix) ||
        grepl("[/\\\\]", output_prefix)) {
      stop("output_prefix must be one non-empty filename prefix without path separators.",
           call. = FALSE)
    }

    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(out_dir)) {
      stop("Could not create out_dir: ", out_dir, call. = FALSE)
    }

    info_file <- file.path(out_dir, paste0(output_prefix, "_state.rds"))
    metrics_file <- file.path(out_dir, paste0(output_prefix, "_metrics.csv"))
    existing <- c(info_file, metrics_file)[file.exists(c(info_file, metrics_file))]
    if (length(existing) > 0L && !isTRUE(overwrite_outputs)) {
      stop(
        "Output file(s) already exist; set overwrite_outputs = TRUE: ",
        paste(existing, collapse = ", "),
        call. = FALSE
      )
    }

    # Record output locations inside the state before serializing it.
    seurat_obj@misc[[misc_name]]$output_files <- list(
      state = info_file,
      metrics = metrics_file
    )
    SaveAlertSystem(
      seurat_obj,
      file = info_file,
      misc_name = misc_name,
      overwrite = isTRUE(overwrite_outputs)
    )
    utils::write.csv(metrics, file = metrics_file, row.names = FALSE)
    say("Saved Alert System state: ", info_file)
    say("Saved Alert System metrics: ", metrics_file)
  }

  say(
    "Completed: ", length(gene_sets_filtered), " pathways, ",
    length(unique(group_vec)), " group(s), methods = ",
    paste(names(pathway_scores), collapse = ", "), "."
  )

  # ---------------------------------------------------------------------------
  # 10. Return either the updated object or an inspection-friendly result list.
  # ---------------------------------------------------------------------------
  if (!isTRUE(return_details)) {
    return(seurat_obj)
  }
  list(
    seurat = seurat_obj,
    per_cell = per_cell,
    metrics = metrics,
    gene_sets = gene_sets_filtered,
    params = params,
    info_file = info_file,
    metrics_file = metrics_file
  )
}
