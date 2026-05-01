# ============================================================
# COMPACT: Combine multiple perturbation assays (delta-additive)
# COMPACT helpers: provenance + baseline sanity-check + combine
# ============================================================

#' Null-coalescing operator
#'
#' Returns the left-hand side if it is not \code{NULL}, otherwise returns the
#' right-hand side. This is useful for safely accessing optional list elements
#' or providing defaults without repeated \code{is.null()} checks.
#'
#' @param a Left-hand side value.
#' @param b Right-hand side fallback value.
#'
#' @return \code{a} if not \code{NULL}, otherwise \code{b}.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' .create_msg: internal message helper factory
#'
#' Internal helper that returns a message-emitting function controlled by a
#' logical \code{verbose} flag. When \code{verbose = TRUE}, the returned function
#' forwards its arguments to \code{message()}; otherwise it is silent.
#'
#' @param verbose Logical; whether messages should be emitted.
#'
#' @return A function that conditionally emits messages.
#'
#' @keywords internal
#' @noRd
.create_msg <- function(verbose) {
  force(verbose)
  function(...) {
    if (isTRUE(verbose)) message(...)
  }
}


# ------------------------------------------------------------
# 1) Optional: stamp baseline provenance onto a perturbation assay
# ------------------------------------------------------------
#' StampPerturbProvenance: record baseline provenance for a perturb assay
#'
#' This function stores lightweight provenance metadata on a perturbation
#' assay (e.g. created by \code{ModulePerturbation()}) indicating which
#' baseline assay and layer(s) it was derived from.
#'
#' This enables \code{CombinePerturbAssays()} to perform a strict, deterministic
#' check that all perturbation assays being combined were generated from the
#' same baseline representation (e.g. baseline assay \code{RNA} vs \code{SCT},
#' and baseline layer \code{counts}).
#'
#' Provenance is stored in \code{seurat_obj[[perturb_assay]]@misc$COMPACT} with fields:
#' \itemize{
#'   \item \code{baseline_assay}
#'   \item \code{baseline_layers} (e.g. \code{"counts"})
#'   \item \code{method}
#'   \item \code{stamped_at}
#' }
#'
#' @param seurat_obj A Seurat object.
#' @param perturb_assay Character; name of the perturbation assay to stamp.
#' @param baseline_assay Character; baseline assay name used to generate the perturbation.
#' @param method Character; description of how the assay was produced.
#' @param overwrite Logical; overwrite existing provenance fields.
#'
#' @return Seurat object with provenance written to the assay misc.
#' @export
StampPerturbProvenance <- function(
  seurat_obj,
  perturb_assay,
  baseline_assay = "RNA",
  method         = "ModulePerturbation",
  overwrite      = TRUE
) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  stopifnot(perturb_assay %in% Seurat::Assays(seurat_obj))
  stopifnot(baseline_assay %in% Seurat::Assays(seurat_obj))

  if (is.null(seurat_obj[[perturb_assay]]@misc)) seurat_obj[[perturb_assay]]@misc <- list()
  if (is.null(seurat_obj[[perturb_assay]]@misc$COMPACT)) seurat_obj[[perturb_assay]]@misc$COMPACT <- list()

  if (!overwrite) {
    if (!is.null(seurat_obj[[perturb_assay]]@misc$COMPACT$baseline_assay) ||
        !is.null(seurat_obj[[perturb_assay]]@misc$COMPACT$baseline_layers)) {
      return(seurat_obj)
    }
  }

  seurat_obj[[perturb_assay]]@misc$COMPACT$baseline_assay  <- baseline_assay
  # seurat_obj[[perturb_assay]]@misc$COMPACT$baseline_layers <- c("counts", "data")
  seurat_obj[[perturb_assay]]@misc$COMPACT$baseline_layers <- c("counts")
  seurat_obj[[perturb_assay]]@misc$COMPACT$method          <- method
  seurat_obj[[perturb_assay]]@misc$COMPACT$stamped_at      <- as.character(Sys.time())

  seurat_obj
}




# ------------------------------------------------------------
# 2) Optional heuristic: check a perturbation assay is consistent with a baseline
#    (useful when provenance is missing; NOT a proof)
# ------------------------------------------------------------
#' .CheckPerturbBaselineConsistency: heuristic sanity-check against a baseline
#'
#' This function performs a heuristic check that a perturbation expression
#' matrix is numerically consistent with a specified baseline matrix. It is
#' intended for cases where explicit provenance is missing (i.e. no call to
#' \code{StampPerturbProvenance()}) and the user wants a guardrail against
#' combining assays generated from mismatched baselines/normalizations/slots.
#'
#' The check computes \eqn{\Delta = P - B} on overlapping genes and summarizes:
#' \itemize{
#'   \item fraction of entries with \code{|Δ| > eps}
#'   \item \code{median(|Δ|)}
#'   \item \code{|mean(Δ)|} (detects global shifts)
#' }
#' If any summary exceeds user-defined thresholds, the function warns or stops.
#'
#' Note: This is not a proof of shared provenance; it is a numeric sanity check.
#'
#' @param base_mat Baseline expression matrix (genes x cells).
#' @param pert_mat Perturbation expression matrix (genes x cells).
#' @param assay_name Character label used in messages.
#' @param eps Numeric; values with \code{|Δ| <= eps} are treated as unchanged.
#' @param max_frac_changed Numeric in (0,1); warn/stop if too many entries differ.
#' @param max_abs_median Numeric; warn/stop if \code{median(|Δ|)} is too large.
#' @param max_abs_mean Numeric; warn/stop if \code{|mean(Δ)|} is too large.
#' @param strict Logical; if \code{TRUE}, stop on suspicious results; otherwise warn.
#' @param n_genes Integer; number of genes to sample for the check (default 2000).
#' @param n_cells Integer; number of cells to sample for the check (default 500).
#' @param seed Optional integer; set for reproducible sampling.
#'
#' @return Invisibly returns a list of summary statistics and a suspicious flag.
#'
#' @keywords internal
.CheckPerturbBaselineConsistency <- function(
  base_mat,
  pert_mat,
  assay_name = "pert",
  eps = 1e-3,
  max_frac_changed = 0.35,
  max_abs_median = 0.05,
  max_abs_mean   = 0.05,
  strict = FALSE,
  n_genes = 2000,
  n_cells = 500,
  seed = NULL
) {
  if (!identical(colnames(base_mat), colnames(pert_mat))) {
    err_msg <- paste0("[", assay_name, "] Cell columns differ from baseline.")
    if (strict) stop(err_msg) else warning(err_msg, call. = FALSE)
    return(invisible(NULL))
  }

  feats_all <- intersect(rownames(base_mat), rownames(pert_mat))
  if (length(feats_all) < 200) {
    err_msg <- paste0("[", assay_name, "] Too few overlapping genes with baseline (", length(feats_all), ").")
    if (strict) stop(err_msg) else warning(err_msg, call. = FALSE)
    return(invisible(NULL))
  }

  if (!is.null(seed)) set.seed(seed)
  feats <- if (length(feats_all) > n_genes) sample(feats_all, n_genes) else feats_all
  cells_all <- colnames(base_mat)
  cells <- if (length(cells_all) > n_cells) sample(cells_all, n_cells) else cells_all

  B <- base_mat[feats, cells, drop = FALSE]
  P <- pert_mat[feats, cells, drop = FALSE]
  D <- P - B

  d_vec <- as.numeric(D)
  frac_changed <- mean(abs(d_vec) > eps, na.rm = TRUE)
  abs_median   <- median(abs(d_vec), na.rm = TRUE)
  abs_mean     <- abs(mean(d_vec, na.rm = TRUE))

  suspicious <- (frac_changed > max_frac_changed) ||
    (abs_median > max_abs_median) ||
    (abs_mean   > max_abs_mean)

  summary_msg <- paste0(
    "[", assay_name, "] baseline check (sampled ",
    length(feats), " genes x ", length(cells), " cells): ",
    "frac(|Δ|>eps)=", signif(frac_changed, 3),
    ", median(|Δ|)=", signif(abs_median, 3),
    ", |mean(Δ)|=", signif(abs_mean, 3),
    " (eps=", eps, ")"
  )

  if (suspicious) {
    warn_msg <- paste0(
      summary_msg,
      "\n  -> Suspicious vs baseline. Possible wrong baseline assay/layer, different normalization, or different object."
    )
    if (strict) stop(warn_msg) else warning(warn_msg, call. = FALSE)
  } else {
    message(summary_msg)
  }

  invisible(list(
    frac_changed = frac_changed,
    abs_median   = abs_median,
    abs_mean     = abs_mean,
    suspicious   = suspicious,
    n_genes_used = length(feats),
    n_cells_used = length(cells)
  ))
}



# ------------------------------------------------------------
# 3) Combine perturbation assays via delta superposition:
#    P_combined = B + sum_i (P_i - B)
#
# Provenance logic:
# - If provenance exists for all assays: enforce match to requested baseline (stop/warn)
# - If provenance missing: optional heuristic checks; always emit a big warning banner
# ------------------------------------------------------------
#' CombinePerturbAssays: create a composite perturb assay by delta superposition (sparse-safe; COUNTS-only)
#'
#' This function combines multiple perturbation assays into a single composite
#' assay using additive superposition of perturbation deltas relative to a shared
#' baseline on the \code{counts} layer:
#' \deqn{P_{\mathrm{combined}} = B + \sum_i (P_i - B)}
#'
#' Implementation note: this function avoids creating a large dense expanded
#' matrix. It combines perturbations strictly on the baseline feature space
#' (rows = \code{rownames(baseline)}). Features present in a perturb assay but
#' absent from baseline are dropped with a warning.
#'
#' Baseline safety logic:
#' \itemize{
#'   \item If all input assays have provenance (via \code{StampPerturbProvenance()}),
#'         the function verifies recorded \code{baseline_assay} and
#'         \code{baseline_layers} include \code{"counts"}.
#'   \item If provenance is missing for any assay, the function emits a warning banner.
#'         Optionally, it runs heuristic numeric checks on the \code{counts} layer
#'         via \code{.CheckPerturbBaselineConsistency()} for the unstamped assays
#'         (not a proof; may be less informative than checks on normalized data).
#' }
#'
#' Data layer handling:
#' After combining \code{counts}, the function generates the composite assay's
#' \code{data} layer by calling the package-internal \code{log_normalize()} helper
#' (the same one used by \code{ModulePerturbation()}) with \code{col_sums} taken
#' from the BASELINE counts (\code{colSums(base_counts)}), not from the composite
#' counts. This preserves the per-cell denominator across the baseline assay and
#' all perturb / composite assays, so their \code{data} layers remain directly
#' comparable. \code{Seurat::NormalizeData()} is intentionally NOT used here
#' because it would recompute library sizes from the composite counts, breaking
#' that comparability.
#'
#' Combine by delta-superposition (counts only):
#' \itemize{
#'   \item \code{comp_counts = base_counts + \sum_i (pert_counts_i - base_counts)}
#' }
#'
#' Requirements:
#' \itemize{
#'   \item baseline assay must have layer \code{counts}
#'   \item each perturb assay must have layer \code{counts}
#'   \item all assays must share identical cell columns
#' }
#'
#' @param seurat_obj A Seurat object containing baseline and perturbation assays.
#' @param perturb_assays Character vector of perturbation assay names to combine.
#' @param baseline_assay Character; baseline assay name (default \code{"RNA"}).
#' @param new_assay Character; name of the new composite assay.
#' @param require_provenance Logical; if TRUE, missing provenance triggers stop/warn.
#' @param strict_provenance_check Logical; if TRUE, stop on provenance mismatch; else warn.
#' @param run_consistency_check_if_missing Logical; if TRUE, run heuristic checks when provenance missing.
#' @param strict_consistency_check Logical; if TRUE, stop on suspicious heuristic results; else warn.
#' @param eps,max_frac_changed,max_abs_median,max_abs_mean Thresholds for heuristic checks.
#' @param allow_unverified Logical; if provenance missing and heuristic checks disabled, proceed only if TRUE.
#' @param allow_dense Logical; if FALSE (default), stop when any requested layer
#'   returns a dense matrix (to avoid RAM blow-ups).
#' @param normalize_method Deprecated/ignored. The composite \code{data} layer is
#'   always produced via \code{log_normalize()} (LogNormalize / \code{log1p})
#'   to match \code{ModulePerturbation()}. Kept only for backward compatibility.
#' @param scale.factor Scale factor passed to the internal \code{log_normalize()}
#'   helper (default \code{1e4}). Should match the value used when the input
#'   perturb assays were generated.
#'
#' @return A Seurat object with a new composite perturbation assay.
#' @export
#'
#' @examples
#' \dontrun{
#' # After creating perturb assays:
#' seurat_obj_perturb <- ModulePerturbation(
#'   seurat_obj_perturb,
#'   mod = "CBh-M6",
#'   perturb_dir = 5,
#'   perturbation_name = "M6_up",
#'   graph = "RNA_nn",
#'   n_hubs = 10
#' )
#'
#' seurat_obj_perturb <- ModulePerturbation(
#'   seurat_obj_perturb,
#'   mod = "CBh-M2",
#'   perturb_dir = -2,
#'   perturbation_name = "M2_down",
#'   graph = "RNA_nn",
#'   n_hubs = 10
#' )
#'
#' # Optional: stamp provenance (recommended)
#' seurat_obj_perturb <- StampPerturbProvenance(
#'   seurat_obj_perturb,
#'   perturb_assay = "M6_up",
#'   baseline_assay = "RNA"
#' )
#'
#' seurat_obj_perturb <- StampPerturbProvenance(
#'   seurat_obj_perturb,
#'   perturb_assay = "M2_down",
#'   baseline_assay = "RNA"
#' )
#'
#' seurat_obj_perturb <- CombinePerturbAssays(
#'   seurat_obj_perturb,
#'   perturb_assays = c("M6_up", "M2_down"),
#'   baseline_assay = "RNA",
#'   new_assay = "M6up_M2down",
#'   require_provenance = FALSE,
#'   run_consistency_check_if_missing = FALSE,
#'   allow_unverified = TRUE
#' )
#' }

CombinePerturbAssays <- function(
  seurat_obj,
  perturb_assays,
  baseline_assay = "RNA",
  new_assay = paste(perturb_assays, collapse = "_"),

  require_provenance = FALSE,
  strict_provenance_check = TRUE,

  run_consistency_check_if_missing = TRUE,
  strict_consistency_check = FALSE,

  eps = 1e-3,
  max_frac_changed = 0.35,
  max_abs_median = 0.05,
  max_abs_mean = 0.05,

  allow_unverified = TRUE,
  allow_dense = FALSE,

  # NEW: normalization settings for data layer
  normalize_method = "LogNormalize",
  scale.factor = 1e4,

  verbose = TRUE
) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  stopifnot(baseline_assay %in% Seurat::Assays(seurat_obj))
  stopifnot(all(perturb_assays %in% Seurat::Assays(seurat_obj)))

  msg_prefix <- "COMPACT:Combine"
  is_v5 <- hdWGCNA::CheckSeurat5()

  msg <- .create_msg(verbose)
  .p <- function(...) msg(paste0("[", msg_prefix, "] ", ...))

  .p("Starting (COUNTS-only). baseline_assay='", baseline_assay, "', new_assay='", new_assay, "'.")
  .p("Perturb assays: ", paste(perturb_assays, collapse = ", "))
  .p("Seurat v5 detected: ", isTRUE(is_v5), ".")

  as_sparse_or_stop <- function(m, label) {
    if (inherits(m, "sparseMatrix")) return(m)
    if (!allow_dense) {
      stop(
        "CombinePerturbAssays(): returned a dense matrix for '", label, "'. ",
        "To avoid RAM blow-ups, keep layers sparse or set allow_dense=TRUE."
      )
    }
    Matrix::Matrix(m, sparse = TRUE)
  }

  has_layer <- function(assay, layer) {
    if (is_v5) {
      layer %in% SeuratObject::Layers(seurat_obj[[assay]])
    } else {
      layer %in% c("counts", "data", "scale.data")
    }
  }

  get_layer <- function(assay, layer) {
    if (!has_layer(assay, layer)) {
      stop("CombinePerturbAssays(): assay '", assay, "' is missing required layer '", layer, "'.")
    }
    if (is_v5) {
      m <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
    } else {
      m <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
    }
    as_sparse_or_stop(m, paste0(assay, ":", layer))
  }

  # ---- hard requirement: counts only ----
  .p("Checking required layer exists: counts (baseline + all perturb assays)...")
  for (a in c(baseline_assay, perturb_assays)) {
    if (!has_layer(a, "counts")) {
      stop("CombinePerturbAssays(): assay '", a, "' is missing required layer 'counts'.")
    }
  }

  # ---- baseline counts ----
  .p("Loading baseline COUNTS from assay '", baseline_assay, "'...")
  base_counts <- get_layer(baseline_assay, "counts")
  cells <- colnames(base_counts)
  .p("Baseline counts loaded: ", nrow(base_counts), " features x ", length(cells), " cells.")

  # ---- perturb counts ----
  .p("Loading perturb COUNTS for ", length(perturb_assays), " assays...")
  pert_counts <- lapply(perturb_assays, function(a) {
    m <- get_layer(a, "counts")
    if (!identical(colnames(m), cells)) stop("Cell columns differ across assays (", a, ") in counts.")
    m
  })
  names(pert_counts) <- perturb_assays
  .p("Perturb counts loaded.")

  # ---- provenance inspection (counts-only) ----
  get_prov <- function(a) seurat_obj[[a]]@misc$COMPACT %||% NULL
  prov_list <- lapply(perturb_assays, get_prov)
  has_prov  <- vapply(prov_list, function(x) !is.null(x), logical(1))
  all_have_prov <- all(has_prov)

  .p("Checking provenance stamps (assay@misc$COMPACT)...")

  if (!all_have_prov) {
    missing <- perturb_assays[!has_prov]
    .p("Provenance missing for: ", paste(missing, collapse = ", "))

    if (require_provenance) {
      prov_msg <- paste0(
        "COMPACT PROVENANCE ERROR: Missing assay@misc$COMPACT provenance for: ",
        paste(missing, collapse = ", "),
        ".\nEither stamp provenance (StampPerturbProvenance) or set require_provenance=FALSE."
      )
      if (strict_provenance_check) stop(prov_msg) else warning(prov_msg, call. = FALSE)
    }

    banner <- paste0(
      "COMPACT WARNING (provenance missing): Unstamped assays: ",
      paste(missing, collapse = ", "),
      ". Cannot verify they were derived from baseline_assay='", baseline_assay, "'."
    )

    if (!run_consistency_check_if_missing) {
      if (!allow_unverified) stop(paste0(banner, "\nSet allow_unverified=TRUE to proceed without checks."))
      warning(paste0(banner, "\nProceeding WITHOUT heuristic checks (researcher must verify)."), call. = FALSE)
    } else {
      warning(paste0(
        banner,
        "\nRunning heuristic baseline-consistency checks on COUNTS (not a proof; may be less informative than data)."
      ), call. = FALSE)

      .p("Running heuristic baseline-consistency checks on COUNTS for: ", paste(missing, collapse = ", "))
      for (a in missing) {
        .CheckPerturbBaselineConsistency(
          base_mat = base_counts,
          pert_mat = pert_counts[[a]],
          assay_name = a,
          eps = eps,
          max_frac_changed = max_frac_changed,
          max_abs_median = max_abs_median,
          max_abs_mean = max_abs_mean,
          strict = strict_consistency_check
        )
      }
    }
  } else {
    .p("All assays have provenance. Verifying baseline_assay='", baseline_assay, "' and baseline_layers include 'counts'...")

    prov_base   <- vapply(prov_list, function(x) x$baseline_assay %||% NA_character_, character(1))
    prov_layers <- lapply(prov_list, function(x) x$baseline_layers %||% character(0))

    base_mismatch <- !is.na(prov_base) & (prov_base != baseline_assay)
    layers_bad <- vapply(prov_layers, function(v) !("counts" %in% v), logical(1))

    if (any(base_mismatch) || any(layers_bad)) {
      badA <- perturb_assays[base_mismatch | layers_bad]
      mismatch_msg <- paste0(
        "COMPACT PROVENANCE MISMATCH: These assays do not match baseline_assay='",
        baseline_assay, "' and/or are not stamped with baseline_layers including 'counts': ",
        paste(badA, collapse = ", "),
        ".\nRecorded provenance:\n",
        paste0(
          "  - ", perturb_assays,
          ": baseline_assay=", prov_base,
          ", baseline_layers=",
          vapply(prov_layers, function(v) paste(v, collapse = ","), character(1)),
          collapse = "\n"
        )
      )
      if (strict_provenance_check) stop(mismatch_msg) else warning(mismatch_msg, call. = FALSE)
    }
  }

  # ---- combine counts via delta-superposition ----
  combine_counts <- function(base_mat, pert_list) {
    base_feats <- rownames(base_mat)
    comp <- base_mat

    for (a in names(pert_list)) {
      p <- pert_list[[a]]
      common <- intersect(base_feats, rownames(p))
      if (length(common) == 0) {
        warning("CombinePerturbAssays(): no overlapping features for assay '", a, "' in counts. Skipping.", call. = FALSE)
        next
      }

      extra <- setdiff(rownames(p), base_feats)
      if (length(extra) > 0) {
        warning("CombinePerturbAssays(): assay '", a, "' has ", length(extra),
                " features not in baseline (counts). Dropping.", call. = FALSE)
      }

      delta <- p[common, , drop = FALSE] - base_mat[common, , drop = FALSE]
      comp[common, ] <- comp[common, , drop = FALSE] + delta
    }

    comp
  }

  .p("Combining COUNTS layer via delta-superposition...")
  comp_counts <- combine_counts(base_counts, pert_counts)

  # ---- create new assay ----
  .p("Writing composite assay '", new_assay, "' (counts) and normalizing data layer ",
     "via log_normalize() with baseline colSums (matching ModulePerturbation)...")

  if (is_v5) {
    seurat_obj[[new_assay]] <- SeuratObject::CreateAssay5Object(counts = comp_counts)
  } else {
    seurat_obj[[new_assay]] <- Seurat::CreateAssayObject(counts = comp_counts)
  }

  # Normalize using BASELINE col_sums to stay consistent with ModulePerturbation(),
  # which calls log_normalize(exp_simulated, colSums(exp)). Using Seurat::NormalizeData()
  # here would instead use colSums(comp_counts), giving a different denominator and
  # breaking comparability with the baseline 'data' layer and with each P_i@data.
  base_col_sums <- Matrix::colSums(base_counts)
  comp_data <- compact:::log_normalize(comp_counts, base_col_sums, scale.factor = scale.factor)

  if (is_v5) {
    seurat_obj <- SeuratObject::SetAssayData(
      seurat_obj, assay = new_assay, layer = "data", new.data = comp_data
    )
  } else {
    seurat_obj <- Seurat::SetAssayData(
      seurat_obj, assay = new_assay, slot = "data", new.data = comp_data
    )
  }

  # provenance for composite
  if (is.null(seurat_obj[[new_assay]]@misc)) seurat_obj[[new_assay]]@misc <- list()
  seurat_obj[[new_assay]]@misc$COMPACT <- list(
    baseline_assay   = baseline_assay,
    baseline_layers  = c("counts"),
    method           = "CombinePerturbAssays (counts combined; data = log_normalize w/ baseline colSums)",
    components       = perturb_assays,
    normalize_method = "LogNormalize",
    scale.factor     = scale.factor,
    created_at       = as.character(Sys.time())
  )

  .p("Done. Added assay '", new_assay,
     "'. (counts combined; data normalized against baseline colSums)")
  seurat_obj
}






# ##############
# CombinePerturbAssays <- function(
#   seurat_obj,
#   perturb_assays,
#   baseline_assay = "RNA",
#   new_assay = paste(perturb_assays, collapse = "_"),
#
#   require_provenance = FALSE,
#   strict_provenance_check = TRUE,
#
#   run_consistency_check_if_missing = TRUE,
#   strict_consistency_check = FALSE,
#
#   eps = 1e-3,
#   max_frac_changed = 0.35,
#   max_abs_median = 0.05,
#   max_abs_mean = 0.05,
#
#   allow_unverified = TRUE,
#   allow_dense = FALSE,
#   verbose = TRUE
# ) {
#   stopifnot(inherits(seurat_obj, "Seurat"))
#   stopifnot(baseline_assay %in% Seurat::Assays(seurat_obj))
#   stopifnot(all(perturb_assays %in% Seurat::Assays(seurat_obj)))
#
#   msg_prefix <- "COMPACT:Combine"
#
#   is_v5 <- hdWGCNA::CheckSeurat5()
#
#   msg <- .create_msg(verbose)
#   .p <- function(...) msg(paste0("[", msg_prefix, "] ", ...))
#
#   .p("Starting. baseline_assay='", baseline_assay, "', new_assay='", new_assay, "'.")
#   .p("Perturb assays: ", paste(perturb_assays, collapse = ", "))
#   .p("Seurat v5 detected: ", isTRUE(is_v5), ".")
#
#
#   as_sparse_or_stop <- function(m, label) {
#     if (inherits(m, "sparseMatrix")) return(m)
#     if (!allow_dense) {
#       stop(
#         "CombinePerturbAssays(): returned a dense matrix for '", label, "'. ",
#         "To avoid RAM blow-ups, keep layers sparse or set allow_dense=TRUE."
#       )
#     }
#     Matrix::Matrix(m, sparse = TRUE)
#   }
#
#   has_layer <- function(assay, layer) {
#     if (is_v5) {
#       layer %in% SeuratObject::Layers(seurat_obj[[assay]])
#     } else {
#       layer %in% c("counts", "data", "scale.data")
#     }
#   }
#
#   get_layer <- function(assay, layer) {
#     if (!has_layer(assay, layer)) {
#       stop("CombinePerturbAssays(): assay '", assay, "' is missing required layer '", layer, "'.")
#     }
#     if (is_v5) {
#       m <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
#     } else {
#       m <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
#     }
#     as_sparse_or_stop(m, paste0(assay, ":", layer))
#   }
#
#   # ---- hard requirement: baseline and all pert assays have BOTH layers ----
#   .p("Checking required layers (counts + data) exist for baseline + all perturb assays...")
#
#   required_layers <- c("counts", "data")
#   for (a in c(baseline_assay, perturb_assays)) {
#     miss <- required_layers[!vapply(required_layers, function(L) has_layer(a, L), logical(1))]
#     if (length(miss) > 0) {
#       stop(
#         "CombinePerturbAssays(): assay '", a, "' is missing required layer(s): ",
#         paste(miss, collapse = ", "),
#         ". This combine requires BOTH counts+data for baseline and each perturb assay."
#       )
#     }
#   }
#
#   # ---- baseline ----
#   .p("Loading baseline layers (counts, data) from assay '", baseline_assay, "'...")
#
#   base_counts <- get_layer(baseline_assay, "counts")
#   base_data   <- get_layer(baseline_assay, "data")
#
#   if (!identical(colnames(base_counts), colnames(base_data))) {
#     stop("CombinePerturbAssays(): baseline counts/data have different cell columns (unexpected).")
#   }
#   cells <- colnames(base_data)
#
#   .p("Baseline loaded: ", nrow(base_data), " features x ", length(cells), " cells.")
#
#   # ---- perturb layers ----
#   .p("Loading perturb layers (counts + data) for ", length(perturb_assays), " assays...")
#
#   pert_counts <- lapply(perturb_assays, function(a) {
#     m <- get_layer(a, "counts")
#     if (!identical(colnames(m), cells)) stop("Cell columns differ across assays (", a, ") in counts.")
#     m
#   })
#   names(pert_counts) <- perturb_assays
#
#   pert_data <- lapply(perturb_assays, function(a) {
#     m <- get_layer(a, "data")
#     if (!identical(colnames(m), cells)) stop("Cell columns differ across assays (", a, ") in data.")
#     m
#   })
#   names(pert_data) <- perturb_assays
#
#   .p("Perturb layers loaded.")
#
#   # ---- provenance inspection (UPDATED to baseline_layers) ----
#   get_prov <- function(a) seurat_obj[[a]]@misc$COMPACT %||% NULL
#   prov_list <- lapply(perturb_assays, get_prov)
#   has_prov  <- vapply(prov_list, function(x) !is.null(x), logical(1))
#   all_have_prov <- all(has_prov)
#
#   .p("Checking provenance stamps (assay@misc$COMPACT)...")
#
#   if (!all_have_prov) {
#     missing <- perturb_assays[!has_prov]
#     .p("Provenance missing for: ", paste(missing, collapse = ", "))
#
#     if (require_provenance) {
#       prov_msg <- paste0(
#         "COMPACT PROVENANCE ERROR: Missing assay@misc$COMPACT provenance for: ",
#         paste(missing, collapse = ", "),
#         ".\nEither stamp provenance (StampPerturbProvenance) or set require_provenance=FALSE."
#       )
#       if (strict_provenance_check) stop(prov_msg) else warning(prov_msg, call. = FALSE)
#     }
#
#     banner <- paste0(
#       "COMPACT WARNING (provenance missing): Unstamped assays: ",
#       paste(missing, collapse = ", "),
#       ". Cannot verify they were derived from baseline_assay='", baseline_assay, "'."
#     )
#
#     if (!run_consistency_check_if_missing) {
#       if (!allow_unverified) stop(paste0(banner, "\nSet allow_unverified=TRUE to proceed without checks."))
#       warning(paste0(banner, "\nProceeding WITHOUT heuristic checks (researcher must verify)."), call. = FALSE)
#     } else {
#       warning(paste0(banner, "\nRunning heuristic baseline-consistency checks on DATA layer (not a proof)."), call. = FALSE)
#
#       .p("Running heuristic baseline-consistency checks on DATA layer for: ", paste(missing, collapse = ", "))
#
#       for (a in missing) {
#         .CheckPerturbBaselineConsistency(
#           base_mat = base_data,
#           pert_mat = pert_data[[a]],
#           assay_name = a,
#           eps = eps,
#           max_frac_changed = max_frac_changed,
#           max_abs_median = max_abs_median,
#           max_abs_mean = max_abs_mean,
#           strict = strict_consistency_check
#         )
#       }
#     }
#   } else {
#     .p("All assays have provenance. Verifying baseline_assay='", baseline_assay, "' and baseline_layers include counts,data...")
#
#     prov_base <- vapply(prov_list, function(x) x$baseline_assay %||% NA_character_, character(1))
#     prov_layers <- lapply(prov_list, function(x) x$baseline_layers %||% character(0))
#
#     base_mismatch <- !is.na(prov_base) & (prov_base != baseline_assay)
#     layers_bad <- vapply(prov_layers, function(v) !all(c("counts", "data") %in% v), logical(1))
#
#     if (any(base_mismatch) || any(layers_bad)) {
#       badA <- perturb_assays[base_mismatch | layers_bad]
#       mismatch_msg <- paste0(
#         "COMPACT PROVENANCE MISMATCH: These assays do not match baseline_assay='",
#         baseline_assay, "' and/or are not stamped with baseline_layers=c('counts','data'): ",
#         paste(badA, collapse = ", "),
#         ".\nRecorded provenance:\n",
#         paste0(
#           "  - ", perturb_assays,
#           ": baseline_assay=", prov_base,
#           ", baseline_layers=",
#           vapply(prov_layers, function(v) paste(v, collapse = ","), character(1)),
#           collapse = "\n"
#         )
#       )
#       if (strict_provenance_check) stop(mismatch_msg) else warning(mismatch_msg, call. = FALSE)
#     }
#   }
#
#   # ---- combine helper ----
#   combine_layer <- function(base_mat, pert_list, layer_label) {
#     base_feats <- rownames(base_mat)
#     comp <- base_mat
#
#     for (a in names(pert_list)) {
#       p <- pert_list[[a]]
#       common <- intersect(base_feats, rownames(p))
#       if (length(common) == 0) {
#         warning("CombinePerturbAssays(): no overlapping features for assay '", a,
#                 "' in layer '", layer_label, "'. Skipping.", call. = FALSE)
#         next
#       }
#
#       extra <- setdiff(rownames(p), base_feats)
#       if (length(extra) > 0) {
#         warning("CombinePerturbAssays(): assay '", a, "' has ", length(extra),
#                 " features not in baseline (layer '", layer_label, "'). Dropping.", call. = FALSE)
#       }
#
#       delta <- p[common, , drop = FALSE] - base_mat[common, , drop = FALSE]
#       comp[common, ] <- comp[common, , drop = FALSE] + delta
#     }
#
#     comp
#   }
#
#   .p("Combining COUNTS layer via delta-superposition...")
#
#   comp_counts <- combine_layer(base_counts, pert_counts, "counts")
#
#   .p("Combining DATA layer via delta-superposition...")
#
#   comp_data   <- combine_layer(base_data,   pert_data,   "data")
#
#   # ---- store new assay with BOTH layers ----
#   .p("Writing composite assay '", new_assay, "' (counts + data)...")
#
#   if (is_v5) {
#     seurat_obj[[new_assay]] <- SeuratObject::CreateAssay5Object(
#       counts = comp_counts,
#       data   = comp_data
#     )
#   } else {
#     seurat_obj[[new_assay]] <- Seurat::CreateAssayObject(
#       counts = comp_counts,
#       data   = comp_data
#     )
#   }
#
#   if (is.null(seurat_obj[[new_assay]]@misc)) seurat_obj[[new_assay]]@misc <- list()
#   seurat_obj[[new_assay]]@misc$COMPACT <- list(
#     baseline_assay  = baseline_assay,
#     baseline_layers = c("counts", "data"),
#     method          = "CombinePerturbAssays",
#     components      = perturb_assays,
#     created_at      = as.character(Sys.time())
#   )
#
#   .p("Done. Added assay '", new_assay, "'.")
#
#   seurat_obj
# }






#
