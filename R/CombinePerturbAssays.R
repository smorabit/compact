# ============================================================
# COMPACT: combine mixed-layer perturbation assays
# Experimental replacement for CombinePerturbAssays.R
# ============================================================

#' Combine perturbation assays created by the mixed-layer workflow
#'
#' Experimental replacement for \code{CombinePerturbAssays()}. This function is
#' designed for perturbation assays in which direct targets (hub genes or
#' user-selected genes) are changed in the \code{counts} layer and secondary
#' genes are changed by network propagation in the log-normalized \code{data}
#' layer.
#'
#' Direct and propagated effects are deliberately combined differently:
#' \enumerate{
#'   \item Count-space deltas are added to the baseline counts.
#'   \item The combined counts are log-normalized with the baseline library
#'         sizes, reproducing the primary/direct component in data space.
#'   \item For each input assay, its data-space delta is added only for genes
#'         that were not directly changed in that assay. These are the
#'         secondary effects produced by propagation.
#' }
#'
#' In symbols, for input assay \eqn{i}, baseline counts \eqn{B}, reconstructed
#' baseline log-expression \eqn{L}, perturbed counts \eqn{P_i}, perturbed data
#' \eqn{D_i}, and direct-gene set \eqn{H_i} inferred from \eqn{P_i-B}:
#' \deqn{C = max(0, B + sum_i(P_i-B))}
#' \deqn{E = max(0, logNorm(C; colSums(B)) +
#'                   sum_i 1[g not in H_i](D_i-L))}
#'
#' This component-wise exclusion is important. If a gene is a direct target in
#' one assay but a secondary target in another, its secondary effect from the
#' other assay is retained.
#'
#' The function is intentionally kept under a new name so the existing
#' \code{CombinePerturbAssays()} implementation remains available for comparison.
#'
#' @param seurat_obj A Seurat object containing the baseline and perturbation
#'   assays.
#' @param perturb_assays Character vector naming the perturbation assays.
#' @param baseline_assay Name of the assay from which the perturbations were
#'   generated. Default is \code{"RNA"}.
#' @param new_assay Name for the combined assay.
#' @param scale.factor Scale factor used by the perturbation workflow's
#'   log-normalization. Default is \code{1e4}.
#' @param direct_tolerance Non-negative numeric tolerance used to infer direct
#'   genes from count-space deltas. A gene is direct for an assay when any cell
#'   has \code{abs(perturbed_count - baseline_count) > direct_tolerance}.
#' @param check_provenance Logical; inspect \code{assay@misc$COMPACT} metadata
#'   written by \code{StampPerturbProvenance()} and verify that every available
#'   stamp names \code{baseline_assay} and includes the \code{"counts"} layer.
#'   Missing stamps warn but do not stop unless \code{require_provenance = TRUE}.
#'   Default is \code{TRUE}.
#' @param require_provenance Logical; require complete provenance metadata for
#'   every perturbation assay. Default is \code{FALSE}, allowing older or
#'   unstamped assays to be combined after a warning.
#' @param strict_provenance_check Logical; stop when an available provenance
#'   stamp conflicts with \code{baseline_assay} or does not include
#'   \code{"counts"}. If \code{FALSE}, issue a warning and continue. Default is
#'   \code{TRUE}.
#' @param floor_counts Logical; floor negative combined counts at zero. This can
#'   be needed when several knock-down deltas affect the same gene. Default is
#'   \code{TRUE}.
#' @param floor_data Logical; floor negative combined log-expression at zero,
#'   matching \code{ApplyPropagation()}. Default is \code{TRUE}.
#' @param allow_dense Logical; permit dense input matrices. Default is
#'   \code{FALSE} to protect against accidental memory expansion.
#' @param overwrite Logical; replace an existing assay named \code{new_assay}.
#'   Default is \code{FALSE}.
#' @param verbose Logical; emit progress messages.
#'
#' @return The Seurat object with the new combined assay.
#'
#' @examples
#' \dontrun{
#' object <- CombinePerturbAssays(
#'   object,
#'   perturb_assays = c("M6_up", "M2_down"),
#'   baseline_assay = "RNA",
#'   new_assay = "M6up_M2down_"
#' )
#' }
#'
#' @import Matrix
#' @export
CombinePerturbAssays <- function(
  seurat_obj,
  perturb_assays,
  baseline_assay = "RNA",
  new_assay = paste0(paste(perturb_assays, collapse = "_"), "_combined"),
  scale.factor = 1e4,
  direct_tolerance = 0,
  check_provenance = TRUE,
  require_provenance = FALSE,
  strict_provenance_check = TRUE,
  floor_counts = TRUE,
  floor_data = TRUE,
  allow_dense = FALSE,
  overwrite = FALSE,
  verbose = TRUE
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (!is.character(perturb_assays) || length(perturb_assays) < 1L ||
      anyNA(perturb_assays) || any(!nzchar(perturb_assays))) {
    stop("perturb_assays must contain at least one non-empty assay name.")
  }
  if (anyDuplicated(perturb_assays)) {
    stop("perturb_assays contains duplicate assay names.")
  }
  if (!is.character(baseline_assay) || length(baseline_assay) != 1L ||
      is.na(baseline_assay) || !nzchar(baseline_assay)) {
    stop("baseline_assay must be one non-empty assay name.")
  }
  if (!is.character(new_assay) || length(new_assay) != 1L ||
      is.na(new_assay) || !nzchar(new_assay)) {
    stop("new_assay must be one non-empty assay name.")
  }
  if (!is.numeric(scale.factor) || length(scale.factor) != 1L ||
      !is.finite(scale.factor) || scale.factor <= 0) {
    stop("scale.factor must be one finite positive number.")
  }
  if (!is.numeric(direct_tolerance) || length(direct_tolerance) != 1L ||
      !is.finite(direct_tolerance) || direct_tolerance < 0) {
    stop("direct_tolerance must be one finite non-negative number.")
  }
  logical_args <- list(
    check_provenance = check_provenance,
    require_provenance = require_provenance,
    strict_provenance_check = strict_provenance_check
  )
  invalid_logical <- vapply(
    logical_args,
    function(x) !is.logical(x) || length(x) != 1L || is.na(x),
    logical(1)
  )
  if (any(invalid_logical)) {
    stop(
      paste(names(logical_args)[invalid_logical], collapse = ", "),
      " must be a single non-missing logical value."
    )
  }
  if (isTRUE(require_provenance) && !isTRUE(check_provenance)) {
    stop("require_provenance = TRUE requires check_provenance = TRUE.")
  }

  assay_names <- Seurat::Assays(seurat_obj)
  missing_assays <- setdiff(c(baseline_assay, perturb_assays), assay_names)
  if (length(missing_assays) > 0L) {
    stop(
      "The following assays are missing: ",
      paste(missing_assays, collapse = ", "), "."
    )
  }
  if (baseline_assay %in% perturb_assays) {
    stop("baseline_assay cannot also appear in perturb_assays.")
  }
  if (new_assay %in% assay_names && !isTRUE(overwrite)) {
    stop(
      "Assay '", new_assay, "' already exists. Choose another name or set ",
      "overwrite = TRUE."
    )
  }

  say <- function(...) {
    if (isTRUE(verbose)) message("[COMPACT:Combine] ", ...)
  }

  provenance_summary <- list(
    checked = isTRUE(check_provenance),
    complete = character(0),
    missing_or_incomplete = character(0),
    mismatched = character(0)
  )

  if (isTRUE(check_provenance)) {
    say("Checking perturbation-assay provenance.")

    provenance <- lapply(perturb_assays, function(assay) {
      misc <- seurat_obj[[assay]]@misc
      if (is.null(misc) || !is.list(misc)) return(NULL)
      compact_meta <- misc$COMPACT
      if (is.null(compact_meta) || !is.list(compact_meta)) return(NULL)
      compact_meta
    })
    names(provenance) <- perturb_assays

    is_complete <- vapply(provenance, function(x) {
      if (is.null(x)) return(FALSE)
      recorded_baseline <- x$baseline_assay
      recorded_layers <- x$baseline_layers
      is.character(recorded_baseline) && length(recorded_baseline) == 1L &&
        !is.na(recorded_baseline) && nzchar(recorded_baseline) &&
        is.character(recorded_layers) && length(recorded_layers) > 0L
    }, logical(1))

    provenance_summary$complete <- perturb_assays[is_complete]
    provenance_summary$missing_or_incomplete <- perturb_assays[!is_complete]

    if (any(!is_complete)) {
      missing_names <- perturb_assays[!is_complete]
      missing_message <- paste0(
        "Provenance is missing or incomplete for perturbation assay(s): ",
        paste(missing_names, collapse = ", "),
        ". Cannot verify that they were derived from baseline_assay='",
        baseline_assay, "'. Use StampPerturbProvenance() to record it."
      )
      if (isTRUE(require_provenance)) {
        stop(missing_message)
      }
      warning(missing_message, call. = FALSE)
    }

    if (any(is_complete)) {
      complete_provenance <- provenance[is_complete]
      recorded_baselines <- vapply(
        complete_provenance, `[[`, character(1), "baseline_assay"
      )
      has_counts <- vapply(
        complete_provenance,
        function(x) "counts" %in% x$baseline_layers,
        logical(1)
      )
      mismatched <- recorded_baselines != baseline_assay | !has_counts
      provenance_summary$mismatched <- names(complete_provenance)[mismatched]

      if (any(mismatched)) {
        details <- vapply(names(complete_provenance)[mismatched], function(assay) {
          x <- complete_provenance[[assay]]
          paste0(
            assay, " (baseline_assay=", x$baseline_assay,
            ", baseline_layers=", paste(x$baseline_layers, collapse = ","), ")"
          )
        }, character(1))
        mismatch_message <- paste0(
          "Provenance mismatch for requested baseline_assay='", baseline_assay,
          "': ", paste(details, collapse = "; "),
          ". Only assays derived from the same baseline counts should be combined."
        )
        if (isTRUE(strict_provenance_check)) {
          stop(mismatch_message)
        }
        warning(mismatch_message, call. = FALSE)
      }
    }
  }

  is_v5 <- inherits(seurat_obj[[baseline_assay]], "Assay5")

  has_layer <- function(assay, layer) {
    if (is_v5) {
      layer %in% SeuratObject::Layers(seurat_obj[[assay]])
    } else {
      layer %in% c("counts", "data", "scale.data")
    }
  }

  as_sparse_or_stop <- function(x, label) {
    if (inherits(x, "sparseMatrix")) {
      return(methods::as(x, "CsparseMatrix"))
    }
    if (!isTRUE(allow_dense)) {
      stop(
        "Layer '", label, "' is dense. Set allow_dense = TRUE only if the ",
        "matrix is small enough to convert safely."
      )
    }
    Matrix::Matrix(x, sparse = TRUE)
  }

  get_layer <- function(assay, layer) {
    if (!has_layer(assay, layer)) {
      stop("Assay '", assay, "' is missing required layer '", layer, "'.")
    }
    value <- if (is_v5) {
      Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
    } else {
      Seurat::GetAssayData(seurat_obj, assay = assay, slot = layer)
    }
    if (nrow(value) == 0L || ncol(value) == 0L) {
      stop("Assay '", assay, "' has an empty '", layer, "' layer.")
    }
    as_sparse_or_stop(value, paste0(assay, ":", layer))
  }

  align_to_baseline <- function(x, assay, layer, features, cells) {
    missing_features <- setdiff(features, rownames(x))
    extra_features <- setdiff(rownames(x), features)
    if (length(missing_features) > 0L) {
      stop(
        "Assay '", assay, "' is missing ", length(missing_features),
        " baseline feature(s) in layer '", layer, "'."
      )
    }
    if (length(extra_features) > 0L) {
      warning(
        "Assay '", assay, "' has ", length(extra_features),
        " extra feature(s) in layer '", layer,
        "'; they will be omitted because they are absent from the baseline.",
        call. = FALSE
      )
    }
    if (!setequal(colnames(x), cells)) {
      stop("Assay '", assay, "' has different cells in layer '", layer, "'.")
    }
    x[features, cells, drop = FALSE]
  }

  log_normalize_fixed <- function(x, library_sizes) {
    if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
      stop(
        "Baseline counts contain cells with non-positive or non-finite library ",
        "sizes; fixed-denominator log-normalization is undefined."
      )
    }
    out <- methods::as(x, "CsparseMatrix")
    multipliers <- scale.factor / library_sizes
    out@x <- out@x * rep(multipliers, diff(out@p))
    out@x <- log1p(out@x)
    Matrix::drop0(out)
  }

  say("Loading baseline counts and reconstructing the workflow's baseline data.")
  base_counts <- get_layer(baseline_assay, "counts")
  features <- rownames(base_counts)
  cells <- colnames(base_counts)
  library_sizes <- Matrix::colSums(base_counts)
  base_log <- log_normalize_fixed(base_counts, library_sizes)

  # Each component records its direct count-space genes and its remaining
  # propagation-only data delta. Keeping this per assay preserves secondary
  # effects when a gene is direct in one assay and secondary in another.
  components <- vector("list", length(perturb_assays))
  names(components) <- perturb_assays

  for (assay in perturb_assays) {
    say("Reading component assay '", assay, "'.")
    pert_counts <- align_to_baseline(
      get_layer(assay, "counts"), assay, "counts", features, cells
    )
    pert_data <- align_to_baseline(
      get_layer(assay, "data"), assay, "data", features, cells
    )

    count_delta <- methods::as(pert_counts - base_counts, "CsparseMatrix")
    direct_mask <- Matrix::rowSums(abs(count_delta) > direct_tolerance) > 0
    direct_genes <- features[direct_mask]

    propagated_delta <- methods::as(pert_data - base_log, "CsparseMatrix")
    if (any(direct_mask)) {
      # ApplyPropagation pins direct rows to their primary log delta. That
      # primary component will instead be rebuilt from the combined counts.
      propagated_delta[direct_mask, ] <- 0
      propagated_delta <- Matrix::drop0(propagated_delta)
    }

    components[[assay]] <- list(
      count_delta = count_delta,
      propagated_delta = propagated_delta,
      direct_genes = direct_genes
    )
    say(
      "Inferred ", length(direct_genes), " direct gene(s) in '", assay,
      "'; retaining data-space deltas for secondary genes."
    )
  }

  say("Combining direct perturbations in count space.")
  combined_counts <- base_counts
  for (assay in perturb_assays) {
    combined_counts <- combined_counts + components[[assay]]$count_delta
  }
  combined_counts <- methods::as(combined_counts, "CsparseMatrix")

  negative_counts <- sum(combined_counts@x < 0)
  if (negative_counts > 0L) {
    if (!isTRUE(floor_counts)) {
      stop(
        "Combined count deltas produced ", negative_counts,
        " negative stored value(s). Use floor_counts = TRUE or revise the ",
        "overlapping perturbations."
      )
    }
    warning(
      "Combined count deltas produced ", negative_counts,
      " negative stored value(s); flooring them to zero.",
      call. = FALSE
    )
    combined_counts@x[combined_counts@x < 0] <- 0
    combined_counts <- Matrix::drop0(combined_counts)
  }

  say("Rebuilding direct data effects and adding propagated secondary effects.")
  combined_data <- log_normalize_fixed(combined_counts, library_sizes)
  for (assay in perturb_assays) {
    combined_data <- combined_data + components[[assay]]$propagated_delta
  }
  combined_data <- methods::as(combined_data, "CsparseMatrix")

  negative_data <- sum(combined_data@x < 0)
  if (negative_data > 0L) {
    if (!isTRUE(floor_data)) {
      stop(
        "Combined propagated effects produced ", negative_data,
        " negative stored data value(s). Use floor_data = TRUE or revise the ",
        "overlapping perturbations."
      )
    }
    warning(
      "Combined propagated effects produced ", negative_data,
      " negative stored data value(s); flooring them to zero.",
      call. = FALSE
    )
    combined_data@x[combined_data@x < 0] <- 0
    combined_data <- Matrix::drop0(combined_data)
  }

  say("Writing combined assay '", new_assay, "'.")
  if (is_v5) {
    seurat_obj[[new_assay]] <- SeuratObject::CreateAssay5Object(
      counts = combined_counts
    )
    seurat_obj <- SeuratObject::SetAssayData(
      seurat_obj,
      assay = new_assay,
      layer = "data",
      new.data = combined_data
    )
  } else {
    seurat_obj[[new_assay]] <- Seurat::CreateAssayObject(
      counts = combined_counts
    )
    seurat_obj <- Seurat::SetAssayData(
      seurat_obj,
      assay = new_assay,
      slot = "data",
      new.data = combined_data
    )
  }

  if (is.null(seurat_obj[[new_assay]]@misc)) {
    seurat_obj[[new_assay]]@misc <- list()
  }
  seurat_obj[[new_assay]]@misc$COMPACT <- list(
    baseline_assay = baseline_assay,
    baseline_layers = c("counts", "data"),
    method = "CombinePerturbAssays (direct counts + secondary data deltas)",
    components = perturb_assays,
    direct_genes_by_component = lapply(components, `[[`, "direct_genes"),
    scale.factor = scale.factor,
    direct_tolerance = direct_tolerance,
    provenance_check = provenance_summary,
    floor_counts = floor_counts,
    floor_data = floor_data,
    created_at = as.character(Sys.time())
  )

  say("Done. The original CombinePerturbAssays() implementation was not changed.")
  seurat_obj
}
