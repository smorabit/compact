#' Project an assay into an existing PCA reference space
#'
#' @description
#' Projects normalized expression from a query assay into a PCA model fitted on
#' a reference assay. The original PCA feature set, reference gene means and
#' standard deviations, and PCA loadings are reused so that reference and query
#' coordinates share the same geometric ruler.
#'
#' @details
#' `ProjectPCAReference()` is intended for comparisons such as observed versus
#' in-silico perturbed expression stored as separate assays in the same Seurat
#' object. It does not refit PCA and does not modify either assay or the original
#' reduction.
#'
#' By default, reference scaling parameters are reconstructed from the
#' `reference_assay` `data` layer. The reconstructed scaled values and PCA scores
#' must reproduce the stored `scale.data` and reference PCA embeddings within
#' `tolerance`. This validation deliberately rejects PCA models produced with
#' incompatible preprocessing, such as regression in `ScaleData()`, rather than
#' silently creating a non-comparable projection.
#'
#' The query layer must already use the same normalization definition as the
#' reference layer. In particular, the function does not call `NormalizeData()`.
#'
#' @param seurat_object A `Seurat` object containing the reference assay, query
#'   assay, and reference PCA reduction.
#' @param query_assay Name of the assay to project.
#' @param reference_assay Name of the assay used to fit the reference PCA.
#' @param reference_reduction Name of the stored PCA reduction to reuse.
#' @param dims Leading sequence of reference PCs to project, such as `1:30`. If
#'   `NULL`, all stored loading columns are used. Seurat reductions require
#'   consecutive dimension identifiers, so skipped or reordered PCs are rejected.
#' @param layer Expression layer used to reconstruct scaling and project the
#'   query. Usually `"data"` for log-normalized expression.
#' @param scale.max Upper clipping value used when the reference assay was
#'   scaled. This must match the original `ScaleData()` call. Use `Inf` for no
#'   clipping. Seurat clips values greater than `scale.max`.
#' @param reduction.name Name for the projected reduction. Defaults to
#'   `paste0(query_assay, "_pca_fixed")`.
#' @param reduction.key Optional key for the projected dimensions. If `NULL`, a
#'   valid key is derived from `reduction.name`.
#' @param unchanged.cells Optional character vector of cells known to be
#'   unchanged between the reference and query assays. Their expression and
#'   projected coordinates are checked within `tolerance`.
#' @param tolerance Maximum allowed absolute reconstruction or unchanged-cell
#'   error.
#' @param block.size Number of query cells projected per block.
#' @param overwrite Logical; overwrite an existing `reduction.name` when `TRUE`.
#' @param verbose Logical; print progress and validation messages.
#'
#' @return The input Seurat object with the fixed-reference query coordinates
#'   stored as a new dimensional reduction. Projection parameters and validation
#'   results are stored in the reduction's `misc$fixed_reference` field.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' object <- ProjectPCAReference(
#'   object,
#'   reference_assay = "RNA",
#'   query_assay = "M6_up",
#'   reference_reduction = "pca",
#'   dims = 1:30,
#'   unchanged.cells = WhichCells(object, expression = treatment == "Bortezomib")
#' )
#'
#' pre <- ComputeDistance(
#'   object, "treatment_subcluster", "pca", method = "edist", dims = 1:30
#' )
#' post <- ComputeDistance(
#'   object, "treatment_subcluster", "M6_up_pca_fixed",
#'   method = "edist", dims = 1:30
#' )
#' HeatmapDistance(pre, post)
#' }
ProjectPCAReference <- function(
    seurat_object,
    query_assay,
    reference_assay = "RNA",
    reference_reduction = "pca",
    dims = NULL,
    layer = "data",
    scale.max = 10,
    reduction.name = NULL,
    reduction.key = NULL,
    unchanged.cells = NULL,
    tolerance = 1e-6,
    block.size = 5000L,
    overwrite = FALSE,
    verbose = TRUE) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("`seurat_object` must be a Seurat object.", call. = FALSE)
  }

  .pca_reference_validate_string(query_assay, "query_assay")
  .pca_reference_validate_string(reference_assay, "reference_assay")
  .pca_reference_validate_string(reference_reduction, "reference_reduction")
  .pca_reference_validate_string(layer, "layer")

  assays <- names(seurat_object@assays)
  if (!(reference_assay %in% assays)) {
    stop("Reference assay '", reference_assay, "' was not found.", call. = FALSE)
  }
  if (!(query_assay %in% assays)) {
    stop("Query assay '", query_assay, "' was not found.", call. = FALSE)
  }
  if (!(reference_reduction %in% names(seurat_object@reductions))) {
    stop(
      "Reference reduction '", reference_reduction, "' was not found.",
      call. = FALSE
    )
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      is.na(tolerance) || tolerance < 0) {
    stop("`tolerance` must be one non-negative number.", call. = FALSE)
  }
  if (!is.numeric(scale.max) || length(scale.max) != 1L ||
      is.na(scale.max) || scale.max <= 0) {
    stop("`scale.max` must be one positive number or Inf.", call. = FALSE)
  }
  if (!is.numeric(block.size) || length(block.size) != 1L ||
      is.na(block.size) || block.size < 1 || block.size %% 1 != 0) {
    stop("`block.size` must be one positive integer.", call. = FALSE)
  }
  block.size <- as.integer(block.size)

  reference_dr <- seurat_object[[reference_reduction]]
  reference_loadings <- Seurat::Loadings(reference_dr, projected = FALSE)
  if (is.null(reference_loadings) || nrow(reference_loadings) == 0L ||
      ncol(reference_loadings) == 0L) {
    stop(
      "Reference reduction '", reference_reduction,
      "' does not contain feature loadings.",
      call. = FALSE
    )
  }

  if (is.null(dims)) {
    dims <- seq_len(ncol(reference_loadings))
  }
  if (!is.numeric(dims) || length(dims) == 0L || anyNA(dims) ||
      any(dims %% 1 != 0) || any(dims < 1) ||
      any(dims > ncol(reference_loadings)) || anyDuplicated(dims)) {
    stop(
      "`dims` must contain unique integer PC indices between 1 and ",
      ncol(reference_loadings), ".",
      call. = FALSE
    )
  }
  dims <- as.integer(dims)
  if (!identical(dims, seq_len(length(dims)))) {
    stop(
      "`dims` must be a leading consecutive PC sequence such as 1:30. ",
      "Skipped or reordered dimensions cannot be stored safely as a Seurat ",
      "reduction.",
      call. = FALSE
    )
  }

  features <- rownames(reference_loadings)
  if (is.null(features) || anyNA(features) || any(!nzchar(features)) ||
      anyDuplicated(features)) {
    stop("Reference PCA loadings must have unique feature names.", call. = FALSE)
  }

  reference_data <- .pca_reference_get_data(
    seurat_object, reference_assay, layer
  )
  query_data <- .pca_reference_get_data(seurat_object, query_assay, layer)
  reference_scaled_stored <- .pca_reference_get_data(
    seurat_object, reference_assay, "scale.data"
  )

  missing_reference <- setdiff(features, rownames(reference_data))
  missing_query <- setdiff(features, rownames(query_data))
  missing_scaled <- setdiff(features, rownames(reference_scaled_stored))
  if (length(missing_reference) > 0L) {
    stop(
      length(missing_reference), " PCA feature(s) are absent from reference assay '",
      reference_assay, "'. First missing feature: ", missing_reference[[1]],
      call. = FALSE
    )
  }
  if (length(missing_query) > 0L) {
    stop(
      length(missing_query), " PCA feature(s) are absent from query assay '",
      query_assay, "'. First missing feature: ", missing_query[[1]],
      call. = FALSE
    )
  }
  if (length(missing_scaled) > 0L) {
    stop(
      length(missing_scaled), " PCA feature(s) are absent from the reference ",
      "scale.data layer, so the original PCA scaling model cannot be safely ",
      "validated. This commonly happens when a later ScaleData() call (for ",
      "example, during hdWGCNA) overwrites scale.data.\n\n",
      "Recommended action: keep the current object unchanged, create a copy, ",
      "recover the original PCA genes with ",
      "rownames(Seurat::Loadings(object[[reference_reduction]])), rerun ",
      "ScaleData() on the reference assay using exactly those genes and the ",
      "original scaling settings, and save a new PCA reduction such as ",
      "'pca_fixed_reference'. Then rerun ProjectPCAReference() with ",
      "reference_reduction = 'pca_fixed_reference'. Do not bypass this check ",
      "by increasing tolerance.",
      call. = FALSE
    )
  }

  reference_data <- reference_data[features, , drop = FALSE]
  query_data <- query_data[features, , drop = FALSE]
  reference_scaled_stored <- reference_scaled_stored[features, , drop = FALSE]

  reference_cells <- colnames(reference_data)
  query_cells <- colnames(query_data)
  if (is.null(reference_cells) || is.null(query_cells) ||
      !setequal(reference_cells, query_cells)) {
    stop(
      "Reference and query assays must contain the same named cells.",
      call. = FALSE
    )
  }
  query_data <- query_data[, reference_cells, drop = FALSE]
  reference_scaled_stored <- reference_scaled_stored[, reference_cells, drop = FALSE]

  stored_embeddings <- Seurat::Embeddings(
    seurat_object, reduction = reference_reduction
  )
  missing_embedding_cells <- setdiff(reference_cells, rownames(stored_embeddings))
  if (length(missing_embedding_cells) > 0L) {
    stop(
      "The reference PCA is missing ", length(missing_embedding_cells),
      " assay cell(s).",
      call. = FALSE
    )
  }
  stored_embeddings <- stored_embeddings[
    reference_cells, dims, drop = FALSE
  ]

  reference_assay_used <- tryCatch(
    methods::slot(reference_dr, "assay.used"),
    error = function(e) NULL
  )
  if (!is.null(reference_assay_used) && length(reference_assay_used) == 1L &&
      !is.na(reference_assay_used) && nzchar(reference_assay_used) &&
      !identical(reference_assay_used, reference_assay)) {
    stop(
      "Reference reduction '", reference_reduction, "' was fitted using assay '",
      reference_assay_used, "', not '", reference_assay, "'.",
      call. = FALSE
    )
  }

  if (verbose) {
    message(
      "Reconstructing reference scaling for ", length(features),
      " PCA features..."
    )
  }
  reference_stats <- .pca_reference_row_stats(reference_data)
  zero_scale <- !is.finite(reference_stats$scale) | reference_stats$scale <= 0
  if (any(zero_scale)) {
    stop(
      sum(zero_scale), " reference PCA feature(s) have zero or non-finite ",
      "standard deviation in the selected layer.",
      call. = FALSE
    )
  }

  selected_loadings <- as.matrix(
    reference_loadings[features, dims, drop = FALSE]
  )
  reference_projection <- .pca_reference_project_matrix(
    expression = reference_data,
    center = reference_stats$center,
    scale = reference_stats$scale,
    loadings = selected_loadings,
    scale.max = scale.max,
    block.size = block.size
  )

  scaling_error <- .pca_reference_scaling_error(
    expression = reference_data,
    stored_scaled = reference_scaled_stored,
    center = reference_stats$center,
    scale = reference_stats$scale,
    scale.max = scale.max,
    block.size = block.size
  )
  reconstruction_error <- max(
    abs(reference_projection - stored_embeddings), na.rm = TRUE
  )

  if (!is.finite(scaling_error) || scaling_error > tolerance) {
    stop(
      "Could not reproduce the reference scale.data layer (maximum absolute ",
      "error = ", signif(scaling_error, 6), ", tolerance = ", tolerance,
      "). The reference PCA may use regression, nondefault centering/scaling, ",
      "a different layer, or a different scale.max value.",
      call. = FALSE
    )
  }
  if (!is.finite(reconstruction_error) || reconstruction_error > tolerance) {
    stop(
      "Could not reproduce the stored reference PCA embeddings (maximum ",
      "absolute error = ", signif(reconstruction_error, 6),
      ", tolerance = ", tolerance, ").",
      call. = FALSE
    )
  }

  if (verbose) message("Projecting query assay into the fixed PCA reference...")
  query_projection <- .pca_reference_project_matrix(
    expression = query_data,
    center = reference_stats$center,
    scale = reference_stats$scale,
    loadings = selected_loadings,
    scale.max = scale.max,
    block.size = block.size
  )

  unchanged_expression_error <- NA_real_
  unchanged_coordinate_error <- NA_real_
  if (!is.null(unchanged.cells)) {
    if (!is.character(unchanged.cells) || length(unchanged.cells) == 0L ||
        anyNA(unchanged.cells) ||
        anyDuplicated(unchanged.cells)) {
      stop(
        "`unchanged.cells` must be a non-empty, unique character vector.",
        call. = FALSE
      )
    }
    missing_unchanged <- setdiff(unchanged.cells, reference_cells)
    if (length(missing_unchanged) > 0L) {
      stop(
        length(missing_unchanged), " unchanged cell(s) were not found. First ",
        "missing cell: ", missing_unchanged[[1]],
        call. = FALSE
      )
    }
    unchanged_expression_error <- max(
      abs(as.matrix(
        query_data[, unchanged.cells, drop = FALSE] -
          reference_data[, unchanged.cells, drop = FALSE]
      )),
      na.rm = TRUE
    )
    unchanged_coordinate_error <- max(
      abs(
        query_projection[unchanged.cells, , drop = FALSE] -
          stored_embeddings[unchanged.cells, , drop = FALSE]
      ),
      na.rm = TRUE
    )
    if (unchanged_expression_error > tolerance) {
      stop(
        "Cells supplied in `unchanged.cells` are not unchanged in the selected ",
        "expression layer (maximum absolute error = ",
        signif(unchanged_expression_error, 6), ").",
        call. = FALSE
      )
    }
    if (unchanged_coordinate_error > tolerance) {
      stop(
        "Projected coordinates for `unchanged.cells` differ from the reference ",
        "coordinates (maximum absolute error = ",
        signif(unchanged_coordinate_error, 6), ").",
        call. = FALSE
      )
    }
  }

  if (is.null(reduction.name)) {
    reduction.name <- paste0(query_assay, "_pca_fixed")
  }
  .pca_reference_validate_string(reduction.name, "reduction.name")
  if (reduction.name %in% names(seurat_object@reductions) && !isTRUE(overwrite)) {
    stop(
      "Reduction '", reduction.name,
      "' already exists. Choose another name or set `overwrite = TRUE`.",
      call. = FALSE
    )
  }

  if (is.null(reduction.key)) {
    key_base <- gsub("[^A-Za-z0-9]", "", reduction.name)
    if (!nzchar(key_base)) key_base <- "PCFIXED"
    if (!grepl("^[A-Za-z]", key_base)) key_base <- paste0("X", key_base)
    reduction.key <- paste0(toupper(key_base), "_")
  }
  .pca_reference_validate_string(reduction.key, "reduction.key")
  if (!grepl("^[A-Za-z][A-Za-z0-9]*_$", reduction.key)) {
    stop(
      "`reduction.key` must begin with a letter, contain only letters or ",
      "numbers, and end with an underscore.",
      call. = FALSE
    )
  }

  dimension_names <- paste0(reduction.key, dims)
  colnames(query_projection) <- dimension_names
  colnames(selected_loadings) <- dimension_names

  reference_stdev <- tryCatch(
    methods::slot(reference_dr, "stdev")[dims],
    error = function(e) rep(NA_real_, length(dims))
  )
  if (length(reference_stdev) != length(dims)) {
    reference_stdev <- rep(NA_real_, length(dims))
  }

  provenance <- list(
    method = "fixed-reference PCA projection",
    reference_assay = reference_assay,
    query_assay = query_assay,
    reference_reduction = reference_reduction,
    layer = layer,
    features = features,
    dims = dims,
    center = reference_stats$center,
    scale = reference_stats$scale,
    scale.max = scale.max,
    validation = list(
      scaling_error = scaling_error,
      reference_reconstruction_error = reconstruction_error,
      unchanged_expression_error = unchanged_expression_error,
      unchanged_coordinate_error = unchanged_coordinate_error,
      tolerance = tolerance,
      unchanged_cells = unchanged.cells
    )
  )

  seurat_object[[reduction.name]] <- Seurat::CreateDimReducObject(
    embeddings = query_projection,
    loadings = selected_loadings,
    stdev = reference_stdev,
    key = reduction.key,
    assay = query_assay,
    misc = list(fixed_reference = provenance)
  )

  if (verbose) {
    message(
      "Stored fixed-reference projection as reduction '", reduction.name,
      "' (maximum reference reconstruction error: ",
      signif(reconstruction_error, 4), ")."
    )
  }

  seurat_object
}


.pca_reference_validate_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be one non-empty character string.", call. = FALSE)
  }
  invisible(TRUE)
}


.pca_reference_get_data <- function(seurat_object, assay, layer) {
  tryCatch(
    suppressWarnings(
      Seurat::GetAssayData(seurat_object, assay = assay, slot = layer)
    ),
    error = function(slot_error) {
      tryCatch(
        Seurat::GetAssayData(seurat_object, assay = assay, layer = layer),
        error = function(layer_error) {
          stop(
            "Could not read layer/slot '", layer, "' from assay '", assay,
            "'. Slot error: ", conditionMessage(slot_error),
            "; layer error: ", conditionMessage(layer_error),
            call. = FALSE
          )
        }
      )
    }
  )
}


.pca_reference_row_stats <- function(expression) {
  n <- ncol(expression)
  if (n < 2L) {
    stop("At least two reference cells are required.", call. = FALSE)
  }

  center <- Matrix::rowMeans(expression)
  sum_squares <- Matrix::rowSums(expression ^ 2)
  variance <- (sum_squares - n * center ^ 2) / (n - 1)
  variance[variance < 0 & variance > -sqrt(.Machine$double.eps)] <- 0
  scale <- sqrt(variance)
  names(center) <- names(scale) <- rownames(expression)
  list(center = center, scale = scale)
}


.pca_reference_scale_block <- function(expression, center, scale, scale.max) {
  scaled <- as.matrix(expression)
  scaled <- sweep(scaled, 1L, center, FUN = "-")
  scaled <- sweep(scaled, 1L, scale, FUN = "/")
  scaled[!is.finite(scaled)] <- 0
  if (is.finite(scale.max)) scaled[scaled > scale.max] <- scale.max
  scaled
}


.pca_reference_project_matrix <- function(
    expression, center, scale, loadings, scale.max, block.size) {
  cells <- colnames(expression)
  output <- matrix(
    NA_real_,
    nrow = length(cells),
    ncol = ncol(loadings),
    dimnames = list(cells, colnames(loadings))
  )

  starts <- seq.int(1L, length(cells), by = block.size)
  for (start in starts) {
    index <- start:min(start + block.size - 1L, length(cells))
    scaled <- .pca_reference_scale_block(
      expression[, index, drop = FALSE], center, scale, scale.max
    )
    output[index, ] <- t(scaled) %*% loadings
  }
  output
}


.pca_reference_scaling_error <- function(
    expression, stored_scaled, center, scale, scale.max, block.size) {
  n_cells <- ncol(expression)
  maximum_error <- 0
  starts <- seq.int(1L, n_cells, by = block.size)
  for (start in starts) {
    index <- start:min(start + block.size - 1L, n_cells)
    reconstructed <- .pca_reference_scale_block(
      expression[, index, drop = FALSE], center, scale, scale.max
    )
    current_error <- max(
      abs(reconstructed - as.matrix(stored_scaled[, index, drop = FALSE])),
      na.rm = TRUE
    )
    maximum_error <- max(maximum_error, current_error)
  }
  maximum_error
}
