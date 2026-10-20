# ============================================================
# COMPACT Alert System helpers — standalone preview
# ============================================================
#
# This file intentionally contains only public I/O and pathway-loading helpers.
# The obsolete internal scoring helpers from AlertSystemScore_OLD.R are not
# restored because their work is performed directly by the updated scorer.
#
# Suggested package location after review: R/AlertSystemHelper.R


# -----------------------------------------------------------------------------
# Internal validation shared by SaveAlertSystem() and LoadAlertSystem().
# -----------------------------------------------------------------------------
# An Alert System state is expected to contain the long-format per-cell table,
# pathway-level metrics, run parameters, and (for the updated schema) the
# filtered gene sets used for scoring.
.validate_alert_system_state <- function(state, require_gene_sets = FALSE) {
  if (!is.list(state)) {
    stop("Alert System state must be a list.", call. = FALSE)
  }

  required <- c("per_cell", "metrics", "params")
  if (isTRUE(require_gene_sets)) {
    required <- c(required, "gene_sets")
  }
  missing <- setdiff(required, names(state))
  if (length(missing) > 0L) {
    stop(
      "Invalid Alert System state; missing element(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.data.frame(state$per_cell)) {
    stop("Alert System state 'per_cell' must be a data frame.", call. = FALSE)
  }
  if (!is.data.frame(state$metrics)) {
    stop("Alert System state 'metrics' must be a data frame.", call. = FALSE)
  }
  if (!is.list(state$params)) {
    stop("Alert System state 'params' must be a list.", call. = FALSE)
  }
  if (!is.null(state$gene_sets) &&
      (!is.list(state$gene_sets) || is.null(names(state$gene_sets)))) {
    stop("Alert System state 'gene_sets' must be a named list.", call. = FALSE)
  }

  invisible(state)
}


#' Save an Alert System analysis state
#'
#' Saves the complete state stored in `seurat_obj@misc[[misc_name]]`. The saved
#' RDS can later be attached to another compatible Seurat object using
#' [LoadAlertSystem()] without recomputing pathway scores.
#'
#' @param seurat_obj A Seurat object containing an Alert System state.
#' @param file Destination `.rds` path.
#' @param misc_name Name of the `@misc` entry. Default is `"AlertSystem"`.
#' @param overwrite Whether an existing file may be replaced.
#'
#' @return Invisibly returns the normalized output path.
#' @export
SaveAlertSystem <- function(
  seurat_obj,
  file,
  misc_name = "AlertSystem",
  overwrite = FALSE
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }
  if (!is.character(misc_name) || length(misc_name) != 1L ||
      is.na(misc_name) || !nzchar(misc_name)) {
    stop("misc_name must be one non-empty character value.", call. = FALSE)
  }
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("file must be one non-empty path.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.", call. = FALSE)
  }
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop("File already exists; set overwrite = TRUE to replace it: ", file,
         call. = FALSE)
  }

  state <- seurat_obj@misc[[misc_name]]
  if (is.null(state)) {
    stop(
      "No Alert System state found in seurat_obj@misc[['", misc_name, "']].",
      call. = FALSE
    )
  }
  # Allow legacy states without gene_sets, but validate the common core.
  .validate_alert_system_state(state, require_gene_sets = FALSE)

  parent <- dirname(file)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(parent)) {
    stop("Could not create output directory: ", parent, call. = FALSE)
  }

  saveRDS(state, file = file)
  invisible(normalizePath(file, mustWork = TRUE))
}


#' Load an Alert System analysis state
#'
#' Reads an RDS created by [SaveAlertSystem()] and attaches it to a Seurat
#' object. This restores analysis results only; it does not add per-cell score
#' columns to `seurat_obj@meta.data`.
#'
#' @param seurat_obj A Seurat object.
#' @param file Path to an Alert System `.rds` state.
#' @param misc_name Destination entry in `seurat_obj@misc`.
#' @param overwrite Whether an existing state at `misc_name` may be replaced.
#' @param require_gene_sets Require the updated schema's `gene_sets` element.
#'   Keep `FALSE` to permit legacy Alert System states.
#' @param verbose Whether to report where the state was attached.
#'
#' @return The updated Seurat object.
#' @export
LoadAlertSystem <- function(
  seurat_obj,
  file,
  misc_name = "AlertSystem",
  overwrite = TRUE,
  require_gene_sets = FALSE,
  verbose = TRUE
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("Alert System RDS file not found: ", file, call. = FALSE)
  }
  if (!is.character(misc_name) || length(misc_name) != 1L ||
      is.na(misc_name) || !nzchar(misc_name)) {
    stop("misc_name must be one non-empty character value.", call. = FALSE)
  }
  if (!isTRUE(overwrite) && !is.null(seurat_obj@misc[[misc_name]])) {
    stop(
      "seurat_obj@misc[['", misc_name,
      "']] already exists and overwrite = FALSE.",
      call. = FALSE
    )
  }

  state <- readRDS(file)
  .validate_alert_system_state(
    state,
    require_gene_sets = isTRUE(require_gene_sets)
  )
  seurat_obj@misc[[misc_name]] <- state

  if (isTRUE(verbose)) {
    message(
      "Loaded Alert System state into seurat_obj@misc[['",
      misc_name, "']]."
    )
  }
  seurat_obj
}


# -----------------------------------------------------------------------------
# Validate and normalize pathway collections from any source.
# -----------------------------------------------------------------------------
.validate_alert_gene_sets <- function(gene_sets, source_label = "gene-set file") {
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || length(gene_sets) == 0L) {
    stop(source_label, " must contain a non-empty named list.", call. = FALSE)
  }
  if (anyNA(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
    stop(source_label, " contains missing or empty pathway names.", call. = FALSE)
  }
  if (anyDuplicated(names(gene_sets))) {
    stop(source_label, " contains duplicated pathway names.", call. = FALSE)
  }

  normalized <- lapply(gene_sets, function(genes) {
    genes <- unique(as.character(genes))
    genes[!is.na(genes) & nzchar(genes)]
  })
  empty <- names(normalized)[lengths(normalized) == 0L]
  if (length(empty) > 0L) {
    stop(
      source_label, " contains empty pathway(s): ",
      paste(empty, collapse = ", "), ".",
      call. = FALSE
    )
  }
  normalized
}


#' Load bundled or user-supplied Alert System pathways
#'
#' Lookup order is: an explicitly supplied file/directory, the installed
#' package's `inst/extdata`, and—only when allowed—a GitHub fallback.
#'
#' @param path Optional `.rds` file or directory containing the expected file.
#' @param genetype Either `"mouse"` or `"human"`.
#' @param database Filename prefix; default `"MSigDBhallmark"`.
#' @param package Package containing bundled pathway files.
#' @param allow_download Permit the GitHub fallback when no local/package file
#'   exists.
#' @param github_base Base raw-content URL used for the optional fallback.
#' @param verbose Whether to report the selected source.
#'
#' @return A validated named list of pathway gene vectors.
#' @export
LoadAlertSystemPathways <- function(
  path = NULL,
  genetype = c("mouse", "human"),
  database = "MSigDBhallmark",
  package = "compact",
  allow_download = TRUE,
  github_base = paste0(
    "https://raw.githubusercontent.com/smorabit/compact/",
    "main/inst/extdata"
  ),
  verbose = TRUE
) {
  genetype <- match.arg(genetype)
  if (!is.character(database) || length(database) != 1L ||
      is.na(database) || !nzchar(database)) {
    stop("database must be one non-empty character value.", call. = FALSE)
  }
  filename <- paste0(database, "_gene_sets_", genetype, ".rds")

  # 1. An explicit path always takes precedence.
  if (!is.null(path)) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
      stop("path must be NULL or one non-empty path.", call. = FALSE)
    }
    candidate <- if (grepl("\\.rds$", path, ignore.case = TRUE)) {
      path
    } else {
      file.path(path, filename)
    }
    if (!file.exists(candidate)) {
      stop("Gene-set RDS file not found: ", candidate, call. = FALSE)
    }
    if (isTRUE(verbose)) message("Loading Alert System pathways from: ", candidate)
    return(.validate_alert_gene_sets(readRDS(candidate), candidate))
  }

  # 2. Prefer package-bundled data; this is reproducible and works offline.
  candidate <- system.file(
    "extdata", filename, package = package, mustWork = FALSE
  )
  if (nzchar(candidate) && file.exists(candidate)) {
    if (isTRUE(verbose)) message("Loading Alert System pathways from package: ", candidate)
    return(.validate_alert_gene_sets(readRDS(candidate), candidate))
  }

  # 3. Network access is an explicit fallback, not the primary data source.
  if (!isTRUE(allow_download)) {
    stop(
      "Pathway file was not supplied or bundled, and allow_download = FALSE: ",
      filename,
      call. = FALSE
    )
  }
  url <- paste0(sub("/$", "", github_base), "/", filename)
  if (isTRUE(verbose)) message("Downloading Alert System pathways from: ", url)
  temporary_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temporary_file), add = TRUE)

  tryCatch(
    utils::download.file(url, temporary_file, mode = "wb", quiet = !isTRUE(verbose)),
    error = function(e) {
      stop("Failed to download pathway file: ", conditionMessage(e), call. = FALSE)
    }
  )
  .validate_alert_gene_sets(readRDS(temporary_file), url)
}
