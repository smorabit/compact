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

  seurat_obj[[perturb_assay]]@misc$COMPACT$baseline_assay   <- baseline_assay
  # seurat_obj[[perturb_assay]]@misc$COMPACT$baseline_layers <- c("counts", "data")
  seurat_obj[[perturb_assay]]@misc$COMPACT$baseline_layers  <- c("counts")
  seurat_obj[[perturb_assay]]@misc$COMPACT$method           <- method
  seurat_obj[[perturb_assay]]@misc$COMPACT$stamped_at       <- as.character(Sys.time())

  seurat_obj
}
#
