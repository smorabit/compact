# Description: Computes E-test statistics and Compute Distance for each group in a Seurat object

library(Seurat)
library(dplyr)
library(rdist)
library(energy)

#' @title edist
#'
#' @description
#' Computes pairwise Energy distances between groups (cell states) in a Seurat object.
#'
#' @details
#' Cells are grouped according to \code{groupby}. For each pair of groups, the empirical
#' Energy distance is computed from pairwise distances between individual cell embeddings
#' in the specified \code{reduction} space. The resulting output is a symmetric
#' groups \eqn{\times} groups distance matrix.
#'
#' @param seurat_object A \code{Seurat} object.
#' @param groupby A character scalar giving the name of a metadata column in
#'   \code{seurat_object@meta.data} that defines group labels.
#' @param reduction A character scalar giving the name of a dimensional reduction
#'   (e.g., \code{"pca"}, \code{"harmony"}) used for distance calculations. This should be
#'   a geometry-preserving representation; visualization-oriented nonlinear embeddings
#'   (e.g., \code{"umap"}, \code{"tsne"}) are not recommended for distance computation.
#' @param sample_correct Logical; if \code{TRUE}, uses the sample-size corrected
#'   (U-statistic) normalization for within-group terms (dividing by \eqn{N(N-1)} rather
#'   than \eqn{N^2}). Note that Energy distance can be negative due to finite-sample
#'   estimation noise even though the population quantity is non-negative.
#' @param squared Logical; if \code{TRUE}, uses squared Euclidean distances in the
#'   embedding space (matching \code{metric="sqeuclidean"} in some Python implementations).
#'   If \code{FALSE} (default), uses Euclidean distances.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' @param dims Optional integer vector selecting dimensions from \code{reduction}.
#'   If \code{NULL}, all dimensions are used.
#'
#' @return A symmetric \code{data.frame} (groups \eqn{\times} groups) of pairwise
#'   Energy distances.
#'
#' @importFrom Seurat Embeddings
#' @importFrom rdist cdist
.edist <- function(seurat_object, groupby, reduction,
                   sample_correct = TRUE,
                   squared = FALSE,
                   verbose = TRUE,
                   dims = NULL) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("The first argument must be a Seurat object.")
  }
  if (!(reduction %in% names(seurat_object@reductions))) {
    stop(
      sprintf(
        "Reduction '%s' was not found in the Seurat object. ",
        reduction
      ),
      "Please provide a valid geometry-preserving dimensional reduction ",
      "(e.g., 'pca' or 'harmony') stored in seurat_object@reductions. ",
      "Visualization-oriented nonlinear embeddings (e.g., UMAP or t-SNE) ",
      "are not suitable for distance-based calculations."
    )
  }


  # ---- extract + clean labels / embeddings ----
  labels_raw <- seurat_object[[]][[groupby]]
  if (is.null(labels_raw)) {
    stop("groupby column not found in Seurat metadata.")
  }

  labels <- as.character(labels_raw)
  keep <- !is.na(labels)
  labels <- labels[keep]
  emb <- Seurat::Embeddings(seurat_object, reduction = reduction)[keep, , drop = FALSE]
  emb <- .distance_select_dims(emb, dims)

  groups <- sort(unique(labels))

  df <- matrix(NA_real_, length(groups), length(groups),
               dimnames = list(groups, groups))

  if (verbose) message("Computing Energy distances between groups...")

  for (i in seq_along(groups)) {
    gx <- groups[i]
    x <- as.matrix(emb[labels == gx, , drop = FALSE])
    N <- nrow(x)

    if (N < 2) next
    df[gx, gx] <- 0

    for (j in (i + 1):length(groups)) {
      if (j > length(groups)) break

      gy <- groups[j]
      y <- as.matrix(emb[labels == gy, , drop = FALSE])
      M <- nrow(y)

      if (M < 2) next

      dxy <- rdist::cdist(x, y)
      dxx <- rdist::cdist(x, x)
      dyy <- rdist::cdist(y, y)

      if (squared) {
        dxy <- dxy^2
        dxx <- dxx^2
        dyy <- dyy^2
      }

      delta_xy <- sum(dxy) / (N * M)

      if (sample_correct) {
        delta_x <- sum(dxx) / (N * (N - 1))   # diag is 0 so OK
        delta_y <- sum(dyy) / (M * (M - 1))
      } else {
        delta_x <- sum(dxx) / (N * N)
        delta_y <- sum(dyy) / (M * M)
      }

      ed <- 2 * delta_xy - delta_x - delta_y
      df[gx, gy] <- df[gy, gx] <- ed


      # NOTE cdist vs pdist
      # rdist::cdist(x,x) gives you ordered pairs (full matrix)
      # rdist::pdist(x) gives you unordered pairs (half the matrix)

      # dist_xy <- rdist::cdist(x,y)
      # dist_x <- rdist::pdist(x)   ## unordered pairs: k<l
      # dist_y <- rdist::pdist(y)   ## unordered pairs: k<l
      #
      # if (sample_correction) {
      #     ed <- 2 * (sum(dist_xy) / (N * M)) - (sum(dist_x) / (N * (N - 1))) - (sum(dist_y) / (M * (M - 1)))
      # } else {
      #     ed <- 2 * (sum(dist_xy) / (N * M)) - (sum(dist_x) / (N * N)) - (sum(dist_y) / (M * M))
      # }


    }
  }

  as.data.frame(df)

}




# .edist <- function(seurat_object, groupby, reduction,
#                   sample_correction = FALSE, verbose = TRUE) {
#     if (!inherits(seurat_object, "Seurat")) {
#         stop("The first argument must be a Seurat object.")
#     }
#     if (!(reduction %in% names(seurat_object@reductions))) {
#         if (verbose) {
#             message("The specified reduction was not found in the Seurat object. Ensure that the applied dimensional reduction method is linear in nature, such as PCA or Harmony on the selected linear model.")
#         }
#         stop("The specified reduction was not found in the Seurat object.")
#     }
#     labels <- seurat_object[[]][[groupby]]
#     groups <- unique(labels)
#     emb <- Seurat::Embeddings(seurat_object, reduction = reduction)
#
#     df <- setNames(data.frame(matrix(ncol = length(groups),
#                                      nrow = length(groups)),
#                               row.names = groups),
#                    groups)
#     if (verbose) {
#         print("Computing E-test statistics for each group.")
#     }
#     completed_groups <- c()
#     for (groupx in groups) {
#         for (groupy in groups){
#             if (groupy %in% completed_groups) {
#                 next
#             }  # skip if already computed
#             x <- as.matrix(emb)[labels == groupx, ]
#             y <- as.matrix(emb)[labels == groupy, ]
#
#             N <- nrow(x)
#             M <- nrow(y)
#
#             dist_xy <- rdist::cdist(x,y)
#             dist_x <- rdist::pdist(x) # bug in scPerturb
#             dist_y <- rdist::pdist(y)
#
#             if (sample_correction) {
#                 ed <- 2 * (sum(dist_xy) / (N * M)) - (sum(dist_x) / (N * (N - 1))) - (sum(dist_y) / (M * (M - 1)))
#             } else {
#                 ed <- 2 * (sum(dist_xy) / (N * M)) - (sum(dist_x) / (N * N)) - (sum(dist_y) / (M * M))
#             }
#
#             df[groupx, groupy] <- df[groupy, groupx] <- ed
#         }
#         completed_groups <- c(completed_groups, groupx)
#     }
#     return(df)
# }


#' @title eucldist
#'
#' @description
#' Computes pairwise Euclidean distances between group centroids (mean embedding vectors)
#' in a specified low-dimensional representation of a Seurat object.
#'
#' @details
#' For each group defined by \code{groupby}, the centroid is computed as the column-wise
#' mean of cell embeddings in \code{reduction}. Euclidean distances are then computed
#' between all pairs of group centroids.
#'
#' @param seurat_object A \code{Seurat} object.
#' @param groupby A character scalar giving the name of a metadata column in
#'   \code{seurat_object@meta.data} that defines group labels.
#' @param reduction A character scalar giving the name of a dimensional reduction
#'   (e.g., \code{"pca"}, \code{"harmony"}) used for distance calculations. This should be
#'   a geometry-preserving representation; visualization-oriented nonlinear embeddings
#'   (e.g., \code{"umap"}, \code{"tsne"}) are not recommended for distance computation.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' @param dims Optional integer vector selecting dimensions from \code{reduction}.
#'   If \code{NULL}, all dimensions are used.
#'
#' @return A symmetric \code{data.frame} (groups \eqn{\times} groups) of Euclidean distances
#'   between group centroids.
#'
#' @importFrom Seurat Embeddings
#' @importFrom stats dist
.eucldist <- function(seurat_object, groupby, reduction, verbose = TRUE,
                      dims = NULL) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("The first argument must be a Seurat object.")
  }
  if (!(reduction %in% names(seurat_object@reductions))) {
    if (verbose) {
      message("Reduction not found. Use a linear representation (e.g., PCA/Harmony), not UMAP/tSNE.")
    }
    stop("The specified reduction was not found in the Seurat object.")
  }

  labels <- seurat_object[[]][[groupby]]
  if (is.null(labels)) stop("groupby column not found in Seurat metadata.")
  labels <- as.character(labels)

  # drop NA groups
  keep <- !is.na(labels)
  labels <- labels[keep]

  emb <- Seurat::Embeddings(seurat_object, reduction = reduction)[keep, , drop = FALSE]
  emb <- .distance_select_dims(emb, dims)

  groups <- sort(unique(labels))

  means <- sapply(groups, function(g) {
    colMeans(emb[labels == g, , drop = FALSE])
  })

  df <- as.data.frame(as.matrix(stats::dist(t(means), method = "euclidean")))
  rownames(df) <- colnames(df) <- groups

  if (verbose) message("Computed Euclidean distances between group means.")
  df
}






#' @title spearmandist
#'
#' @description
#' Computes pairwise Spearman distances between group centroids (mean embedding vectors)
#' in a Seurat object.
#'
#' @details
#' Cells are grouped according to \code{groupby}. For each group, a centroid is computed
#' as the column-wise mean of cell embeddings in the specified \code{reduction} space.
#' Spearman's rank correlation coefficient is then computed between all pairs of group
#' centroids and converted to a distance measure as \eqn{1 - \rho}.
#'
#' This metric captures monotonic similarity between group-level expression profiles and
#' is insensitive to absolute scale, but does not account for within-group heterogeneity
#' or higher-order distributional differences.
#'
#' @param seurat_object A \code{Seurat} object.
#' @param groupby A character scalar giving the name of a metadata column in
#'   \code{seurat_object@meta.data} that defines group labels.
#' @param reduction A character scalar giving the name of a dimensional reduction
#'   (e.g., \code{"pca"}, \code{"harmony"}) used for distance calculations. This should be
#'   a geometry-preserving representation; visualization-oriented nonlinear embeddings
#'   (e.g., \code{"umap"}, \code{"tsne"}) are not recommended for distance computation.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' @param dims Optional integer vector selecting dimensions from \code{reduction}.
#'   If \code{NULL}, all dimensions are used.
#'
#' @return A symmetric \code{data.frame} (groups \eqn{\times} groups) containing
#'   pairwise Spearman distances between group centroids.
#'
#' @importFrom Seurat Embeddings
#' @importFrom stats cor

.spearmandist <- function(seurat_object, groupby, reduction, verbose = TRUE,
                          dims = NULL) {

  if (!inherits(seurat_object, "Seurat")) {
    stop("The first argument must be a Seurat object.")
  }

  if (!(reduction %in% names(seurat_object@reductions))) {
    stop(
      sprintf(
        "Reduction '%s' not found. Use a geometry-preserving reduction ",
        reduction
      ),
      "(e.g., PCA or Harmony) rather than visualization embeddings such as UMAP or t-SNE."
    )
  }

  # extract labels
  labels_raw <- seurat_object[[]][[groupby]]
  if (is.null(labels_raw)) {
    stop("groupby column not found in Seurat metadata.")
  }

  labels <- as.character(labels_raw)

  # drop NA labels
  keep <- !is.na(labels)
  labels <- labels[keep]
  emb <- Seurat::Embeddings(seurat_object, reduction = reduction)[keep, , drop = FALSE]
  emb <- .distance_select_dims(emb, dims)

  groups <- sort(unique(labels))

  # compute group means
  means <- sapply(groups, function(g) {
    colMeans(emb[labels == g, , drop = FALSE])
  })

  # initialize distance matrix
  df <- matrix(NA_real_, length(groups), length(groups),
               dimnames = list(groups, groups))

  if (verbose) message("Computing Spearman distances between group means...")

  # pairwise Spearman distance
  for (i in seq_along(groups)) {
    df[i, i] <- 0
    for (j in (i + 1):length(groups)) {
      if (j > length(groups)) break

      rho <- suppressWarnings(
        stats::cor(means[, i], means[, j], method = "spearman")
      )

      d <- 1 - rho
      df[i, j] <- df[j, i] <- d
    }
  }

  as.data.frame(df)
}




#' @title ComputeDistance
#'
#' @description
#' Computes pairwise distances between groups (cell states) in a Seurat object
#' using one of several distance metrics.
#'
#' @details
#' Groups are defined by \code{groupby}. Distances are computed in a specified
#' low-dimensional representation (\code{reduction}) that preserves meaningful
#' transcriptional geometry (e.g., PCA or Harmony). Depending on \code{method},
#' distances are computed between group centroids (Euclidean, Spearman) or between
#' full cell distributions (Energy distance).
#'
#' Visualization-oriented nonlinear embeddings (e.g., UMAP or t-SNE) are not
#' recommended for distance-based calculations.
#'
#' @param seurat_object A \code{Seurat} object.
#' @param groupby A character scalar giving the name of a metadata column in
#'   \code{seurat_object@meta.data} that defines group labels.
#' @param reduction A character scalar giving the name of a dimensional reduction
#'   (e.g., \code{"pca"}, \code{"harmony"}) used for distance calculations.
#' @param method A character scalar specifying the distance metric to use.
#'   Supported options are \code{"euclidean"}, \code{"edist"} (Energy distance),
#'   and \code{"spearman"}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#' @param dims Optional integer vector selecting dimensions from \code{reduction}.
#'   If \code{NULL} (default), all stored dimensions are used. Use the same
#'   dimensions for matrices that will be compared directly.
#'
#' @return A symmetric \code{data.frame} (groups \eqn{\times} groups) containing
#'   pairwise distances between groups.
#'
#' @examples
#' # Euclidean distance between group centroids
#' ComputeDistance(seurat_object, groupby = "Diagnosis",
#'                 reduction = "pca", method = "euclidean")
#'
#' # Energy distance between group distributions
#' ComputeDistance(seurat_object, groupby = "Diagnosis",
#'                 reduction = "pca", method = "edist")
#'
#' # Spearman distance between group centroids
#' ComputeDistance(seurat_object, groupby = "Diagnosis",
#'                 reduction = "pca", method = "spearman")
#'
#' @export

ComputeDistance <- function(seurat_object, groupby, reduction,
                            method = c("edist", "euclidean", "spearman"),
                            verbose = TRUE,
                            dims = NULL) {

  method <- match.arg(method)

  if (method == "euclidean") {
    return(.eucldist(
      seurat_object, groupby, reduction, verbose = verbose, dims = dims
    ))
  }

  if (method == "spearman") {
    return(.spearmandist(
      seurat_object, groupby, reduction, verbose = verbose, dims = dims
    ))
  }

  if (method == "edist") {
    return(.edist(
      seurat_object, groupby, reduction, verbose = verbose, dims = dims
    ))
  }
}


.distance_select_dims <- function(embeddings, dims = NULL) {
  if (is.null(dims)) return(embeddings)

  if (!is.numeric(dims) || length(dims) == 0L || anyNA(dims) ||
      any(dims %% 1 != 0) || any(dims < 1) ||
      any(dims > ncol(embeddings)) || anyDuplicated(dims)) {
    stop(
      "`dims` must contain unique integer indices between 1 and ",
      ncol(embeddings), ".",
      call. = FALSE
    )
  }

  embeddings[, as.integer(dims), drop = FALSE]
}
