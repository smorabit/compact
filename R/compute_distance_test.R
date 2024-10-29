# Description: Computes E-test statistics and Compute Distance for each group in a Seurat object

library(Seurat)
library(dplyr)
library(rdist)
library(energy)


#' @title edist
#'
#' @description Computes pairwise E-distances on a Seurat object.
#'     Computes E-distance between all groups in a Seurat object in space
#'     given by reduction.
#' @note This function is for internal use and is not exported. Computes pairwise E-distances on a Seurat object.
#'     Computes E-distance between all groups in a Seurat object in space
#'     given by reduction.
#' @param seurat_object An object of class Seurat.
#' @param groupby An object of class character. Points to the column in the
#'     Seurat object's meta data that contains the group labels.
#' @param reduction An object of class character. The reduction / embedding in
#'     seurat_object that is used to compute the E-distance in. Can be "pca", "harmony", or any linear dimensionality reduction.
#' @param sample_correction An object of class logical. If TRUE, the
#'     E-distances are corrected for sample size. Will make it not a proper
#'     distance, leads to negative values.
#' @param verbose An object of class logical. If TRUE, prints messages.
#'    Default is TRUE.
#' @return Returns an object of class data.frame. For each group contains the
#'     E-test p-value and the E-distance to control group.
#' @examples
#'     # Add some code illustrating how to use the function
#' @importFrom Seurat Embeddings VariableFeatures RunPCA
#' @importFrom rdist cdist pdist
#' @importFrom dplyr select
.edist <- function(seurat_object, groupby, reduction,
                  sample_correction = FALSE, verbose = TRUE) {
    if (!inherits(seurat_object, "Seurat")) {
        stop("The first argument must be a Seurat object.")
    }
    if (!(reduction %in% names(seurat_object@reductions))) {
        if (verbose) {
            message("The specified reduction was not found in the Seurat object. Ensure that the applied dimensional reduction method is linear in nature, such as PCA or Harmony on the selected linear model.")
        }
        stop("The specified reduction was not found in the Seurat object.")
    }
    labels <- seurat_object[[]][[groupby]]
    groups <- unique(labels)
    emb <- Seurat::Embeddings(seurat_object, reduction = reduction)

    df <- setNames(data.frame(matrix(ncol = length(groups),
                                     nrow = length(groups)),
                              row.names = groups),
                   groups)
    if (verbose) {
        print("Computing E-test statistics for each group.")
    }
    completed_groups <- c()
    for (groupx in groups) {
        for (groupy in groups){
            if (groupy %in% completed_groups) {
                next
            }  # skip if already computed
            x <- as.matrix(emb)[labels == groupx, ]
            y <- as.matrix(emb)[labels == groupy, ]

            N <- nrow(x)
            M <- nrow(y)

            dist_xy <- rdist::cdist(x,y)
            dist_x <- rdist::pdist(x)
            dist_y <- rdist::pdist(y)

            if (sample_correction) {
                ed <- 2 * (sum(dist_xy) / (N * M)) - (sum(dist_x) / (N * (N - 1))) - (sum(dist_y) / (M * (M - 1)))
            } else {
                ed <- 2 * (sum(dist_xy) / (N * M)) - (sum(dist_x) / (N * N)) - (sum(dist_y) / (M * M))
            }

            df[groupx, groupy] <- df[groupy, groupx] <- ed
        }
        completed_groups <- c(completed_groups, groupx)
    }
    return(df)
}

#' @title eucldist
#' @description Computes the Euclidean distance between the means of two groups in a Seurat object.
#' @note This function is for internal use and is not exported. Computes the Euclidean distance between the means of two groups in a Seurat object.
#' @param seurat_object An object of class Seurat.
#' @param groupby An object of class character. Points to the column in the
#' Seurat object's meta data that contains the group labels.
#' @param reduction An object of class character. The reduction / embedding in
#' seurat_object that is used to compute the means. Can be "pca", "harmony", or any linear dimensionality reduction.
#' @param verbose An object of class logical. If TRUE, prints messages.
#' Default is TRUE.
#' @return Returns an object of class data.frame containing the Euclidean distance between each pair of groups.
#' @importFrom Seurat Embeddings
#' @importFrom stats dist
.eucldist <- function(seurat_object, groupby, reduction, verbose = TRUE) {
    if (!inherits(seurat_object, "Seurat")) {
        stop("The first argument must be a Seurat object.")
    }
    if (!(reduction %in% names(seurat_object@reductions))) {
        if (verbose) {
            message("The specified reduction was not found in the Seurat object. Ensure that the applied dimensional reduction method is linear in nature, such as PCA or Harmony on the selected linear model.")
        }
        stop("The specified reduction was not found in the Seurat object.")
    }
    labels <- seurat_object[[]][[groupby]]
    groups <- unique(labels)
    emb <- Seurat::Embeddings(seurat_object, reduction = reduction)

    means <- sapply(groups, function(group) {
        colMeans(emb[labels == group, ])
    })

    df <- as.data.frame(as.matrix(stats::dist(t(means))))
    rownames(df) <- colnames(df) <- groups

    if (verbose) {
        print("Computed Euclidean distances between group means.")
    }
    return(df)
}

#' @title ComputeDistance
#' @description Computes either the Euclidean or E-distance between groups in a Seurat object, depending on the chosen method.
#' @param seurat_object An object of class Seurat.
#' @param groupby An object of class character. Points to the column in the
#' Seurat object's meta data that contains the group labels.
#' @param reduction An object of class character. The reduction / embedding in
#' seurat_object that is used to compute the distance. Can be "pca", "harmony", or any linear dimensionality reduction.
#' @param method An object of class character. The distance calculation method, either "euclidean" or "edist".
#' @param verbose An object of class logical. If TRUE, prints messages.
#' Default is TRUE.
#' @return Returns an object of class data.frame containing the distance between each pair of groups.
#' @examples
#' compdist <- function(seurat_object, groupby = "Diagnosis", reduction = "pca", method = "euclidean", verbose = TRUE)
#' @export
ComputeDistance <- function(seurat_object, groupby, reduction, method = "edist", verbose = TRUE) {

    if (method == "euclidean") {
    if (verbose) {
        message("Computing euclidean distance.")
    }
        return(.eucldist(seurat_object, groupby, reduction, verbose))
    } else if (method == "edist") {
    if (verbose) {
        message("Computing edist distance.")
    }
        return(.edist(seurat_object, groupby, reduction, verbose = verbose))
    } else {
        stop("Invalid method. Choose either 'euclidean' or 'edist'.")
    }
}
