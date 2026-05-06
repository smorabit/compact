#' Calculate Log2 Fold Change of In-Silico Perturbation
#'
#' @description
#' Computes the log2 fold change (Log2FC), mean delta, and average baseline
#' expression for genes following an in-silico perturbation. It compares the
#' perturbed assay against the baseline observed expression.
#'
#' @param seurat_obj A Seurat object containing both the baseline and perturbed data.
#' @param perturbation_name Character. The name of the assay containing the perturbed expression data.
#' @param assay_obs Character. The name of the baseline assay. Default is \code{'RNA'}.
#' @param module Character. Optional. The specific hdWGCNA module name to analyze. If provided, the output will include a column indicating whether each gene is a top hub.
#' @param n_hubs Integer. The number of top hub genes to flag if \code{module} is provided. Default is \code{10}.
#' @param pseudocount Numeric. Pseudocount added to average expression before log2 transformation to prevent infinite values. Default is \code{1}.
#' @param layer_norm Character. The layer/layer containing normalized data. Default is \code{'data'}.
#' @param layer_counts Character. The layer/layer containing raw counts. Default is \code{'counts'}.
#' @param wgcna_name Character. The name of the hdWGCNA experiment in \code{seurat_obj@misc}. If NULL, uses the active WGCNA experiment.
#'
#' @return A data.frame containing gene names, mean delta (raw counts difference), average baseline expression, Log2FC, and optionally module hub status.
#'
#' @import Seurat
#' @importFrom Matrix rowMeans rowSums
#' @importFrom dplyr left_join rename mutate
#' @export
PerturbationLog2FC <- function(
    seurat_obj,
    perturbation_name,
    assay_obs = 'RNA',
    module = NULL,
    n_hubs = 10,
    pseudocount = 1,
    layer_norm = 'data',
    layer_counts = 'counts',
    wgcna_name = NULL
) {

    if(is.null(wgcna_name)){ wgcna_name <- seurat_obj@misc$active_wgcna }

    # validate assays
    if(!(perturbation_name %in% names(seurat_obj@assays))){
        stop(paste0("Perturbation assay '", perturbation_name, "' not found in seurat_obj. Available assays: ",
                    paste(names(seurat_obj@assays), collapse = ', ')))
    }
    if(!(assay_obs %in% names(seurat_obj@assays))){
        stop(paste0("Baseline assay '", assay_obs, "' not found in seurat_obj. Available assays: ",
                    paste(names(seurat_obj@assays), collapse = ', ')))
    }

    # 1. Determine features to analyze
    if (!is.null(module)) {

        # get all genes in this module; use direct indexing to avoid subset()
        # scoping bug where `module == module` compares the column to itself
        hub_df <- GetHubGenes(seurat_obj, n_hubs = Inf, wgcna_name = wgcna_name)
        hub_df <- hub_df[hub_df$module == module, ]

        if(nrow(hub_df) == 0){
            stop(paste0("Module '", module, "' not found. Check available modules with GetModules()."))
        }

        top_hubs <- GetHubGenes(seurat_obj, n_hubs = n_hubs, wgcna_name = wgcna_name)
        top_hubs <- top_hubs[top_hubs$module == module, "gene_name"]

        hub_df$hub <- ifelse(hub_df$gene_name %in% top_hubs, 'hub', 'other')
        features <- hub_df$gene_name
    } else {
        # if no module provided, use all genes present in the perturbed assay
        features <- rownames(seurat_obj[[perturbation_name]])
    }

    # 2. Extract matrices
    X_counts <- GetAssayData(seurat_obj, assay = assay_obs, layer = layer_counts)[features, ]
    X_norm <- GetAssayData(seurat_obj, assay = assay_obs, layer = layer_norm)[features, ]

    X_per_counts <- GetAssayData(seurat_obj, assay = perturbation_name, layer = layer_counts)[features, ]
    X_per_norm <- GetAssayData(seurat_obj, assay = perturbation_name, layer = layer_norm)[features, ]

    # 3. Fast Matrix calculations
    delta_matrix <- X_per_counts - X_counts
    mean_delta <- Matrix::rowMeans(delta_matrix)
    avg_exp <- Matrix::rowMeans(X_norm)
    avg_per_exp <- Matrix::rowMeans(X_per_norm)

    # 4. Calculate Log2FC
    l2fc <- log2((avg_per_exp + pseudocount) / (avg_exp + pseudocount))

    # 5. Construct final dataframe
    plot_df <- data.frame(
        gene_name = features,
        mean_delta = mean_delta,
        avg_exp = avg_exp,
        log2fc = l2fc,
        stringsAsFactors = FALSE
    )

    # merge hub info if a module was specified
    if (!is.null(module)) {
        plot_df <- dplyr::left_join(plot_df, hub_df, by = "gene_name")
    }

    return(plot_df)
}

#' Check for Signal Saturation After Network Propagation
#'
#' After running \code{ApplyPropagation}, checks whether the propagated
#' expression has saturated at the biological floor (0) or ceiling (max observed
#' expression). Heavy saturation can cause up- and down-regulation perturbations
#' to produce indistinguishable transition vector fields.
#'
#' @param exp_obs A features-by-cells matrix of baseline (unperturbed) expression counts.
#' @param exp_prop A features-by-cells matrix of the final propagated expression counts
#'   (output of \code{ApplyPropagation}).
#' @param delta_scale Numeric. The \code{delta_scale} used during propagation; reported in the warning.
#' @param n_iters Integer. The \code{n_iters} used during propagation; reported in the warning.
#' @param apply_ceiling Logical. Whether a ceiling constraint was applied during propagation.
#'   Default \code{FALSE}.
#' @param max_obs Numeric vector. Per-gene ceiling values (length equal to \code{nrow(exp_prop)}).
#'   Only used when \code{apply_ceiling = TRUE}.
#' @param saturation_thresh Numeric. Fraction of perturbed gene-cell entries that must be
#'   saturated to trigger a warning. Default \code{0.5}.
#'
#' @return Invisibly returns a named list with \code{floor_frac} and, if
#'   \code{apply_ceiling = TRUE}, \code{ceil_frac} — the computed saturation fractions.
#'
#' @details
#' Perturbed genes are identified as rows where the sum of \code{|exp_prop - exp_obs|}
#' is non-zero. Among those rows, floor saturation is the fraction of gene-cell entries
#' where propagated expression equals 0. If \code{apply_ceiling = TRUE}, ceiling saturation
#' is the fraction of entries at or above \code{max_obs}. A warning is issued for each
#' fraction exceeding \code{saturation_thresh}.
#'
#' @export
CheckSaturation <- function(
    exp_obs,
    exp_prop,
    delta_scale,
    n_iters,
    apply_ceiling = FALSE,
    max_obs = NULL,
    saturation_thresh = 0.5
) {
    perturbed_genes <- which(Matrix::rowSums(abs(exp_prop - exp_obs)) > 0)
    result <- list()

    if (length(perturbed_genes) == 0) {
        return(invisible(result))
    }

    prop_sub <- as.matrix(exp_prop[perturbed_genes, , drop = FALSE])

    # floor saturation: entries driven to expression = 0
    floor_frac <- mean(prop_sub == 0)
    result$floor_frac <- floor_frac

    if (floor_frac > saturation_thresh) {
        warning(sprintf(
            "%.0f%% of perturbed gene-cell entries saturated at the expression floor (= 0) after propagation. Up- and down-regulation perturbations may produce indistinguishable results. Consider reducing `delta_scale` (current: %.2f) or `n_iters` (current: %d), or setting `row_normalize = TRUE`.",
            floor_frac * 100, delta_scale, n_iters
        ), call. = FALSE)
    }

    # ceiling saturation: entries driven to max_obs (only when apply_ceiling = TRUE)
    if (apply_ceiling && !is.null(max_obs)) {
        ceil_check <- sweep(prop_sub, 1, max_obs[perturbed_genes], ">=")
        ceil_frac <- mean(ceil_check)
        result$ceil_frac <- ceil_frac

        if (ceil_frac > saturation_thresh) {
            warning(sprintf(
                "%.0f%% of perturbed gene-cell entries saturated at the expression ceiling after propagation. Up- and down-regulation perturbations may produce indistinguishable results. Consider reducing `delta_scale` (current: %.2f) or `n_iters` (current: %d), or setting `row_normalize = TRUE`.",
                ceil_frac * 100, delta_scale, n_iters
            ), call. = FALSE)
        }
    }

    invisible(result)
}

#' Check for Propagation Signal Decay
#'
#' After running \code{ApplyPropagation}, checks whether the propagated signal
#' reaching downstream (non-directly-perturbed) genes has decayed below a
#' detectable threshold. With \code{row_normalize = TRUE} and
#' \code{delta_scale <= 1}, the row-normalized network acts as a stochastic
#' matrix — repeated multiplication can only redistribute signal, never amplify
#' it — so the propagated delta decays geometrically toward zero and is often
#' lost to integer rounding.
#'
#' @param delta_initial A features-by-cells matrix of the initial perturbation
#'   delta (\code{exp_per - exp}) before any propagation iterations.
#' @param delta_final A features-by-cells matrix of the propagated delta after
#'   \code{n_iters} iterations (i.e., the delta added to \code{exp_per} to
#'   produce the final propagated expression).
#' @param row_normalize Logical. Whether the network was row-normalized during
#'   propagation.
#' @param delta_scale Numeric. The \code{delta_scale} used during propagation;
#'   reported in the warning.
#' @param n_iters Integer. The \code{n_iters} used during propagation; reported
#'   in the warning.
#' @param min_mean_delta Numeric. Minimum mean absolute delta in downstream
#'   genes considered a detectable signal. Defaults to \code{0.001}, which is
#'   appropriate for log-normalized expression space.
#'
#' @return Invisibly returns a named list with:
#'   \describe{
#'     \item{\code{mean_delta_initial}}{Mean absolute delta across directly
#'       perturbed gene-cell entries.}
#'     \item{\code{mean_delta_downstream}}{Mean absolute delta across downstream
#'       gene-cell entries after propagation.}
#'   }
#'
#' @details
#' Directly perturbed genes are identified as rows where the sum of
#' \code{|delta_initial|} is non-zero. All remaining rows are treated as
#' downstream. A warning is issued when the mean absolute downstream delta
#' is below \code{min_mean_delta}, with specific guidance when
#' \code{row_normalize = TRUE} and \code{delta_scale <= 1} — the parameter
#' combination that guarantees decay regardless of \code{n_iters}.
#'
#' @export
CheckSignalDecay <- function(
    delta_initial,
    delta_final,
    row_normalize,
    delta_scale,
    n_iters,
    min_mean_delta = 0.001
) {
    result <- list()

    # identify directly perturbed genes (non-zero initial delta) and downstream genes
    directly_perturbed <- which(Matrix::rowSums(abs(delta_initial)) > 0)
    downstream_genes   <- setdiff(seq_len(nrow(delta_initial)), directly_perturbed)

    if (length(downstream_genes) == 0 || length(directly_perturbed) == 0) {
        return(invisible(result))
    }

    mean_delta_initial    <- mean(abs(as.matrix(delta_initial[directly_perturbed, , drop = FALSE])))
    mean_delta_downstream <- mean(abs(as.matrix(delta_final[downstream_genes,   , drop = FALSE])))

    result$mean_delta_initial    <- mean_delta_initial
    result$mean_delta_downstream <- mean_delta_downstream

    # nothing to check if the initial perturbation itself was zero
    if (mean_delta_initial == 0) {
        return(invisible(result))
    }

    if (mean_delta_downstream < min_mean_delta) {
        if (row_normalize && delta_scale <= 1) {
            warning(sprintf(
                paste0(
                    "Mean absolute delta in downstream genes (%.6f) is below the detectable ",
                    "threshold (%.4f). With row_normalize = TRUE, a delta_scale <= 1 ",
                    "(current: %.2f) guarantees signal decay over iterations because the ",
                    "row-normalized network can only redistribute signal, not amplify it. ",
                    "Consider setting row_normalize = FALSE to preserve propagation signal."
                ),
                mean_delta_downstream, min_mean_delta, delta_scale
            ), call. = FALSE)
        } else {
            warning(sprintf(
                paste0(
                    "Mean absolute delta in downstream genes (%.6f) is below the detectable ",
                    "threshold (%.4f). Consider increasing delta_scale (current: %.2f) or ",
                    "n_iters (current: %d)."
                ),
                mean_delta_downstream, min_mean_delta, delta_scale, n_iters
            ), call. = FALSE)
        }
    }

    invisible(result)
}