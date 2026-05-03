
#' Propagate In-Silico Perturbation Signal Through a Gene Network
#'
#' This function models the downstream secondary effects of an in-silico perturbation 
#' by propagating the initial expression changes through a gene-gene co-expression 
#' or regulatory network. 
#' 
#' @param seurat_obj A Seurat object containing the dataset.
#' @param exp A features-by-cells matrix containing the baseline (unperturbed) observed expression data.
#' @param exp_per A features-by-cells matrix containing the initial perturbation results (e.g., output from \code{ApplyPerturbation}).
#' @param network A gene-by-gene matrix representing the network structure (e.g., a Topological Overlap Matrix or gene regulatory network) used to route the signal.
#' @param perturb_dir Numeric. The direction and magnitude of the original perturbation, used to extract the sign for propagation scaling.
#' @param n_iters Numeric/Integer. The number of iterative steps to apply the signal propagation. Default is \code{3}.
#' @param delta_scale Numeric. A penalization/dampening factor applied to the delta matrix at each iteration to prevent biologically unrealistic exponential explosions in expression values. Default is \code{0.2}.
#' @param row_normalize Logical. If \code{TRUE}, normalizes the network such that each row sums to 1. This converts the matrix multiplication from a cumulative sum into a weighted average, helping to prevent signal explosion in highly connected background genes. Default is \code{FALSE}.
#' @param apply_ceiling Logical. If \code{TRUE}, constrains the maximum propagated expression for each gene so it does not exceed its maximum observed value in the baseline data (adjusted by the \code{ceiling_multiplier}). Default is \code{FALSE}.
#' @param ceiling_multiplier Numeric. A scaling factor applied to the expression ceiling if \code{apply_ceiling = TRUE}. For example, \code{1.10} allows simulated expression to reach 10\% higher than the absolute maximum observed baseline value. Default is \code{1.0}.
#' @param prune_network Logical. If \code{TRUE}, removes weak background edges from the network prior to propagation based on a calculated percentile threshold. Highly recommended for dense networks like WGCNA TOMs. Default is \code{FALSE}.
#' @param prune_percentile Numeric. The percentile threshold (between 0 and 1) used to prune the network if \code{prune_network = TRUE}. For example, \code{0.95} drops the weakest 95\% of edges and retains only the top 5\% strongest connections. Default is \code{0.95}.
#' 
#' @details 
#' The function calculates the initial difference (\code{delta}) between the perturbed 
#' and baseline expression matrices. Over \code{n_iters} iterations, this \code{delta} 
#' is multiplied by the gene-gene \code{network} matrix via dot product. 
#' This mathematical diffusion simulates how a change in a hub gene ripples out to its 
#' connected targets.
#' 
#' To maintain stability and biological realism, the propagated signal is penalized at each step by the 
#' \code{delta_scale} parameter. Users can further stabilize the diffusion by pruning weak network edges 
#' (\code{prune_network}), row-normalizing the network (\code{row_normalize}), and enforcing a strict 
#' biological ceiling (\code{apply_ceiling}) to prevent genes from entering biologically impossible 
#' high-expression states. Finally, total expression is strictly floored at 0 to 
#' prevent negative read counts, and the final matrix is rounded to represent valid count data.
#'
#' @return A matrix (or \code{dgCMatrix}) containing the fully updated expression 
#'   matrix representing the global cell state after the perturbation signal has 
#'   propagated through the network.
#' 
#' @import Seurat
#' @importFrom Matrix rowSums
#' @export
#' 
ApplyPropagation <- function(
    seurat_obj,
    exp,
    exp_per,
    network,
    perturb_dir,
    n_iters = 3,
    delta_scale = 0.2,
    row_normalize = FALSE,
    apply_ceiling = FALSE,
    ceiling_multiplier = 1.0,
    prune_network = FALSE,
    prune_percentile = 0.95
){

    # n_iters must be at least 1; seq_len(0) would silently skip the loop but leave
    # delta as (exp_per - exp), producing exp_prop = 2*exp_per - exp, which is wrong
    if(!is.numeric(n_iters) || n_iters < 1){
        stop("n_iters must be a positive integer (>= 1)")
    }

    # note: perturb_dir is accepted for API compatibility with callers but is not
    # used in the propagation math; the sign of the effect is carried by exp_per.

    # 1. OPTIONAL: Prune the network based on edge weight percentiles
    if(prune_network){
        # Calculate the threshold value at the requested percentile
        threshold <- quantile(as.vector(network), prune_percentile, na.rm = TRUE)
        network[network < threshold] <- 0
    }
    
    # 2. OPTIONAL: Row-normalize the network to prevent summation explosion
    if(row_normalize){
        row_sums <- Matrix::rowSums(network)
        row_sums[row_sums == 0] <- 1
        network <- network / row_sums
    }

    # Ensure network is a sparse matrix for fast multiplication
    network <- methods::as(network, "CsparseMatrix")

    # Compute the initial difference between perturbation and observed expression
    delta <- exp_per - exp
    delta <- methods::as(delta, "CsparseMatrix")

    # 3. OPTIONAL: Define the biological ceiling with a dynamic multiplier
    if(apply_ceiling){
        # multiply the max observed by the buffer (e.g., 1.05 for a 5% buffer)
        max_obs <- apply(exp, 1, max) * ceiling_multiplier
    }

    # Run the signaling processing step iteratively
    for(i in seq_len(n_iters)){

        # Compute the dot product between the TOM coefficients and the exp matrix
        delta <- network %*% delta

        # Penalize the delta, or else the values rapidly get too large
        delta <- delta * delta_scale

        # Update the expression matrix
        exp_update <- exp_per + delta
        exp_update[exp_update < 0] <- 0

        # Ensure sparse format for memory efficiency and for the ceiling subsetting below
        exp_update <- methods::as(exp_update, "CsparseMatrix")

        # 4. OPTIONAL: Apply biological constraints to the non-zero values
        if(apply_ceiling){
            exp_update@x <- pmin(exp_update@x, max_obs[exp_update@i + 1])
        }

        # Update the delta for the next iteration
        delta <- exp_update - exp_per

    }
    
    # Return the matrix with the signal propagation applied
    exp_prop <- exp_per + delta
    
    # Floor negative values to 0 and round to maintain integer count structure
    exp_prop[exp_prop < 0] <- 0
    exp_prop <- round(exp_prop)
    
    return(exp_prop)
}






# ApplyPropagation <- function(
#     seurat_obj,
#     exp,
#     exp_per,
#     network,
#     perturb_dir = perturb_dir,
#     n_iters = 3,
#     delta_scale = 0.2
# ){

#     # todo: some checks for the network?
#     # check that the network has the right genes

#     # experimental: 
#     # network[network < 0.05] <- 0
#     # network <- methods::as(network, "CsparseMatrix")

#     # row-normalize 
#     row_sums <- Matrix::rowSums(network)
#     row_sums[row_sums == 0] <- 1
#     network <- network / row_sums

#     # compute the difference between the perturbation and the observed expression
#     delta <- exp_per - exp
#     delta <- methods::as(delta, "CsparseMatrix")

#     # backups
#     delta_orig <- delta 
#     exp_orig <- exp 
#     exp_per_orig <- exp_per

#     # define the "ceiling" of gene expression
#     max_obs <- apply(exp, 1, max)

#     # what is the sign of the perturbation dir?
#     perturb_sign <- sign(perturb_dir)

#     # run the signaling processing step iteratively
#     for(i in 1:n_iters){  

#         # compute the dot product between the TOM coefficients and the exp matrix
#         delta <- network %*% delta

#         # penalize the delta, or else the values rapidly get too large
#         delta <- delta * delta_scale

#         # update the expression matrix:
#         exp_update <- exp_per + delta # (delta * perturb_sign)
#         exp_update[exp_update < 0] <- 0

#         # apply biological constraints (round such that we still have integers)
#         # exp_update <- round(exp_update)
#         exp_update <- methods::as(exp_update, "CsparseMatrix")
#         exp_update@x <- pmin(exp_update@x, max_obs[exp_update@i + 1])

#         # update the delta
#         delta <- exp_update - exp_per

#     }
    
#     # return the matrix with the signal propagation applied
#     exp_prop <- exp_per + delta
#     exp_prop[exp_prop < 0] <- 0
#     exp_prop <- round(exp_prop)
#     exp_prop
# }
