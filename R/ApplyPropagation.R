
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
#' 
#' @details 
#' The function calculates the initial difference (\code{delta}) between the perturbed 
#' and baseline expression matrices. Over \code{n_iters} iterations, this \code{delta} 
#' is multiplied by the gene-gene \code{network} matrix via dot product . 
#' This mathematical diffusion simulates how a change in a hub gene ripples out to its 
#' connected targets.
#' 
#' To maintain stability, the propagated signal is penalized at each step by the 
#' \code{delta_scale} parameter, and total expression is strictly floored at 0 to 
#' prevent negative read counts. The final matrix is rounded to represent valid 
#' count data.
#'
#' @return A matrix (or \code{dgCMatrix}) containing the fully updated expression 
#'   matrix representing the global cell state after the perturbation signal has 
#'   propagated through the network.
#' 
#' @import Seurat
#' @export
ApplyPropagation <- function(
    seurat_obj,
    exp,
    exp_per,
    network,
    perturb_dir = perturb_dir,
    n_iters = 3,
    delta_scale = 0.2
){

    # todo: some checks for the network?
    # check that the network has the right genes
    print('here')
    # experimental: 
    # network[network < 0.05] <- 0
    # network <- methods::as(network, "CsparseMatrix")

    # row-normalize 
    row_sums <- Matrix::rowSums(network)
    row_sums[row_sums == 0] <- 1
    network <- network / row_sums

    # compute the difference between the perturbation and the observed expression
    delta <- exp_per - exp
    delta <- methods::as(delta, "CsparseMatrix")

    # backups
    delta_orig <- delta 
    exp_orig <- exp 
    exp_per_orig <- exp_per

    # define the "ceiling" of gene expression
    max_obs <- apply(exp, 1, max)

    # what is the sign of the perturbation dir?
    perturb_sign <- sign(perturb_dir)

    # run the signaling processing step iteratively
    for(i in 1:n_iters){  

        # compute the dot product between the TOM coefficients and the exp matrix
        delta <- network %*% delta

        # penalize the delta, or else the values rapidly get too large
        delta <- delta * delta_scale

        # update the expression matrix:
        exp_update <- exp_per + delta # (delta * perturb_sign)
        exp_update[exp_update < 0] <- 0

        # apply biological constraints (round such that we still have integers)
        # exp_update <- round(exp_update)
        exp_update <- methods::as(exp_update, "CsparseMatrix")
        exp_update@x <- pmin(exp_update@x, max_obs[exp_update@i + 1])

        # update the delta
        delta <- exp_update - exp_per

    }
    
    # return the matrix with the signal propagation applied
    exp_prop <- exp_per + delta
    exp_prop[exp_prop < 0] <- 0
    exp_prop <- round(exp_prop)
    exp_prop
}
