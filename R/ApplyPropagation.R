
#' ApplyPropagation
#'
#' This function applies an in-silico perturbation to selected features in a Seurat object.
#' 
#' @return A dgCMatrix object containing the updated expression matrix with the applied perturbations
#' 
#' @param seurat_obj A Seurat object
#' @param exp A features by cells matrix containing the observed expression matrix.
#' @param exp_per A features by cells matrix containing the perturbation results from ApplyPerturbation.
#' @param network A gene-gene network to apply the signal propagation.
#' @param n_iters The number of times to apply the signal propagation.
#' 
#' @details 
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

    # compute the difference between the perturbation and the observed expression
    delta <- exp_per - exp

    # backups
    delta_orig <- delta 
    exp_orig <- exp 
    exp_per_orig <- exp_per

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

        exp_update <- round(exp_update)

        # update the delta
        delta <- exp_update - exp_per

    }
    
    # return the matrix with the signal propagation applied
    exp_prop <- exp_per + delta
    exp_prop[exp_prop < 0] <- 0
    exp_prop
}
