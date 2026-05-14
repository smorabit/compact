
#' Propagate In-Silico Perturbation Signal Through a Gene Network
#'
#' Propagates a log-space perturbation delta through a gene-gene co-expression
#' or regulatory network without count-space floor constraints, ensuring
#' symmetric behavior between up- and down-regulation perturbations.
#'
#' @param log_obs_mod A genes-by-cells log-normalized expression matrix for the
#'   module (or network) genes at baseline. Typically a row subset of the full
#'   log-normalized observed expression matrix.
#' @param delta_log A genes-by-cells log-space delta matrix for the same genes.
#'   Only the directly-perturbed (hub) gene rows are non-zero at input; the
#'   propagation step fills in downstream gene rows.
#' @param network A gene-by-gene matrix (e.g., TOM or regulatory adjacency)
#'   defining how signal routes between genes.
#' @param n_iters Number of iterative propagation steps. Default is \code{3}.
#' @param delta_scale Dampening factor applied at each iteration to prevent
#'   signal explosion. Default is \code{0.2}.
#' @param row_normalize Logical. If \code{TRUE}, each row of the network is
#'   scaled to sum to 1 before propagation (weighted-average diffusion rather
#'   than summation). Default is \code{FALSE}.
#' @param prune_network Logical. If \code{TRUE}, zero out edges below the
#'   \code{prune_percentile} threshold before propagation. Default is \code{FALSE}.
#' @param prune_percentile Numeric in \eqn{(0,1)}. Percentile threshold for edge
#'   pruning when \code{prune_network = TRUE}. Default is \code{0.95}.
#'
#' @details
#' This function propagates perturbation signal in log-normalized expression
#' space. At each iteration the delta matrix is multiplied by the network and
#' damped:
#' \deqn{\delta^{(t+1)} = W \cdot \delta^{(t)} \times \texttt{delta\_scale}}
#' where \eqn{W} is the (optionally row-normalized) network matrix.
#' The intermediate log-simulated expression for the module genes is
#' \deqn{\texttt{log\_sim} = \texttt{log\_obs\_mod} + \delta^{(n\_iters)}}}
#' and the result is floored at zero before returning, since log-normalized
#' expression cannot be negative.
#'
#' \strong{Signal amplification warning}: when \code{row_normalize = FALSE}
#' (default), the per-iteration effective multiplier is
#' \eqn{\bar{r} \times \texttt{delta\_scale}}, where \eqn{\bar{r}} is the
#' mean network row sum. If this product exceeds 1 the signal amplifies rather
#' than dampens across iterations, which can produce very large deltas that are
#' subsequently clipped by the floor. A warning is issued when this condition
#' holds so that the user can consider reducing \code{delta\_scale} or enabling
#' \code{row\_normalize = TRUE}.
#'
#' @return A genes-by-cells matrix of log-normalized simulated expression for
#'   the module genes.
#'
#' @import Matrix
#' @importFrom methods as
#' @export
ApplyPropagation <- function(
    log_obs_mod,
    delta_log,
    network,
    n_iters          = 3,
    delta_scale      = 0.2,
    row_normalize    = FALSE,
    prune_network    = FALSE,
    prune_percentile = 0.95
){

    if(!is.numeric(n_iters) || n_iters < 1){
        stop("n_iters must be a positive integer (>= 1)")
    }

    if(prune_network){
        threshold <- quantile(as.vector(network), prune_percentile, na.rm = TRUE)
        network[network < threshold] <- 0
    }

    if(row_normalize){
        row_sums <- Matrix::rowSums(network)
        row_sums[row_sums == 0] <- 1
        network <- network / row_sums
    }

    network   <- methods::as(network,   "CsparseMatrix")
    delta_log <- methods::as(delta_log, "CsparseMatrix")

    # warn when the effective per-iteration multiplier exceeds 1 — the signal
    # will amplify rather than dampen, and the floor may clip large fractions
    # of the propagated delta. consider reducing delta_scale or enabling
    # row_normalize = TRUE.
    if (!row_normalize) {
        mean_row_sum <- mean(Matrix::rowSums(network))
        if (mean_row_sum * delta_scale > 1) {
            warning(sprintf(
                paste0(
                    "Signal amplification detected: mean network row sum (%.2f) x ",
                    "delta_scale (%.2f) = %.2f > 1. The propagated delta will grow ",
                    "across iterations rather than dampen, and large negative values ",
                    "will be clipped to zero by the expression floor. Consider ",
                    "reducing delta_scale or setting row_normalize = TRUE."
                ),
                mean_row_sum, delta_scale, mean_row_sum * delta_scale
            ), call. = FALSE)
        }
    }

    delta_initial <- delta_log

    for(i in seq_len(n_iters)){
        delta_log <- network %*% delta_log * delta_scale
    }

    CheckSignalDecay(
        delta_initial = delta_initial,
        delta_final   = delta_log,
        row_normalize = row_normalize,
        delta_scale   = delta_scale,
        n_iters       = n_iters
    )

    # floor at zero: log-normalized expression cannot be negative
    result <- log_obs_mod + delta_log
    result[result < 0] <- 0
    result
}
