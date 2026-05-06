
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
#' space. Unlike the previous count-space implementation, there is no floor at
#' zero and no integer rounding, so down-regulation produces negative
#' log-deltas in downstream genes for \emph{all} cells — including
#' zero-expressing cells that would have been clipped to zero in count space.
#' This symmetry between up- and down-regulation is essential for
#' \code{PerturbationTransitions} to produce correctly opposite vector fields.
#'
#' At each iteration the delta matrix is multiplied by the network and damped:
#' \deqn{\delta^{(t+1)} = W \cdot \delta^{(t)} \times \texttt{delta\_scale}}
#' where \eqn{W} is the (optionally row-normalized) network matrix.
#' The final log-simulated expression for the module genes is
#' \deqn{\texttt{log\_sim} = \texttt{log\_obs\_mod} + \delta^{(n\_iters)}}
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

    log_obs_mod + delta_log
}
