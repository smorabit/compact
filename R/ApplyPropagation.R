
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
#' space, treating the directly-perturbed (hub) genes as a \strong{persistent
#' signal source} rather than a one-time impulse. This models a constitutive
#' perturbation (e.g. a stable CRISPR KO/OE or a sustained stimulus), in which
#' the hub genes are held at their new expression level indefinitely while
#' downstream genes respond iteratively. Hub genes are identified as the rows
#' of \code{delta_log} that are non-zero at input.
#'
#' At each iteration the delta matrix is multiplied by the network and damped,
#' then hub gene rows are reset to their initial value:
#' \deqn{\delta^{(t)} = W \cdot \delta^{(t-1)} \times \texttt{delta\_scale},
#'   \quad \text{then} \quad \delta^{(t)}_{\mathcal{H}} \leftarrow
#'   \delta^{(0)}_{\mathcal{H}}}
#' where \eqn{W} is the (optionally row-normalized) network matrix and
#' \eqn{\mathcal{H}} is the set of hub gene rows. Without this reset, hub
#' genes would be overwritten by the network multiplication like any other
#' gene, and their signal would weaken with each iteration rather than persist.
#' The intermediate log-simulated expression for the module genes is
#' \deqn{\texttt{log\_sim} = \texttt{log\_obs\_mod} + \delta^{(n\_iters)}}
#' and the result is floored at zero before returning, since log-normalized
#' expression cannot be negative.
#'
#' \strong{Signal amplification warning}: when \code{row_normalize = FALSE}
#' (default), the per-iteration effective multiplier for the non-hub
#' (downstream) genes is \eqn{\bar{r} \times \texttt{delta\_scale}}, where
#' \eqn{\bar{r}} is the mean row sum of the network restricted to non-hub
#' rows. If this product exceeds 1 the downstream signal amplifies rather
#' than dampens across iterations, which can produce very large deltas that
#' are subsequently clipped by the floor. A warning is issued when this
#' condition holds so that the user can consider reducing \code{delta\_scale}
#' or enabling \code{row\_normalize = TRUE}. Hub rows are excluded from this
#' check because they are pinned every iteration and cannot amplify.
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

    # snapshot the primary perturbation and identify hub (directly-perturbed)
    # genes as rows with a non-zero initial delta. hub genes are treated as a
    # persistent signal source: they are pinned back to delta_initial after
    # every iteration below, so only non-hub rows evolve or can amplify.
    delta_initial <- delta_log
    hub_mask <- Matrix::rowSums(abs(delta_initial)) > 0

    # warn when the effective per-iteration multiplier exceeds 1 — the signal
    # will amplify rather than dampen, and the floor may clip large fractions
    # of the propagated delta. consider reducing delta_scale or enabling
    # row_normalize = TRUE. only the non-hub subgraph is checked: hub rows
    # are pinned every iteration and cannot contribute to runaway amplification.
    if (!row_normalize && any(!hub_mask)) {
        mean_row_sum <- mean(Matrix::rowSums(network[!hub_mask, , drop = FALSE]))
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

    # iterative diffusion with a persistent hub source: hub genes model a
    # constitutive perturbation (e.g. stable CRISPR KO/OE) that holds them at
    # their new expression level throughout propagation, rather than a
    # one-time impulse that would otherwise get diffused away and replaced by
    # a weak, washed-out echo of itself. after each network multiplication,
    # hub rows are reset to their original delta so they continue
    # broadcasting at full strength while downstream (non-hub) genes
    # accumulate contributions from progressively longer network paths.
    for(i in seq_len(n_iters)){
        delta_log <- network %*% delta_log * delta_scale
        delta_log[hub_mask, ] <- delta_initial[hub_mask, ]
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
