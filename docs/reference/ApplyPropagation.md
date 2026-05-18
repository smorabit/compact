# Propagate In-Silico Perturbation Signal Through a Gene Network

Propagates a log-space perturbation delta through a gene-gene
co-expression or regulatory network without count-space floor
constraints, ensuring symmetric behavior between up- and down-regulation
perturbations.

## Usage

``` r
ApplyPropagation(
  log_obs_mod,
  delta_log,
  network,
  n_iters = 3,
  delta_scale = 0.2,
  row_normalize = FALSE,
  prune_network = FALSE,
  prune_percentile = 0.95
)
```

## Arguments

- log_obs_mod:

  A genes-by-cells log-normalized expression matrix for the module (or
  network) genes at baseline. Typically a row subset of the full
  log-normalized observed expression matrix.

- delta_log:

  A genes-by-cells log-space delta matrix for the same genes. Only the
  directly-perturbed (hub) gene rows are non-zero at input; the
  propagation step fills in downstream gene rows.

- network:

  A gene-by-gene matrix (e.g., TOM or regulatory adjacency) defining how
  signal routes between genes.

- n_iters:

  Number of iterative propagation steps. Default is `3`.

- delta_scale:

  Dampening factor applied at each iteration to prevent signal
  explosion. Default is `0.2`.

- row_normalize:

  Logical. If `TRUE`, each row of the network is scaled to sum to 1
  before propagation (weighted-average diffusion rather than summation).
  Default is `FALSE`.

- prune_network:

  Logical. If `TRUE`, zero out edges below the `prune_percentile`
  threshold before propagation. Default is `FALSE`.

- prune_percentile:

  Numeric in \\(0,1)\\. Percentile threshold for edge pruning when
  `prune_network = TRUE`. Default is `0.95`.

## Value

A genes-by-cells matrix of log-normalized simulated expression for the
module genes.

## Details
