# colDeltaCor_all

This function calculates cell-to-cell correlations based on observed and
perturbed gene expression matrices, with an option to limit computations
to pairs of cells connected in a KNN graph. Correlations are computed in
parallel to optimize performance.

## Usage

``` r
colDeltaCor_knn(
  observed_matrix,
  delta_matrix,
  knn_graph = NULL,
  n_cores = detectCores() - 1
)
```

## Arguments

- observed_matrix:

  A matrix of observed gene expression values (genes by cells).

- delta_matrix:

  A matrix of perturbed gene expression values (genes by cells).

- knn_graph:

  Optional; an adjacency matrix representing the k-nearest-neighbor
  (KNN) graph. If provided, only correlations for pairs of cells
  connected in the KNN graph will be computed.

- n_cores:

  The number of threads to use for parallel processing. Defaults to
  `detectCores() - 1`.

## Value

A matrix of cell-to-cell correlation values.

## Details

`colDeltaCor_all` calculates cell-to-cell correlations based on
variance-stabilized transformed differences between observed and
perturbed gene expression matrices. The function computes correlations
either for all cell pairs or, when a KNN graph is provided, only for
cells connected in the graph.

The primary steps of this analysis are:

1.  Variance-stabilizing transformation: For each cell, observed and
    perturbed expression values are transformed to reduce variance
    effects.

2.  Cell-to-cell correlation computation: For each cell, the function
    calculates correlation values with either all other cells or only
    the nearest neighbors, as defined by the KNN graph.

3.  Parallel processing: Correlations are computed in parallel for
    optimized performance.
