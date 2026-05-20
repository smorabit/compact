# MacrostateTransitions

Coarse-grains the cell-level perturbation transition probability (TP)
matrix into a group-level summary matrix (K × K), analogous to the
macrostate transition probability heatmap in CellRank (Lange et al.
2022, Nature Methods, Fig. 2c). For each pair of user-defined cell
groups, computes the average transition probability flowing from cells
in the source group to cells in the destination group under the
specified perturbation.

The diagonal of the resulting matrix gives a *stability index* per
group: the mean probability that a cell in that group transitions to
another cell *within the same group*. High stability (close to 1)
indicates attractor-like behavior — the perturbation does not move cells
out of this group. Low stability indicates a transient source state from
which cells are being redirected toward other groups.

Unlike `PlotTransitionVectors`, this analysis is entirely independent of
the low-dimensional embedding. All computations operate on the
full-dimensional transition matrix, avoiding the distortions introduced
by UMAP/t-SNE projections.

## Usage

``` r
MacrostateTransitions(
  seurat_obj,
  perturbation_name,
  graph,
  group.by,
  normalize_rows = TRUE,
  min_group_size = 5,
  store_result = TRUE,
  result_name = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object with a computed perturbation transition probability
  matrix stored in `seurat_obj@graphs` (produced by
  `ModulePerturbation`, `TFPerturbation`, or `PerturbationTransitions`).

- perturbation_name:

  Character. The base name of the perturbation. The function looks up
  the TP graph as `paste0(perturbation_name, '_tp')`.

- graph:

  Character. Name of the KNN graph stored in `seurat_obj@graphs` (e.g.,
  `"RNA_nn"`). Used to mask transition probabilities to the biological
  manifold, identical to the masking applied in `PredictAttractors` and
  `PredictPerturbationTime`.

- group.by:

  Character. Metadata column containing group labels used for
  coarse-graining (e.g., `"functional_anno"`, `"seurat_clusters"`).

- normalize_rows:

  Logical. If `TRUE` (default), each row of the coarse-grained matrix Q
  is normalized to sum to 1, making it row-stochastic. If `FALSE`, rows
  sum to the fraction of total outflow captured within the non-dropped
  groups (useful for diagnosing how much flow escapes excluded small
  groups).

- min_group_size:

  Integer. Groups with fewer than this many cells are excluded before
  computing Q, with a warning. Default: `5`.

- store_result:

  Logical. If `TRUE` (default), the result list is stored in
  `seurat_obj@misc$MacrostateTransitions[[result_name]]`.

- result_name:

  Character. Key used for `@misc` storage. Defaults to
  `perturbation_name`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

## Value

The Seurat object (if `store_result = TRUE`) with the following list
stored in `seurat_obj@misc$MacrostateTransitions[[result_name]]`:

- `Q`:

  K × K numeric matrix. Coarse-grained (row-stochastic) transition
  probability matrix. `Q[k1, k2]` is the average probability of a cell
  in group `k1` transitioning to a cell in group `k2`.

- `stability`:

  Named numeric vector of length K. Diagonal of Q; the self-transition
  probability (stability index) for each group.

- `group_sizes`:

  Named integer vector. Number of cells per group after small-group
  filtering.

- `group.by`:

  Character. The metadata column used.

- `perturbation_name`:

  Character. The perturbation name.

If `store_result = FALSE`, returns the result list directly.

## Details

**Mathematical summary:**

Let P be the n × n row-stochastic cell-level transition matrix (after
transposing the stored column-stochastic TP, applying KNN masking, and
row-normalizing). Let M be the n × K sparse indicator matrix where
`M[i, k] = 1` if cell i belongs to group k. Then:

\$\$Q = \mathrm{diag}(1/n_k) \cdot M^\top \cdot P \cdot M\$\$

where \\n_k\\ is the number of cells in group k. Q is row-stochastic
(each row sums to 1) because P is row-stochastic and M sums each row
over all K groups.

**Complexity:** O(n · k_nn · K). The bottleneck is `P %*% M` (sparse n ×
n times sparse n × K). For n = 150,000, k_nn = 20, K = 20 this requires
~3M floating-point operations and ~24 MB of sparse storage. No dense n ×
n matrix is ever materialized.

## References

Lange, M. et al. CellRank for directed single-cell fate mapping. *Nature
Methods* **19**, 159–170 (2022).
[doi:10.1038/s41592-021-01346-6](https://doi.org/10.1038/s41592-021-01346-6)

## See also

[`PlotMacrostateTransitions`](https://smorabit.github.io/compact/reference/PlotMacrostateTransitions.md),
[`PredictAttractors`](https://smorabit.github.io/compact/reference/PredictAttractors.md),
[`PredictFates`](https://smorabit.github.io/compact/reference/PredictFates.md),
[`VectorFieldCoherence`](https://smorabit.github.io/compact/reference/VectorFieldCoherence.md)
