# CustomPerturbation

This function enables in-silico gene expression perturbation analysis
for user-selected genes using a co-expression network. It applies
primary perturbations to the selected genes, propagates the signal
throughout their co-expression neighborhood, and computes cell-cell
transition probabilities.

## Usage

``` r
CustomPerturbation(
  seurat_obj,
  selected_features,
  perturb_dir,
  perturbation_name,
  graph,
  group.by = NULL,
  group_name = NULL,
  n_connections = NULL,
  perturb_mode = "zinb",
  random_connections = FALSE,
  exclude_grey_genes = FALSE,
  n_iters = 3,
  delta_scale = 0.2,
  row_normalize = FALSE,
  prune_network = FALSE,
  prune_percentile = 0.95,
  corr_sigma = 0.05,
  n_threads = 4,
  use_velocyto = TRUE,
  use_graph_tp = FALSE,
  layer = "counts",
  slot = "counts",
  assay = "RNA",
  n_workers = 1,
  custom_network = NULL,
  custom_modules = NULL,
  custom_weights = NULL,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object containing the gene expression and co-expression data.

- selected_features:

  A character vector of gene names to perturb (e.g. candidate driver
  genes).

- perturb_dir:

  A numeric value determining the type of perturbation to apply:

  - negative for knock-down,

  - positive for knock-in,

  - 0 for knock-out.

- perturbation_name:

  A string representing the name of the in-silico perturbation. This
  will be stored as a new assay in the Seurat object.

- graph:

  Name of the cell-cell graph in `Graphs(seurat_obj)`, used for
  transition probability calculations.

- group.by:

  Optional. A string specifying the column in `seurat_obj@meta.data`
  used for cell grouping.

- group_name:

  Optional. A string or vector specifying the group(s) within `group.by`
  to use for perturbation. If NULL, perturbation is applied to all
  cells.

- n_connections:

  Number of co-expressed neighborhood genes to include alongside the
  selected features. If NULL, defaults to the median module size.
  Ignored when `custom_weights` is provided.

- perturb_mode:

  Character. Perturbation model passed to `ApplyPerturbation`. One of
  `"zinb"` (default) or `"multiplicative"`. See
  [`ApplyPerturbation`](https://smorabit.github.io/compact/reference/ApplyPerturbation.md)
  for details.

- random_connections:

  Logical. If TRUE, selects random non-selected genes instead of the
  most strongly co-expressed genes.

- exclude_grey_genes:

  Logical. If TRUE, excludes grey (unassigned) genes from the
  co-expression network.

- n_iters:

  Number of iterations for propagating the perturbation signal through
  the network. Default is 3.

- delta_scale:

  A numeric scaling factor controlling the influence of the propagated
  perturbation. Default is 0.2.

- row_normalize:

  Logical. If TRUE, normalizes the network rows before propagation.
  Default is FALSE.

- prune_network:

  Logical. If TRUE, removes weak edges from the network before
  propagation. Default is FALSE.

- prune_percentile:

  Numeric. Percentile threshold for pruning network edges. Default is
  0.95.

- corr_sigma:

  A numeric scaling factor for adjusting the correlation matrix during
  transition probability calculations. Default is 0.05.

- n_threads:

  Number of threads to use for parallel computation during correlation
  calculations. Default is 4.

- use_velocyto:

  Logical. If TRUE, leverages velocyto.R functions for transition
  probabilities. Default is TRUE.

- use_graph_tp:

  Logical. If TRUE, transition probabilities are computed using the
  cell-cell graph specified in `graph`. Default is FALSE.

- layer:

  Layer in the assay used for the perturbation. Default is 'counts'.

- slot:

  Slot to extract data for aggregation (e.g., "counts", "data", or
  "scale.data"). Default is 'counts'.

- assay:

  Name of the assay in `seurat_obj` containing the expression data.
  Default is 'RNA'.

- custom_network:

  Optional. A pre-computed gene-by-gene network matrix (e.g., a
  correlation matrix or TOM). Must be supplied together with
  `custom_modules`. When provided, bypasses `GetTOM`.

- custom_modules:

  Optional. A data.frame with at least `gene_name` and `module` columns,
  analogous to the output of `GetModules`. Must be supplied together
  with `custom_network`.

- custom_weights:

  Optional. Name of a column in `custom_modules` to use for ranking
  neighborhood genes instead of network connectivity (column sums). The
  top `n_connections` genes by this score are selected as the
  propagation neighborhood.

- wgcna_name:

  Optional. Name of the hdWGCNA experiment in `seurat_obj@misc`. If
  NULL, defaults to the active WGCNA experiment.

## Value

A Seurat object containing the in-silico perturbation results as a new
assay.

## Details

`CustomPerturbation` allows in-silico perturbation of any user-specified
set of genes without requiring them to be module hub genes. The analysis
consists of three steps:

1.  **Primary perturbation**: applies a knock-in, knock-down, or
    knock-out to `selected_features`.

2.  **Signal propagation**: the perturbation signal is propagated to the
    `n_connections` most strongly co-expressed neighborhood genes via
    iterative network diffusion.

3.  **Transition probabilities**: cell-cell transition probabilities are
    computed from the simulated expression state to quantify phenotypic
    shifts.

Either the hdWGCNA TOM (via `wgcna_name`) or a user-supplied
`custom_network` / `custom_modules` pair can be used as the underlying
co-expression network.
