# ModulePerturbation

This function enables in-silico gene expression perturbation analysis
using a co-expression network. It applies primary perturbations to hub
genes, propagates the signal throughout the co-expression network, and
computes cell-cell transition probabilities.

## Usage

``` r
ModulePerturbation(
  seurat_obj,
  mod,
  perturb_dir,
  perturbation_name,
  graph = "RNA_nn",
  group.by = NULL,
  group_name = NULL,
  n_hubs = 5,
  perturb_mode = "zinb",
  n_iters = 3,
  expand_module = 0,
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

- mod:

  Name of the co-expression module to perturb.

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

- n_hubs:

  Number of hub genes to perturb in the selected co-expression module.
  Default is 5.

- perturb_mode:

  Character. Perturbation model passed to `ApplyPerturbation`. One of
  `"zinb"` (default, ZINB-based additive model) or `"multiplicative"`
  (cell-specific fold-change scaling). See
  [`ApplyPerturbation`](https://smorabit.github.io/compact/reference/ApplyPerturbation.md)
  for full details.

- n_iters:

  Number of iterations for propagating the perturbation signal through
  the network. Default is 3.

- delta_scale:

  A numeric scaling factor controlling the influence of the propagated
  perturbation. Default is 0.2.

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

- wgcna_name:

  Optional. Name of the hdWGCNA experiment in `seurat_obj@misc`. If
  NULL, defaults to the active WGCNA experiment.

## Value

A Seurat object containing the in-silico perturbation results as a new
assay.

## Details

Following co-expression network analysis with hdWGCNA,
ModulePerturbation performs in-silico gene expression perturbation
analysis in three steps.

1.  **Primary perturbation**: A direct perturbation is applied to the
    hub genes of the selected module via `ApplyPerturbation`. In ZINB
    mode, expression is modeled by a Zero-Inflated Negative Binomial
    distribution; in multiplicative mode, each cell's hub gene counts
    are scaled by `perturb_dir` as a fold change.

2.  **Log-space signal propagation**: The perturbation delta is
    propagated through the co-expression network in log-normalized
    expression space via `ApplyPropagation`. Working in log space avoids
    the count-space floor asymmetry (where down-regulation deltas
    collapse to zero for non-expressing cells), ensuring that both up-
    and down-regulation produce correctly signed, non-zero deltas in all
    downstream genes.

3.  **Transition probability computation**: Cell-cell transition
    probabilities are computed from the log-space simulated expression
    via `PerturbationTransitions`.
