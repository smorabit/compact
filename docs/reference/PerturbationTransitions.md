# PerturbationTransitions

This function computes cell-cell transition probabilities based on an
in-silico perturbation experiment.

## Usage

``` r
PerturbationTransitions(
  seurat_obj,
  perturbation_name,
  features,
  graph,
  use_velocyto = FALSE,
  use_graph_tp = FALSE,
  corr_sigma = 0.05,
  n_threads = 4,
  layer = "data",
  slot = "data",
  assay = "RNA"
)
```

## Arguments

- seurat_obj:

  A Seurat object

- perturbation_name:

  A name for the in-silico perturbation that will be stored in the
  Seurat obejct

- features:

  Selected features to use for the transition probability calculation

- graph:

  Name of the cell-cell graph in the Graphs(seurat_obj)

- corr_sigma:

  A numeric scaling factor for the correlation matrix.

- n_threads:

  Number of threads for the correlation calculation

- slot:

  Slot to extract data for aggregation. Default = 'data'

- assay:

  Assay in seurat_obj containing expression information.

## Value

A Seurat object with the cell-cell transition probabilities stored in
the Graphs slot of the Seurat object.
