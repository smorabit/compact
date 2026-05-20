# TFPerturbation

This function enables in-silico transcription factor (TF) perturbation
analysis using a regulatory network derived from co-expression and
regulon information.

## Usage

``` r
TFPerturbation(
  seurat_obj,
  selected_tf,
  perturb_dir,
  perturbation_name,
  graph = "RNA_nn",
  group.by = NULL,
  group_name = NULL,
  perturb_mode = "zinb",
  n_iters = 1,
  delta_scale = 1,
  row_normalize = FALSE,
  prune_network = FALSE,
  prune_percentile = 0.95,
  corr_sigma = 0.05,
  n_threads = 4,
  use_velocyto = FALSE,
  use_graph_tp = FALSE,
  depth = 2,
  target_type = "both",
  use_regulons = TRUE,
  layer = "counts",
  slot = "counts",
  assay = "RNA",
  n_workers = 1,
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- selected_tf:

  The name of the transcription factor (TF) to perturb.

- perturb_dir:

  A numeric determining the type of perturbation to apply. Negative
  values for knock-down, positive for knock-in, and 0 for knock-out.

- perturbation_name:

  A name for the in-silico perturbation that will be stored in the
  Seurat object.

- graph:

  Name of the cell-cell graph in the Graphs(seurat_obj). Default =
  "RNA_nn".

- perturb_mode:

  Character. Perturbation model passed to `ApplyPerturbation`. One of
  `"zinb"` (default) or `"multiplicative"`. See
  [`ApplyPerturbation`](https://smorabit.github.io/compact/reference/ApplyPerturbation.md)
  for details.

- n_iters:

  The number of times to apply the signal propagation throughout the TF
  regulatory network. Default = 1.

- delta_scale:

  A numeric factor scaling the perturbation during propagation. Default
  = 1.

- corr_sigma:

  A numeric scaling factor for the correlation matrix. Default = 0.05.

- n_threads:

  Number of threads for the correlation calculation. Default = 4.

- use_velocyto:

  Logical indicating whether to compute velocity-based transition
  probabilities. Default = FALSE.

- use_graph_tp:

  Logical indicating whether to use the graph topology for transition
  probabilities. Default = FALSE.

- depth:

  The depth of the regulatory network to use for target identification.
  Default = 2.

- target_type:

  A string specifying the type of targets to include ("upstream",
  "downstream", or "both"). Default = "both".

- use_regulons:

  Logical indicating whether to use regulons for TF-target
  relationships. Default = TRUE.

- layer:

  Layer of the assay containing expression data. Default = "counts".

- slot:

  Slot of the assay containing expression data. Default = "counts".

- assay:

  Assay in seurat_obj containing expression information. Default =
  "RNA".

- wgcna_name:

  The name of the hdWGCNA experiment in the seurat_obj@misc slot. If
  NULL, uses the active WGCNA experiment.

## Value

A Seurat object containing the in-silico TF perturbation results as a
new assay.

## Details

TFPerturbation enables in-silico transcription factor perturbation by
simulating changes in TF activity and propagating these changes
throughout the associated regulatory network. The workflow consists of
the following steps:

1.  **Primary Perturbation**: A primary in-silico perturbation is
    applied to the selected TF. The observed expression levels of the TF
    are adjusted by simulating changes using a perturbation direction
    (`perturb_dir`) and generating a new expression matrix.

2.  **Signal Propagation**: The perturbation signal is propagated
    throughout the TF regulatory network using the adjacency matrix
    derived from TF-target relationships. The propagation process
    considers the network structure and propagates the signal over
    `n_iters` iterations.

3.  **Transition Probability Computation**: Transition probabilities
    between cells are calculated based on the perturbation, optionally
    incorporating velocity or graph topology information.

This method leverages co-expression and regulatory network information
to study the potential downstream effects of transcription factor
perturbations on gene expression and cellular states.
