# Calculate Log2 Fold Change of In-Silico Perturbation

Computes the log2 fold change (Log2FC), mean delta, and average baseline
expression for genes following an in-silico perturbation. It compares
the perturbed assay against the baseline observed expression.

## Usage

``` r
PerturbationLog2FC(
  seurat_obj,
  perturbation_name,
  assay_obs = "RNA",
  module = NULL,
  n_hubs = 10,
  pseudocount = 1,
  layer_norm = "data",
  layer_counts = "counts",
  wgcna_name = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object containing both the baseline and perturbed data.

- perturbation_name:

  Character. The name of the assay containing the perturbed expression
  data.

- assay_obs:

  Character. The name of the baseline assay. Default is `'RNA'`.

- module:

  Character. Optional. The specific hdWGCNA module name to analyze. If
  provided, the output will include a column indicating whether each
  gene is a top hub.

- n_hubs:

  Integer. The number of top hub genes to flag if `module` is provided.
  Default is `10`.

- pseudocount:

  Numeric. Pseudocount added to average expression before log2
  transformation to prevent infinite values. Default is `1`.

- layer_norm:

  Character. The layer/layer containing normalized data. Default is
  `'data'`.

- layer_counts:

  Character. The layer/layer containing raw counts. Default is
  `'counts'`.

- wgcna_name:

  Character. The name of the hdWGCNA experiment in `seurat_obj@misc`. If
  NULL, uses the active WGCNA experiment.

## Value

A data.frame containing gene names, mean delta (raw counts difference),
average baseline expression, Log2FC, and optionally module hub status.
