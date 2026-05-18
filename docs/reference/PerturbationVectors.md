# Calculate Perturbation Transition Vectors for Cell Embeddings

This function calculates transition vectors for cells in a Seurat object
based on a specified perturbation, providing insight into cell-state
shifts on a dimensionality-reduced embedding.

## Usage

``` r
PerturbationVectors(
  seurat_obj,
  perturbation_name,
  reduction = "umap",
  n_threads = 4,
  arrow_scale = 1,
  max_pct = 0.9,
  use_velocyto = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data, including embeddings and
  graphs.

- perturbation_name:

  A character string specifying the name of the perturbation. This is
  used to retrieve the perturbation-specific transition probability
  graph.

- reduction:

  A character string specifying the name of the dimensional reduction to
  use for the embedding (default: 'umap').

- n_threads:

  An integer specifying the number of threads to use for parallel
  processing (default: 4).

- arrow_scale:

  A numeric value to scale the transition vectors, adjusting their
  length on the plot (default: 1).

- max_pct:

  A numeric value between 0 and 1 indicating the maximum percentile of
  vector lengths to display. This parameter caps the extreme values,
  limiting outlier vector lengths for better visualization (default:
  0.90).

## Value

A list containing two data frames:

- ars:

  A data frame with coordinates for plotting transition vectors,
  including initial points (`x0`, `y0`) and final points (`x1`, `y1`)
  for each cell.

- arsd:

  A data frame with the transition vector distances (`xd`, `yd`) for
  each cell, normalized based on `max_pct`.

## Details

This function retrieves the 2D embedding coordinates for the specified
reduction in the Seurat object. It then extracts the
perturbation-specific transition probability graph from the Seurat
object and uses it to calculate transition vectors with the
`embArrows_velocyto` helper function. Vector lengths are scaled by
`arrow_scale` and capped at the specified `max_pct` to handle extreme
values. Cells with `NA` values in either `ars` or `arsd` are excluded
from the final output.
