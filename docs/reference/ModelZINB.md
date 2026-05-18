# ModelZINB

This function models the expression of a selected feature as a
zero-inflated negative binomial (ZINB) distribution.

## Usage

``` r
ModelZINB(
  seurat_obj,
  feature,
  cells_use = NULL,
  layer = "counts",
  slot = "counts",
  add_zero = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object

- feature:

  List of cell barcodes from the Seurat object

- slot:

  Slot to extract data for aggregation. Default = 'counts'

## Value

The zero-inflated negative binomial model fit to the observed data from
the selected feature.

## Details

ModelZINB is an internal helper function that calls on zeroinfl from the
pscl package to model expression of a selected feature.
