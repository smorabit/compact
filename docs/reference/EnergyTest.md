# EnergyTest

Performs E-testing on a Seurat object. Computes E-test statistics for
each group in a Seurat object, using the E-distance in space given by
reduction to the group defined by control.

## Usage

``` r
EnergyTest(
  seurat_object,
  groupby,
  control,
  reduction = "pca",
  verbose = TRUE,
  permutations = 1000
)
```

## Arguments

- seurat_object:

  An object of class Seurat.

- groupby:

  An object of class character. Points to the column in the Seurat
  object's meta data that contains the group labels.

- control:

  An object of class character. The group that is used as the control
  for the statistical comparison.

- reduction:

  An object of class character. The reduction / embedding in
  seurat_object that is used to compute the E-distance in. Can be "pca",
  "harmony", or any linear dimensionality reduction.

- verbose:

  An object of class logical. If TRUE, prints messages. Default is TRUE.

- permutations:

  An object of class integer. The number of permutations used to compute
  the p-value. Default is 1000.

## Value

Returns an object of class data.frame. For each group contains the
E-test p-value and the E-distance to control group.

## Examples

``` r
    # Add some code illustrating how to use the function
```
