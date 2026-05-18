# HeatmapDistance: Generate Heatmaps for Original and Perturbed Matrices

This function generates two heatmaps from two matrices (original and
perturbed) and displays them side by side on the same color scale.

## Usage

``` r
HeatmapDistance(
  df_original,
  df_perturbed,
  custom_palette = NULL,
  title_original = "Original Assay Cluster Similarity Distance",
  title_perturbed = "Perturbed Assay Cluster Similarity Distance",
  custom_order = NULL,
  col_fun = NULL
)
```

## Arguments

- df_original:

  A numeric matrix representing the original (unperturbed) data.

- df_perturbed:

  A numeric matrix representing the perturbed data.

- custom_palette:

  A vector of colors to define the color palette. Defaults to a red/blue
  gradient.

- title_original:

  The title for the original heatmap. Defaults to "Original Assay
  Cluster Similarity Distance".

- title_perturbed:

  The title for the perturbed heatmap. Defaults to "Perturbed Assay
  Cluster Similarity Distance".

- custom_order:

  A character vector specifying the order of clusters for the x and y
  axes.

- col_fun:

  Optional color mapping function created by `circlize::colorRamp2()`.
  If provided, `HeatmapDistance()` will derive a shared color scale from
  the `breaks` and `colors` attributes of `col_fun` and apply it to both
  heatmaps. When `col_fun` is supplied, `custom_palette` is ignored.

## Value

A patchwork object combining the two heatmaps.

## Examples

``` r
p <- HeatmapDistance(df_edist_observed, df_edist_perturbed) # , custom_order = custom_order
#> Error: object 'df_edist_observed' not found
```
