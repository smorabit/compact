# Reload SHAP Output from Disk

This function reloads previously saved SHAP analysis results from a
specified directory (e.g., created by
[`FindShapKeyDriver()`](https://smorabit.github.io/compact/reference/FindShapKeyDriver.md))
and attaches them to a Seurat object.

## Usage

``` r
ReloadShapOutput(seurat_obj, shap_dir, enforce_pins = FALSE)
```

## Arguments

- seurat_obj:

  A `Seurat` object to attach SHAP results to.

- shap_dir:

  Path to the directory containing SHAP outputs (e.g.,
  "results/shap_output_split").

## Value

A `Seurat` object with SHAP results loaded into `seurat_obj@misc$shap`.

## Details

It is useful when SHAP results were not stored in `@misc` or when
working across sessions.

## Examples

``` r
seurat_obj <- ReloadShapOutput(seurat_obj, shap_dir = "results/shap_output_split")
#> Error: object 'seurat_obj' not found
```
