# Plot SHAP Beeswarm Summary for Top Driver Genes

Generates a beeswarm-style summary plot of SHAP values for the top `n`
driver genes, as computed from
[`FindShapKeyDriver()`](https://smorabit.github.io/compact/reference/FindShapKeyDriver.md).
This visualization helps interpret gene-level contributions to model
predictions, showing both SHAP value and feature value per cell.

## Usage

``` r
BeeswarmplotShap(
  seurat_obj,
  top_n = 20,
  save_plot = FALSE,
  out_dir = NULL,
  filename = "beeswarm_top_shap_DriverGenes.pdf",
  min_color = "#FFCC33",
  max_color = "#6600CC",
  width = 8,
  height = 6,
  dpi = 300
)
```

## Arguments

- seurat_obj:

  A `Seurat` object containing SHAP results from
  [`FindShapKeyDriver()`](https://smorabit.github.io/compact/reference/FindShapKeyDriver.md),
  stored in `seurat_obj@misc$shap`.

- top_n:

  Number of top driver genes to display (default: 20).

- save_plot:

  Save the Beeswarmplot automatically (default: save_plot = FALSE).

- out_dir:

  Output directory for the saved plot file. Will be created if it
  doesn't exist.

- filename:

  Name of the output plot file (e.g.,
  `"beeswarm_top_shap_DriverGenes.pdf"`). File format is inferred from
  the extension (.pdf or .png).

- min_color:

  Color for low feature values (default: `"#FFCC33"`).

- max_color:

  Color for high feature values (default: `"#6600CC"`).

- width:

  Width of the plot in inches (default: 8).

- height:

  Height of the plot in inches (default: 6).

- dpi:

  Resolution in dots per inch (only used for raster outputs like PNG;
  default: 300).

## Value

A `ggplot` object representing the beeswarm summary plot.

## Details

The function extracts SHAP outputs from the `@misc$shap` slot of a
Seurat object and plots per-cell SHAP values for the most important
genes. The plot is saved as a PDF (with rasterized points) or PNG
depending on the filename extension.

## Examples

``` r
# Assuming SHAP analysis was performed using FindShapKeyDriver()
BeeswarmplotShap(seurat_obj, out_dir = "figures/", top_n = 30)
#> Error: object 'seurat_obj' not found

# To save as a raster-aware PDF
BeeswarmplotShap(seurat_obj,save_plot = TRUE,filename = "beeswarm_top_shap_DriverGenes.pdf",out_dir = out_dir)
#> Error: object 'seurat_obj' not found
```
