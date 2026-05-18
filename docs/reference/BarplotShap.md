# Plot Bar Chart of Top SHAP Driver Genes

Generates a horizontal bar plot of the top `n` genes ranked by mean
absolute SHAP value, as computed by
[`FindShapKeyDriver()`](https://smorabit.github.io/compact/reference/FindShapKeyDriver.md).
The plot provides an interpretable summary of the most important
features (genes) driving model predictions.

## Usage

``` r
BarplotShap(
  seurat_obj,
  top_n = 20,
  bar_color = "steelblue",
  save_plot = FALSE,
  out_dir = NULL,
  filename = "barplot_top_shap_DriverGenes.pdf",
  width = 8,
  height = 6,
  dpi = 300
)
```

## Arguments

- seurat_obj:

  A `Seurat` object containing SHAP results from
  [`FindShapKeyDriver()`](https://smorabit.github.io/compact/reference/FindShapKeyDriver.md),
  with a `shap_summary` table stored in `seurat_obj@misc$shap`.

- top_n:

  Number of top genes to display (default: 20).

- bar_color:

  Fill color for the bars (default: `"steelblue"`).

- save_plot:

  Save the Barplot automatically (default: save_plot = FALSE).

- out_dir:

  Output directory for the saved plot. Created if it doesn't exist
  (default: NULL).

- filename:

  Name of the plot file to save (e.g.,
  `"barplot_top_shap_DriverGenes.pdf"`). File format is determined by
  the extension (`.pdf` or `.png`).

- width:

  Plot width in inches (default: 8).

- height:

  Plot height in inches (default: 6).

- dpi:

  Resolution in dots per inch for raster formats like PNG (default:
  300).

## Value

A `ggplot` object containing the bar plot of SHAP summary values.

## Details

The function pulls SHAP summary statistics from
`seurat_obj@misc$shap$shap_summary`, creates a bar plot of the top `n`
genes, and saves it to disk as a PNG or PDF.

## Examples

``` r
# After running FindShapKeyDriver():
BarplotShap(seurat_obj, top_n = 30)
#> Error: object 'seurat_obj' not found

# Save as a PDF with a custom color
BarplotShap(seurat_obj,save_plot = TRUE,filename = "barplot_top_shap_DriverGenes.pdf", bar_color = "#F96815")
#> Error: object 'seurat_obj' not found
```
