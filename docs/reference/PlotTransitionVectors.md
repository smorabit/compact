# Plot Transition Vectors on a Reduced Dimensional Embedding

This function visualizes cell-state transitions in a Seurat object by
plotting transition vectors on a dimensionality-reduced embedding (e.g.,
UMAP). It computes transition vectors based on a specified perturbation
and overlays them on a scatter plot of the cells, allowing insights into
cell-state shifts following the perturbation. Four visualization modes
are available: classic grid arrows (`"arrows"`), streamlines
(`"streamlines"`), both overlaid (`"both"`), or one arrow per cell
(`"cells"`).

## Usage

``` r
PlotTransitionVectors(
  seurat_obj,
  perturbation_name,
  color.by,
  reduction = "umap",
  n_threads = 4,
  arrow_scale = 1,
  grid_resolution = 25,
  max_pct = 0.9,
  point_alpha = 0.25,
  point_size = 1,
  arrow_size = 0.25,
  raster_dpi = 300,
  arrow_alpha = TRUE,
  use_velocyto = FALSE,
  plot_mode = c("arrows", "streamlines", "both", "cells"),
  stream_n = 15,
  stream_L = 2,
  stream_normalize = TRUE,
  stream_linewidth = 0.5,
  stream_alpha = 0.8,
  stream_color = "black",
  stream_density_thresh = 0
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data, including metadata and
  embeddings.

- perturbation_name:

  A character string specifying the name of the perturbation. This is
  used to calculate transition vectors based on the perturbed assay.

- color.by:

  A character string indicating the metadata column to color the scatter
  plot points by.

- reduction:

  A character string specifying the name of the dimensional reduction to
  use for the embedding (default: `'umap'`).

- n_threads:

  An integer specifying the number of threads to use for parallel
  processing (default: 4).

- arrow_scale:

  A numeric value to scale the transition vectors, adjusting their
  length on the plot (default: 1).

- grid_resolution:

  An integer specifying the resolution of the grid for overlaying
  transition vectors. A higher resolution yields a finer grid, which can
  provide more detailed vector placement (default: 25).

- max_pct:

  A numeric value between 0 and 1 indicating the maximum percentile of
  vector lengths to display. This parameter helps limit the display of
  outlier vector lengths, making the plot easier to interpret (default:
  0.90).

- point_alpha:

  A numeric value controlling the transparency of scatter plot points,
  with 1 being fully opaque and 0 being fully transparent (default:
  0.25).

- point_size:

  A numeric value defining the size of the scatter plot points (default:
  1).

- arrow_size:

  Numeric. Line width of the arrow segments (default: 0.25).

- raster_dpi:

  Integer. DPI for rasterising the scatter plot layer (default: 300).

- arrow_alpha:

  Logical. If `TRUE`, arrow transparency is scaled by vector length so
  shorter arrows are more transparent (default: `TRUE`).

- use_velocyto:

  Logical. Use the velocyto.R C backend for embedding-arrow calculation
  instead of the built-in `SparseEmbArrows` (default: `FALSE`).

- plot_mode:

  Character. One of `"arrows"` (default), `"streamlines"`, `"both"`, or
  `"cells"`. `"arrows"` reproduces the classic grid-arrow plot.
  `"streamlines"` replaces arrows with continuous streamlines traced
  through the vector field. `"both"` overlays streamlines on top of the
  arrow plot. `"cells"` draws one arrow per cell at its exact embedding
  position, bypassing grid aggregation entirely — the arrows themselves
  are coloured by `color.by` and no background scatter plot is drawn.
  Recommended for small datasets (\< ~2,000 cells) where grid bins would
  be too sparse. Streamline rendering requires the ggvfields package.

- stream_n:

  Integer. Number of streamline seed points along each axis of the
  embedding grid passed to
  [`ggvfields::geom_stream_field()`](https://rdrr.io/pkg/ggvfields/man/geom_stream_field.html)
  (default: 15).

- stream_L:

  Numeric. Target length of each streamline, in embedding coordinate
  units, passed to
  [`ggvfields::geom_stream_field()`](https://rdrr.io/pkg/ggvfields/man/geom_stream_field.html)
  (default: 2).

- stream_normalize:

  Logical. If `TRUE` (default), all streamlines are drawn at equal
  length. If `FALSE`, length reflects local vector magnitude.

- stream_linewidth:

  Numeric. Line width of the streamlines (default: 0.5).

- stream_alpha:

  Numeric. Transparency of the streamlines, in `[0, 1]` (default: 0.8).

- stream_color:

  Character. Colour of the streamlines (default: `"black"`).

## Value

A ggplot object showing the scatter plot of cells colored by the
specified metadata and overlaid with transition vectors and/or
streamlines indicating cell-state shifts.

## Details

This function retrieves the embedding coordinates from the specified
dimensional reduction in the Seurat object. It then calculates the
transition vectors using `PerturbationVectors`, which returns vector
coordinates (`ars`) and distances (`arsd`). A grid of arrows is created
using `grid_vectors`, allowing for visual simplification of transitions
by grouping vectors into a grid at the specified resolution.

When `plot_mode` includes streamlines, the grid-aggregated vector field
is converted into an inverse-distance-weighted interpolation function
via the internal helper
[`build_field_fun()`](https://smorabit.github.io/compact/reference/build_field_fun.md),
which is passed directly to
[`ggvfields::geom_stream_field()`](https://rdrr.io/pkg/ggvfields/man/geom_stream_field.html).
The ggvfields package must be installed for streamline rendering.
