# Build a Vector Field Interpolation Function from a Grid Data Frame

Converts the output of `grid_vectors` into a callable function
`f(c(x, y)) -> c(fx, fy)` suitable for use as the `fun` argument of
[`ggvfields::geom_stream_field()`](https://rdrr.io/pkg/ggvfields/man/geom_stream_field.html).
Interpolation is performed using inverse-distance weighting (IDW) over
the four nearest grid points.

## Usage

``` r
build_field_fun(grid_df, emb = NULL, density_thresh = 0)
```

## Arguments

- grid_df:

  A data frame produced by `grid_vectors`, containing columns
  `start.emb_1`, `start.emb_2` (grid point positions) and `end.xd`,
  `end.yd` (grid point endpoints).

- emb:

  A cells-by-2 matrix of embedding coordinates. Required when
  `density_thresh > 0` to compute the KDE.

- density_thresh:

  Numeric in `[0, 1]`. Normalised KDE density below which the field
  returns zero velocity. Default `0` disables the mask.

## Value

A function that accepts a numeric vector `c(x, y)` and returns the
IDW-interpolated vector components `c(fx, fy)`, or `c(0, 0)` when the
query point falls in a low-density region.

## Details

When `density_thresh > 0`, a 2D kernel density estimate (KDE) of the
cell positions is computed via
[`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html) and
normalised to `[0, 1]`. Any query point whose normalised density falls
below `density_thresh` returns a zero vector `c(0, 0)`, causing
`geom_stream_field()` to terminate the streamline at that point.
