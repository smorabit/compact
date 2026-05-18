# Load a previously saved AlertSystem state into a Seurat object

This helper reads an AlertSystem RDS file (created by
`AlertSystemScore(..., out_dir = ...)` or
[`SaveAlertSystem()`](https://smorabit.github.io/compact/reference/SaveAlertSystem.md))
and attaches it to `seurat_obj@misc[[misc_name]]`.

## Usage

``` r
LoadAlertSystem(
  seurat_obj,
  file,
  misc_name = "AlertSystem",
  overwrite = TRUE,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- file:

  Path to an RDS file containing an AlertSystem list (e.g.
  `"AlertSystem_M6_M6_up_AlertSystem_info.rds"`).

- misc_name:

  Name of the misc slot where the state should be stored. Default:
  `"AlertSystem"`.

- overwrite:

  Logical; if `FALSE` and a state already exists at `misc_name`, an
  error is thrown. Default: `TRUE`.

- verbose:

  Logical; print a short message upon successful load.

## Value

The updated Seurat object.

## Details

Typically, this RDS file will be the `*_AlertSystem_info.rds` file
written by `AlertSystemScore`, which contains `per_cell`, `metrics`, and
`params` in a single list. Any associated metrics CSV (e.g.
`*_metrics.csv`) is not required for loading and plotting, but can be
inspected separately.

Note: this does *not* modify `seurat_obj@meta.data`; per-cell scores
remain accessible via `seurat_obj@misc[[misc_name]]$per_cell`.

## Examples

``` r
if (FALSE) { # \dontrun{
test_obj <- LoadAlertSystem(
  seurat_obj = test_obj,
  file       = "PathToAlertSystemResults/AlertSystem_XXXXXX_AlertSystem_info.rds"
)
} # }
```
