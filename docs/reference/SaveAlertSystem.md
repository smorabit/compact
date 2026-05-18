# Save the full AlertSystem state to disk

This helper writes the complete `seurat_obj@misc[[misc_name]]` list to
an RDS file, so that AlertSystem results can be re-attached later
without recomputing. This is equivalent to the RDS file written by
`AlertSystemScore(..., out_dir = ...)` (the `*_AlertSystem_info.rds`
file), but can be called manually at any time.

## Usage

``` r
SaveAlertSystem(seurat_obj, file, misc_name = "AlertSystem")
```

## Arguments

- seurat_obj:

  A Seurat object with AlertSystem results stored in
  `seurat_obj@misc[[misc_name]]`.

- file:

  Path to an .rds file where the state will be saved (e.g.
  `"AlertSystem_M6_M6_up_AlertSystem_info.rds"`).

- misc_name:

  Name of the misc slot containing AlertSystem results. Default:
  `"AlertSystem"`.

## Value

Invisibly returns `file`.
