# StampPerturbProvenance: record baseline provenance for a perturb assay

This function stores lightweight provenance metadata on a perturbation
assay (e.g. created by
[`ModulePerturbation()`](https://smorabit.github.io/compact/reference/ModulePerturbation.md))
indicating which baseline assay and layer(s) it was derived from.

## Usage

``` r
StampPerturbProvenance(
  seurat_obj,
  perturb_assay,
  baseline_assay = "RNA",
  method = "ModulePerturbation",
  overwrite = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- perturb_assay:

  Character; name of the perturbation assay to stamp.

- baseline_assay:

  Character; baseline assay name used to generate the perturbation.

- method:

  Character; description of how the assay was produced.

- overwrite:

  Logical; overwrite existing provenance fields.

## Value

Seurat object with provenance written to the assay misc.

## Details

This enables
[`CombinePerturbAssays()`](https://smorabit.github.io/compact/reference/CombinePerturbAssays.md)
to perform a strict, deterministic check that all perturbation assays
being combined were generated from the same baseline representation
(e.g. baseline assay `RNA` vs `SCT`, and baseline layer `counts`).

Provenance is stored in `seurat_obj[[perturb_assay]]@misc$COMPACT` with
fields:

- `baseline_assay`

- `baseline_layers` (e.g. `"counts"`)

- `method`

- `stamped_at`
