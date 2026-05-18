# CombinePerturbAssays: create a composite perturb assay by delta superposition (sparse-safe; COUNTS-only)

This function combines multiple perturbation assays into a single
composite assay using additive superposition of perturbation deltas
relative to a shared baseline on the `counts` layer:
\$\$P\_{\mathrm{combined}} = B + \sum_i (P_i - B)\$\$

## Usage

``` r
CombinePerturbAssays(
  seurat_obj,
  perturb_assays,
  baseline_assay = "RNA",
  new_assay = paste(perturb_assays, collapse = "_"),
  require_provenance = FALSE,
  strict_provenance_check = TRUE,
  run_consistency_check_if_missing = TRUE,
  strict_consistency_check = FALSE,
  eps = 0.001,
  max_frac_changed = 0.35,
  max_abs_median = 0.05,
  max_abs_mean = 0.05,
  allow_unverified = TRUE,
  allow_dense = FALSE,
  normalize_method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing baseline and perturbation assays.

- perturb_assays:

  Character vector of perturbation assay names to combine.

- baseline_assay:

  Character; baseline assay name (default `"RNA"`).

- new_assay:

  Character; name of the new composite assay.

- require_provenance:

  Logical; if TRUE, missing provenance triggers stop/warn.

- strict_provenance_check:

  Logical; if TRUE, stop on provenance mismatch; else warn.

- run_consistency_check_if_missing:

  Logical; if TRUE, run heuristic checks when provenance missing.

- strict_consistency_check:

  Logical; if TRUE, stop on suspicious heuristic results; else warn.

- eps, max_frac_changed, max_abs_median, max_abs_mean:

  Thresholds for heuristic checks.

- allow_unverified:

  Logical; if provenance missing and heuristic checks disabled, proceed
  only if TRUE.

- allow_dense:

  Logical; if FALSE (default), stop when any requested layer returns a
  dense matrix (to avoid RAM blow-ups).

- normalize_method:

  Deprecated/ignored. The composite `data` layer is always produced via
  `log_normalize()` (LogNormalize / `log1p`) to match
  [`ModulePerturbation()`](https://smorabit.github.io/compact/reference/ModulePerturbation.md).
  Kept only for backward compatibility.

- scale.factor:

  Scale factor passed to the internal `log_normalize()` helper (default
  `1e4`). Should match the value used when the input perturb assays were
  generated.

## Value

A Seurat object with a new composite perturbation assay.

## Details

Implementation note: this function avoids creating a large dense
expanded matrix. It combines perturbations strictly on the baseline
feature space (rows = `rownames(baseline)`). Features present in a
perturb assay but absent from baseline are dropped with a warning.

Baseline safety logic:

- If all input assays have provenance (via
  [`StampPerturbProvenance()`](https://smorabit.github.io/compact/reference/StampPerturbProvenance.md)),
  the function verifies recorded `baseline_assay` and `baseline_layers`
  include `"counts"`.

- If provenance is missing for any assay, the function emits a warning
  banner. Optionally, it runs heuristic numeric checks on the `counts`
  layer via
  [`.CheckPerturbBaselineConsistency()`](https://smorabit.github.io/compact/reference/dot-CheckPerturbBaselineConsistency.md)
  for the unstamped assays (not a proof; may be less informative than
  checks on normalized data).

Data layer handling: After combining `counts`, the function generates
the composite assay's `data` layer by calling the package-internal
`log_normalize()` helper (the same one used by
[`ModulePerturbation()`](https://smorabit.github.io/compact/reference/ModulePerturbation.md))
with `col_sums` taken from the BASELINE counts (`colSums(base_counts)`),
not from the composite counts. This preserves the per-cell denominator
across the baseline assay and all perturb / composite assays, so their
`data` layers remain directly comparable.
[`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html)
is intentionally NOT used here because it would recompute library sizes
from the composite counts, breaking that comparability.

Combine by delta-superposition (counts only):

- `comp_counts = base_counts + _i (pert_counts_i - base_counts)`

Requirements:

- baseline assay must have layer `counts`

- each perturb assay must have layer `counts`

- all assays must share identical cell columns

## Examples

``` r
if (FALSE) { # \dontrun{
# After creating perturb assays:
seurat_obj_perturb <- ModulePerturbation(
  seurat_obj_perturb,
  mod = "CBh-M6",
  perturb_dir = 5,
  perturbation_name = "M6_up",
  graph = "RNA_nn",
  n_hubs = 10
)

seurat_obj_perturb <- ModulePerturbation(
  seurat_obj_perturb,
  mod = "CBh-M2",
  perturb_dir = -2,
  perturbation_name = "M2_down",
  graph = "RNA_nn",
  n_hubs = 10
)

# Optional: stamp provenance (recommended)
seurat_obj_perturb <- StampPerturbProvenance(
  seurat_obj_perturb,
  perturb_assay = "M6_up",
  baseline_assay = "RNA"
)

seurat_obj_perturb <- StampPerturbProvenance(
  seurat_obj_perturb,
  perturb_assay = "M2_down",
  baseline_assay = "RNA"
)

seurat_obj_perturb <- CombinePerturbAssays(
  seurat_obj_perturb,
  perturb_assays = c("M6_up", "M2_down"),
  baseline_assay = "RNA",
  new_assay = "M6up_M2down",
  require_provenance = FALSE,
  run_consistency_check_if_missing = FALSE,
  allow_unverified = TRUE
)
} # }
```
