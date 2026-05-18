# Check for Signal Saturation After Network Propagation

After running `ApplyPropagation`, checks whether the propagated
expression has saturated at the biological floor (0) or ceiling (max
observed expression). Heavy saturation can cause up- and down-regulation
perturbations to produce indistinguishable transition vector fields.

## Usage

``` r
CheckSaturation(
  exp_obs,
  exp_prop,
  delta_scale,
  n_iters,
  apply_ceiling = FALSE,
  max_obs = NULL,
  saturation_thresh = 0.5
)
```

## Arguments

- exp_obs:

  A features-by-cells matrix of baseline (unperturbed) expression
  counts.

- exp_prop:

  A features-by-cells matrix of the final propagated expression counts
  (output of `ApplyPropagation`).

- delta_scale:

  Numeric. The `delta_scale` used during propagation; reported in the
  warning.

- n_iters:

  Integer. The `n_iters` used during propagation; reported in the
  warning.

- apply_ceiling:

  Logical. Whether a ceiling constraint was applied during propagation.
  Default `FALSE`.

- max_obs:

  Numeric vector. Per-gene ceiling values (length equal to
  `nrow(exp_prop)`). Only used when `apply_ceiling = TRUE`.

- saturation_thresh:

  Numeric. Fraction of perturbed gene-cell entries that must be
  saturated to trigger a warning. Default `0.5`.

## Value

Invisibly returns a named list with `floor_frac` and, if
`apply_ceiling = TRUE`, `ceil_frac` — the computed saturation fractions.

## Details

Perturbed genes are identified as rows where the sum of
`|exp_prop - exp_obs|` is non-zero. Among those rows, floor saturation
is the fraction of gene-cell entries where propagated expression equals
0. If `apply_ceiling = TRUE`, ceiling saturation is the fraction of
entries at or above `max_obs`. A warning is issued for each fraction
exceeding `saturation_thresh`.
