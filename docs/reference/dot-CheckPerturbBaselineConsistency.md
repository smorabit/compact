# .CheckPerturbBaselineConsistency: heuristic sanity-check against a baseline

This function performs a heuristic check that a perturbation expression
matrix is numerically consistent with a specified baseline matrix. It is
intended for cases where explicit provenance is missing (i.e. no call to
[`StampPerturbProvenance()`](https://smorabit.github.io/compact/reference/StampPerturbProvenance.md))
and the user wants a guardrail against combining assays generated from
mismatched baselines/normalizations/slots.

## Usage

``` r
.CheckPerturbBaselineConsistency(
  base_mat,
  pert_mat,
  assay_name = "pert",
  eps = 0.001,
  max_frac_changed = 0.35,
  max_abs_median = 0.05,
  max_abs_mean = 0.05,
  strict = FALSE,
  n_genes = 2000,
  n_cells = 500,
  seed = NULL
)
```

## Arguments

- base_mat:

  Baseline expression matrix (genes x cells).

- pert_mat:

  Perturbation expression matrix (genes x cells).

- assay_name:

  Character label used in messages.

- eps:

  Numeric; values with `|Δ| <= eps` are treated as unchanged.

- max_frac_changed:

  Numeric in (0,1); warn/stop if too many entries differ.

- max_abs_median:

  Numeric; warn/stop if `median(|Δ|)` is too large.

- max_abs_mean:

  Numeric; warn/stop if `|mean(Δ)|` is too large.

- strict:

  Logical; if `TRUE`, stop on suspicious results; otherwise warn.

- n_genes:

  Integer; number of genes to sample for the check (default 2000).

- n_cells:

  Integer; number of cells to sample for the check (default 500).

- seed:

  Optional integer; set for reproducible sampling.

## Value

Invisibly returns a list of summary statistics and a suspicious flag.

## Details

The check computes \\\Delta = P - B\\ on overlapping genes and
summarizes:

- fraction of entries with `|Δ| > eps`

- `median(|Δ|)`

- `|mean(Δ)|` (detects global shifts)

If any summary exceeds user-defined thresholds, the function warns or
stops.

Note: This is not a proof of shared provenance; it is a numeric sanity
check.
