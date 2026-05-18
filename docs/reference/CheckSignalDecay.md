# Check for Propagation Signal Decay

After running `ApplyPropagation`, checks whether the propagated signal
reaching downstream (non-directly-perturbed) genes has decayed below a
detectable threshold. With `row_normalize = TRUE` and
`delta_scale <= 1`, the row-normalized network acts as a stochastic
matrix — repeated multiplication can only redistribute signal, never
amplify it — so the propagated delta decays geometrically toward zero
and is often lost to integer rounding.

## Usage

``` r
CheckSignalDecay(
  delta_initial,
  delta_final,
  row_normalize,
  delta_scale,
  n_iters,
  min_mean_delta = 0.001
)
```

## Arguments

- delta_initial:

  A features-by-cells matrix of the initial perturbation delta
  (`exp_per - exp`) before any propagation iterations.

- delta_final:

  A features-by-cells matrix of the propagated delta after `n_iters`
  iterations (i.e., the delta added to `exp_per` to produce the final
  propagated expression).

- row_normalize:

  Logical. Whether the network was row-normalized during propagation.

- delta_scale:

  Numeric. The `delta_scale` used during propagation; reported in the
  warning.

- n_iters:

  Integer. The `n_iters` used during propagation; reported in the
  warning.

- min_mean_delta:

  Numeric. Minimum mean absolute delta in downstream genes considered a
  detectable signal. Defaults to `0.001`, which is appropriate for
  log-normalized expression space.

## Value

Invisibly returns a named list with:

- `mean_delta_initial`:

  Mean absolute delta across directly perturbed gene-cell entries.

- `mean_delta_downstream`:

  Mean absolute delta across downstream gene-cell entries after
  propagation.

## Details

Directly perturbed genes are identified as rows where the sum of
`|delta_initial|` is non-zero. All remaining rows are treated as
downstream. A warning is issued when the mean absolute downstream delta
is below `min_mean_delta`, with specific guidance when
`row_normalize = TRUE` and `delta_scale <= 1` — the parameter
combination that guarantees decay regardless of `n_iters`.
