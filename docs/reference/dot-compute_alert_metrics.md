# Compute concordance and robustness metrics for AlertSystem

For each group × pathway (logFC column) this computes:

- effect : mean or median logFC (matching `fun`)

- concordance : fraction of cells whose logFC sign matches the group
  effect sign

- robust_gt_thresh : fraction of cells whose logFC is both: (i) above
  the magnitude threshold (\|logFC\| ≥ robust_effect_thresh) and (ii)
  has the same sign as the group-level effect.

## Usage

``` r
.compute_alert_metrics(
  seurat_obj,
  score_pattern = "_logFC$",
  group.by = NULL,
  fun = c("mean", "median"),
  compute_concordance = TRUE,
  compute_robustness = TRUE,
  robust_effect_thresh = 0,
  base_tag = NULL,
  pert_tag = NULL
)
```

## Details

This replaces the previous bootstrap-based credibility metric, which
tended to saturate at 0/1 for large cell numbers. Robustness is now
defined directly from the per-cell logFC distribution.
