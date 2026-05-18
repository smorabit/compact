# AlertHeatmap: visualize pathway-level perturbation metrics as a heatmap

This function visualizes the AlertSystem metrics stored in
`seurat_obj@misc$AlertSystem$metrics` as a pathways-by-group heatmap. It
selects the top N most positively and/or negatively changed pathways
based on the mean `effect` (summarized logFC) across groups, and can
optionally rescale values for visualization.

## Usage

``` r
AlertHeatmap(
  seurat_obj,
  misc_name = "AlertSystem",
  top_n = 5,
  side = c("both", "positive", "negative"),
  scale_mode = c("none", "unit", "zscale"),
  unit_mode = c("max", "robust"),
  select_by = c("concordance", "effect", "abs_effect"),
  reorder_by = c("effect", "select"),
  unit_q = 0.9,
  reorder = TRUE,
  fill_limits = NULL,
  show_concordance = TRUE,
  show_robustness = FALSE
)
```

## Arguments

- seurat_obj:

  A Seurat object that has been processed with
  [`AlertSystemScore()`](https://smorabit.github.io/compact/reference/AlertSystemScore.md)
  and contains a metrics table in `seurat_obj@misc$AlertSystem$metrics`.
  This can be either from a fresh run of
  [`AlertSystemScore()`](https://smorabit.github.io/compact/reference/AlertSystemScore.md)
  or after re-attaching a previously saved AlertSystem state (e.g. via
  [`LoadAlertSystem()`](https://smorabit.github.io/compact/reference/LoadAlertSystem.md)
  or `readRDS("..._AlertSystem_info.rds")`).

- misc_name:

  Name of the misc slot where AlertSystem results are stored. Default is
  `"AlertSystem"`.

- top_n:

  Integer; number of top pathways to select per side (positive and/or
  negative). For example, `top_n = 5` and `side = "both"` can yield up
  to 10 pathways (5 most positive, 5 most negative).

- side:

  Character string indicating which side(s) of the effect distribution
  to use. One of `"both"`, `"positive"`, or `"negative"`. Default is
  `"both"`.

- scale_mode:

  Character; one of:

  "none"

  :   Use raw effect values (default).

  "unit"

  :   Scale effect into approximately `c(-1, 1)` by dividing by a global
      magnitude (see `unit_mode`, `unit_q`).

  "zscale"

  :   Z-score effect values within each pathway across groups: \\z =
      (x - \mu)/\sigma\\. If there is only one group, this mode
      automatically falls back to `"none"`.

- unit_mode:

  Character; only used when `scale_mode = "unit"`. One of:

  "max"

  :   Divide by the maximum absolute effect across all selected
      pathway/group combinations. This makes the most extreme value
      exactly ±1.

  "robust"

  :   Divide by a high quantile (see `unit_q`) of the absolute effect
      distribution, which reduces the influence of extreme outliers and
      improves contrast among moderately changed pathways.

- select_by:

  Character; metric used to rank pathways for inclusion within each side
  of the mean `effect` distribution. One of `"concordance"` (default),
  `"effect"`, or `"abs_effect"`.

- reorder_by:

  Character; if `reorder = TRUE`, controls how pathways are ordered on
  the y-axis. One of `"effect"` (default; negative to positive mean
  effect) or `"select"` (order by the selection rank).

- unit_q:

  Numeric between 0 and 1; quantile of `|effect|` used as the scaling
  factor when `unit_mode = "robust"`. Default is 0.9.

- reorder:

  Logical; if `TRUE` (default), reorder pathways on the y-axis by their
  mean effect (from most negative to most positive).

- fill_limits:

  Optional numeric vector of length 2 giving limits for the fill scale.
  If `NULL`, limits are set automatically depending on `scale_mode`.

- show_concordance:

  Logical; if TRUE (default), show the concordance metric as text labels
  when available.

- show_robustness:

  Logical; if TRUE, show the robustness metric (column
  `robust_gt_thresh`) as text labels when available. Default is FALSE,
  so robustness is not shown even if computed, unless explicitly
  requested.

## Value

A ggplot2 heatmap object.

## Details

Tiles are filled by the `effect` column, while additional metrics such
as `concordance` and the robustness score (`robust_gt_thresh`, fraction
of strong, sign-consistent cells among all cells in a group) are
displayed as text overlays (labels) if present and non-missing. In the
plot, this robustness is shown simply as `"robustness"`.
