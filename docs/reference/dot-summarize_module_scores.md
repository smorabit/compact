# Summarize pathway scores across cells or groups

This internal helper summarizes columns in `seurat_obj@meta.data` whose
names match `score_pattern` (e.g. all logFC columns ending in
`"M6_logFC"`). If `group.by` is provided, scores are summarized within
each group; otherwise, a single row summarizing all cells is returned.

## Usage

``` r
.summarize_module_scores(
  seurat_obj,
  score_pattern = "_logFC$",
  group.by = NULL,
  fun = c("mean", "median")
)
```

## Arguments

- seurat_obj:

  A Seurat object with per-cell score columns in `seurat_obj@meta.data`.

- score_pattern:

  Regular expression used to select score columns (e.g. `"_logFC$"`).

- group.by:

  Optional metadata column name used to group cells prior to
  summarization. If `NULL`, all cells are treated as one group.

- fun:

  Aggregation function to use for summarization, one of `"mean"` or
  `"median"`.

## Value

A data frame with one row per group (or a single row if `group.by` is
`NULL`), containing the group label, number of cells, and the aggregated
scores for each selected column.

## Note

Currently unused by AlertSystemScore; kept for potential reuse.
