# AlertSystem: score pathway gene sets for a given module perturbation

This function implements the AlertSystem pathway scoring workflow for a
single module perturbation. It:

1.  Filters the input gene sets to features present in the Seurat
    object.

2.  Computes UCell scores for each gene set on the baseline assay (e.g.
    `"RNA"`) and on the perturbation assay (e.g. `"M6_up"`).

3.  Computes per-pathway log2 fold-changes, \\ \log_2 \left(
    \frac{\mathrm{post} + \varepsilon}{\mathrm{pre} + \varepsilon}
    \right) \\.

4.  Aggregates logFC per group (if `group.by` is provided) or across all
    cells, and computes group-wise metrics including:

    - *effect*: mean/median logFC,

    - *concordance*: fraction of cells whose logFC sign matches the
      group-level effect sign.

    - *robustness* (optional): fraction of cells whose logFC is both (i)
      above the magnitude threshold (\|logFC\| ≥ `robust_effect_thresh`)
      and (ii) has the same sign as the group-level effect.

5.  Extracts all UCell (pre/post) and logFC columns created during
    scoring and stores them, along with the metrics table and run
    parameters, in `seurat_obj@misc$AlertSystem`.

6.  Optionally saves all AlertSystem results as a single RDS file
    (containing `seurat_obj@misc$AlertSystem`) and exports the
    pathway-level metrics table as a CSV file in `out_dir`.

7.  Restores the original `seurat_obj@meta.data`, so that the object is
    not permanently bloated by intermediate score columns.

## Usage

``` r
AlertSystemScore(
  seurat_obj,
  gene_sets,
  perturbation_name,
  base_tag,
  baseline_assay = "RNA",
  group.by = NULL,
  out_dir = NULL,
  min_genes = 5,
  eps = 1e-10,
  summary_fun = c("mean", "median"),
  return_details = FALSE,
  verbose = TRUE,
  compute_concordance = TRUE,
  compute_robustness = TRUE,
  robust_effect_thresh = 0.01
)
```

## Arguments

- seurat_obj:

  A Seurat object that already contains the perturbation assay generated
  by `ModulePerturbation` (e.g. an assay named `"M6_up"`).

- gene_sets:

  Named list of character vectors, each defining a pathway gene set to
  score.

- perturbation_name:

  Name of the perturbation assay in `seurat_obj` (e.g. `"M6_up"`).

- base_tag:

  Character string used to construct suffixes for baseline, perturbation
  and logFC columns. For example, if `base_tag = "M6"` and
  `baseline_assay = "RNA"`, generated columns will end in `"_M6_RNA"`,
  `"_M6_M6up"`, and `"_M6_logFC"`.

- baseline_assay:

  Name of the baseline assay to use for pre-perturbation UCell scoring
  (default `"RNA"`).

- group.by:

  Optional metadata column name used to group cells for summarization.
  If `NULL`, all cells are treated as a single group named
  `"all_cells"`.

- out_dir:

  Optional directory where AlertSystem outputs will be saved. If
  provided, the function writes a single RDS file containing the full
  `seurat_obj@misc$AlertSystem` list (AlertSystem_info) and a CSV file
  with the pathway-level metrics table. If `NULL`, results are not
  written to disk.

- min_genes:

  Minimum number of genes required to retain a gene set after
  intersecting with the Seurat object's feature set.

- eps:

  Small numeric constant, epsilon, added to pre and post scores in the
  logFC computation to avoid division by zero.

- summary_fun:

  Aggregation function used for group-wise effect computation in the
  metrics table, one of `"mean"` or `"median"`.

- return_details:

  Logical; if `FALSE` (default), the function returns only the updated
  Seurat object with results stored under `seurat_obj@misc$AlertSystem`.
  If `TRUE`, a list is returned containing the Seurat object, per-cell
  scores, metrics table, and any file paths used for saving.

- verbose:

  Logical; if `TRUE` (default), progress messages are printed during
  scoring and summarization. Set to `FALSE` to suppress console output.

- compute_concordance:

  Logical; if `TRUE` (default), compute a concordance score per group ×
  pathway as the fraction of cells whose logFC sign matches the
  group-level effect sign.

- compute_robustness:

  Logical; if `TRUE` (default), a robustness metric is computed per
  group × pathway as the fraction of all cells whose logFC is both (i)
  above the magnitude threshold (\|logFC\| ≥ `robust_effect_thresh`)
  and (ii) has the same sign as the group-level effect.

- robust_effect_thresh:

  Numeric threshold on per-cell logFC used to define "strong responders"
  in the robustness metric `robust_gt_thresh`. For each group × pathway,
  `robust_gt_thresh` is the fraction of all cells whose logFC is
  both (i) above this magnitude threshold (\|logFC\| ≥
  `robust_effect_thresh`) and (ii) has the same sign as the group-level
  effect. If `robust_effect_thresh = 0`, this metric is not computed and
  `robust_gt_thresh` is set to `NA`.

## Value

If `return_details = FALSE` (default), returns the Seurat object with
`seurat_obj@misc$AlertSystem` populated and the original metadata
restored. If `return_details = TRUE`, returns a list with elements:

- seurat:

  The updated Seurat object.

- per_cell:

  Data frame of per-cell UCell and logFC scores (also stored in
  `seurat_obj@misc$AlertSystem$per_cell`).

- metrics:

  Long-format data frame of group-wise metrics (effect, concordance,
  robustness) for each pathway (also stored in
  `seurat_obj@misc$AlertSystem$metrics`).

- info_file:

  Path to the RDS file storing the full `seurat_obj@misc$AlertSystem`
  list (or `NULL` if `out_dir` was not provided).

- metrics_file:

  Path to the CSV file containing the metrics table (or `NULL` if
  `out_dir` was not provided).

## Examples

``` r
if (FALSE) { # \dontrun{
# Load pathway gene sets (local directory or GitHub fallback)
alertdatabasePATH <- "/dataPATH/"
genetype <- "mouse"  # or "human"

gene_sets_mouse <- LoadAlertSystemPathways(
  path = alertdatabasePATH,
  genetype= "mouse",
  database = "MSigDBhallmark"
)

# Basic usage: concordance on, robustness on (default)
seurat_obj <- AlertSystemScore(
  seurat_obj         = seurat_obj,
  gene_sets          = gene_sets_mouse,
  perturbation_name  = "M6_up",
  base_tag           = "M6",
  baseline_assay     = "RNA",
  group.by           = "cell_types_broad",
  out_dir            = tempdir(),   # save AlertSystem_info RDS + metrics CSV
  verbose            = TRUE
)

# Access metrics (effect + concordance + robustness)
head(seurat_obj@misc$AlertSystem$metrics)

# Example with a stricter robustness threshold
res <- AlertSystemScore(
  seurat_obj          = seurat_obj,
  gene_sets           = gene_sets_mouse,
  perturbation_name   = "M6_up",
  base_tag            = "M6",
  baseline_assay      = "RNA",
  group.by            = "cell_types_broad",
  compute_robustness = TRUE,
  robust_effect_thresh  = 0.1,   # only count cells with |logFC| ≥ 0.1
  return_details      = TRUE,
  verbose             = TRUE
)

head(res$metrics)
} # }
```
