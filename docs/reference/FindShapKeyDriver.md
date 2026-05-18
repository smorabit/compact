# Identify Key Driver Genes Using SHAP and XGBoost

This function trains an XGBoost classifier to distinguish either:

- a specified case group against all other cells in a metadata column
  (`set_case` vs "rest"), or

- two specific groups within that column (`ident.1` vs `ident.2`,
  pairwise mode),

using gene expression data from a Seurat object.

## Usage

``` r
FindShapKeyDriver(
  seurat_obj,
  conditions,
  set_case = NULL,
  ident.1 = NULL,
  ident.2 = NULL,
  top_n = 50,
  max_depth = 4,
  eta = 0.1,
  nrounds = 50,
  nthread = 2,
  variable_genes = NULL,
  out_dir = NULL,
  mode = c("full", "split", "cv"),
  train_fraction = 0.8,
  nfold = 5,
  seed = 123,
  expr_layer = "data"
)
```

## Arguments

- seurat_obj:

  A `Seurat` object with normalized gene expression and metadata.

- conditions:

  Column name in `seurat_obj@meta.data` indicating the condition labels.

- set_case:

  Character. Case label to compare against all other values in
  `conditions` (case vs rest mode). Must be `NULL` if `ident.1` and
  `ident.2` are used.

- ident.1:

  Character. First group label in `conditions` for pairwise comparison
  (treated as the positive class, label 1). Must be used together with
  `ident.2` and `set_case` must be `NULL`.

- ident.2:

  Character. Second group label in `conditions` for pairwise comparison
  (treated as the negative class, label 0). Must be used together with
  `ident.1` and `set_case` must be `NULL`.

- top_n:

  Integer. Number of top SHAP-ranked genes to return (default: 50).

- max_depth:

  Integer. Maximum depth of each XGBoost tree (default: 4).

- eta:

  Numeric. Learning rate (shrinkage) for XGBoost (default: 0.1).

- nrounds:

  Integer. Number of boosting rounds for XGBoost (default: 50).

- nthread:

  Integer. Number of CPU threads to use for XGBoost (default: 2).

- variable_genes:

  Optional character vector of genes to use. If `NULL`, uses
  `VariableFeatures(seurat_obj)`.

- out_dir:

  Optional directory to save results. If `NULL`, results are not saved
  to disk.

- mode:

  Character. Training mode: `"full"`, `"split"`, or `"cv"` (default:
  `"full"`).

- train_fraction:

  Numeric. Fraction of cells used for training in `"split"` mode
  (default: 0.8).

- nfold:

  Integer. Number of folds for cross-validation (default: 5).

- seed:

  Integer. Random seed for reproducibility (default: 123).

- expr_layer:

  Character. Expression layer/slot to use when building the feature
  matrix. In Seurat v4 this is passed as `slot`; in Seurat v5 this is
  passed as `layer`. Common values are `"data"`, `"counts"`, or
  `"scale.data"` (default: `"data"`).

## Value

A modified `Seurat` object with SHAP results stored in
`seurat_obj@misc$shap`, including:

- `model`: Trained XGBoost model (only in `"full"` and `"split"` modes;
  not stored in `"cv"`).

- `shap_result`: Raw SHAP result object (only in `"full"` and `"split"`
  modes; not stored in `"cv"`).

- `shap_long`: Long-format SHAP values (per-cell, per-gene).

- `shap_summary`: Mean absolute SHAP values per gene.

- `key_drivers`: Top `top_n` SHAP-ranked genes.

- `variable_genes`: Genes used in the model.

- `test_auc`: AUC on held-out test set (only in `"split"` mode; `NA`
  otherwise).

- `auc_per_fold`: Vector of AUC scores per fold (only in `"cv"` mode).

- `mean_auc`: Mean AUC across folds (only in `"cv"` mode; `NA` if AUC
  not computed).

- `mode`: Training mode used.

- `comparison`: A list describing the comparison type and labels, with
  fields: `type` (either `"case_vs_rest"` or `"pairwise"`) and
  `set_case` or `ident.1`/`ident.2` accordingly.

- `params`: A list of key settings used (e.g., `max_depth`, `eta`,
  `nrounds`, `nthread`, `expr_layer`).

- `package_versions`: A named list of package versions used (e.g.,
  `xgboost`, `SHAPforxgboost`).

This implementation enforces exact package versions for reproducibility.
It will stop with an error unless `xgboost == "1.7.11.1"` and
`SHAPforxgboost == "0.1.3"` are installed in the active R library.

## Details

SHAP (SHapley Additive exPlanations) values are computed to quantify the
contribution of each gene to the model predictions. The top SHAP-ranked
genes are returned as key driver candidates.

Three model training modes are supported via the `mode` parameter:

- `"full"`: Train on the entire dataset (no test set, no
  cross-validation).

- `"split"`: Randomly split data into training and testing subsets.

- `"cv"`: Perform k-fold cross-validation.

In `"split"` mode, test set AUC is computed and stored. In `"cv"` mode,
per-fold AUCs and mean AUC are computed. In `"full"` mode, AUC is not
computed to avoid overestimation.

SHAP results and metadata are stored in `seurat_obj@misc$shap`.
Optionally, results are written to disk under `out_dir`, with all files
documented in an `info.txt`.

Exactly one of the following must be supplied:

- `set_case`: for a "case vs rest" comparison, or

- `ident.1` and `ident.2`: for a pairwise comparison between two groups.

## Output Files (if `out_dir` is provided)

- `model.txt`: XGBoost model summary (not in `"cv"` mode).

- `model.rds`: Serialized XGBoost model object (not in `"cv"` mode).

- `shap_result.rds`: Output of
  [`SHAPforxgboost::shap.values()`](https://rdrr.io/pkg/SHAPforxgboost/man/shap.values.html)
  (not in `"cv"` mode).

- `shap_summary.txt`: Mean absolute SHAP value per gene.

- `variable_genes.txt`: List of genes used for model training.

- `key_drivers.txt`: Top SHAP-ranked genes.

- `shap_long.rds`: Full SHAP values in long format.

- `auc.txt`: AUC from held-out test set (only in `"split"` mode, if
  computed).

- `auc_per_fold.txt`: AUC per fold (only in `"cv"` mode, if computed).

- `info.txt`: Summary of saved outputs.

## Examples

``` r
# Case vs rest: AD vs all other diagnoses
seurat_obj <- FindShapKeyDriver(
  seurat_obj,
  conditions = "Diagnosis",
  set_case = "AD",
  top_n = 30
)
#> Error: Package version mismatch for 'xgboost': required 1.7.11.1 but found 3.2.0.1.
#> Fix by installing the required version in this env.

# Case vs rest in split mode with custom output directory
seurat_obj <- FindShapKeyDriver(
  seurat_obj,
  conditions = "Diagnosis",
  set_case = "AD",
  mode = "split",
  out_dir = "results/shap/"
)
#> Error: Package version mismatch for 'xgboost': required 1.7.11.1 but found 3.2.0.1.
#> Fix by installing the required version in this env.

# Case vs rest in 5-fold CV mode
seurat_obj <- FindShapKeyDriver(
  seurat_obj,
  conditions = "Diagnosis",
  set_case = "AD",
  mode = "cv",
  nfold = 5
)
#> Error: Package version mismatch for 'xgboost': required 1.7.11.1 but found 3.2.0.1.
#> Fix by installing the required version in this env.

# Pairwise comparison: AD vs PiD only
seurat_obj <- FindShapKeyDriver(
  seurat_obj,
  conditions = "Diagnosis",
  ident.1 = "AD",
  ident.2 = "PiD",
  mode = "cv",
  nfold = 5
)
#> Error: Package version mismatch for 'xgboost': required 1.7.11.1 but found 3.2.0.1.
#> Fix by installing the required version in this env.
```
