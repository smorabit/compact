library(testthat)
library(compact)
library(hdWGCNA)
library(Seurat)

# ---------------------------------------------------------------------------
# Fixtures — loaded once at file scope to avoid repeated UCell scoring
#
# Uses the pre-computed red_up assay already stored in simulation_branch.rds.
# Hub genes of the red module are directly verified to be elevated in the
# red_up data layer (hub genes are 2x perturbed), so the hub_set gene set
# provides a reliable positive-effect signal for biological direction tests.
# ---------------------------------------------------------------------------

data_path <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/simulation_branch.rds"

skip_if_no_data <- function() {
  skip_if(!file.exists(data_path), "Test dataset not found")
}

if (file.exists(data_path)) {
  seurat_obj <- readRDS(data_path)
  orig_meta  <- seurat_obj@meta.data

  mods          <- hdWGCNA::GetModules(seurat_obj, wgcna_name = "simulation")
  hub_df        <- hdWGCNA::GetHubGenes(seurat_obj, n = 5, wgcna_name = "simulation")
  hub_genes     <- subset(hub_df, module == "red")$gene_name
  red_genes     <- subset(mods, module == "red")$gene_name
  non_mod_genes <- setdiff(rownames(seurat_obj), red_genes)

  gene_sets_ok <- list(
    hub_set    = hub_genes,
    nonmod_set = head(non_mod_genes, 30)
  )

  gene_sets_with_tiny <- c(
    gene_sets_ok,
    list(tiny_set = head(rownames(seurat_obj), 3))
  )

  tmp_dir <- tempfile()
  dir.create(tmp_dir)

  # grouped fixture: return_details=TRUE, group.by="Group", out_dir set
  res_grouped <- AlertSystemScore(
    seurat_obj        = seurat_obj,
    gene_sets         = gene_sets_with_tiny,
    perturbation_name = "red_up",
    base_tag          = "red",
    baseline_assay    = "RNA",
    group.by          = "Group",
    out_dir           = tmp_dir,
    min_genes         = 5,
    return_details    = TRUE,
    verbose           = FALSE
  )

  # ungrouped fixture: group.by=NULL (for all_cells and direction tests)
  res_ungrouped <- AlertSystemScore(
    seurat_obj        = seurat_obj,
    gene_sets         = gene_sets_ok,
    perturbation_name = "red_up",
    base_tag          = "red",
    baseline_assay    = "RNA",
    group.by          = NULL,
    min_genes         = 5,
    return_details    = TRUE,
    verbose           = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. Default return is an updated Seurat object
# ---------------------------------------------------------------------------

test_that("AlertSystemScore returns a Seurat object by default", {
  skip_if_no_data()

  obj_out <- AlertSystemScore(
    seurat_obj        = seurat_obj,
    gene_sets         = gene_sets_ok,
    perturbation_name = "red_up",
    base_tag          = "red",
    min_genes         = 5,
    verbose           = FALSE
  )

  expect_true(inherits(obj_out, "Seurat"),
    label = "default return is a Seurat object")
})

# ---------------------------------------------------------------------------
# 2. @misc$AlertSystem is populated with per_cell, metrics, params
# ---------------------------------------------------------------------------

test_that("misc$AlertSystem contains per_cell, metrics, and params", {
  skip_if_no_data()

  misc <- res_grouped$seurat@misc$AlertSystem

  expect_true(!is.null(misc),                      label = "misc$AlertSystem is not NULL")
  expect_true(!is.null(misc$per_cell),             label = "misc$AlertSystem$per_cell exists")
  expect_true(!is.null(misc$metrics),              label = "misc$AlertSystem$metrics exists")
  expect_true(!is.null(misc$params),               label = "misc$AlertSystem$params exists")
})

# ---------------------------------------------------------------------------
# 3. per_cell contains UCell and logFC columns with expected suffix patterns
# ---------------------------------------------------------------------------

test_that("per_cell contains expected pre, post, and logFC columns", {
  skip_if_no_data()

  pc   <- res_grouped$per_cell
  cols <- colnames(pc)

  # suffixes: pre=_red_RNA, post=_red_redup, logfc=_red_logFC
  has_pre   <- any(grepl("_red_RNA$",   cols))
  has_post  <- any(grepl("_red_redup$", cols))
  has_logfc <- any(grepl("_red_logFC$", cols))

  expect_true(has_pre,   label = "per_cell has baseline UCell columns (_red_RNA)")
  expect_true(has_post,  label = "per_cell has perturbation UCell columns (_red_redup)")
  expect_true(has_logfc, label = "per_cell has logFC columns (_red_logFC)")
  expect_equal(nrow(pc), ncol(seurat_obj), label = "per_cell has one row per cell")
})

# ---------------------------------------------------------------------------
# 4. metrics has required columns
# ---------------------------------------------------------------------------

test_that("metrics data frame has required columns", {
  skip_if_no_data()

  m <- res_grouped$metrics

  expect_true(is.data.frame(m),                       label = "metrics is a data frame")
  expect_true("group"            %in% colnames(m),    label = "metrics has 'group' column")
  expect_true("pathway"          %in% colnames(m),    label = "metrics has 'pathway' column")
  expect_true("score_name"       %in% colnames(m),    label = "metrics has 'score_name' column")
  expect_true("n_cells"          %in% colnames(m),    label = "metrics has 'n_cells' column")
  expect_true("effect"           %in% colnames(m),    label = "metrics has 'effect' column")
  expect_true("concordance"      %in% colnames(m),    label = "metrics has 'concordance' column")
  expect_true("robust_gt_thresh" %in% colnames(m),    label = "metrics has 'robust_gt_thresh' column")
})

# ---------------------------------------------------------------------------
# 5. Original metadata is restored after the call
# ---------------------------------------------------------------------------

test_that("seurat_obj@meta.data is restored to original after AlertSystemScore", {
  skip_if_no_data()

  obj_out <- AlertSystemScore(
    seurat_obj        = seurat_obj,
    gene_sets         = gene_sets_ok,
    perturbation_name = "red_up",
    base_tag          = "red",
    min_genes         = 5,
    verbose           = FALSE
  )

  expect_equal(colnames(obj_out@meta.data), colnames(orig_meta),
    label = "metadata column names restored to original")
  expect_equal(nrow(obj_out@meta.data), nrow(orig_meta),
    label = "metadata row count unchanged")
})

# ---------------------------------------------------------------------------
# 6. group.by = NULL produces a single "all_cells" group
# ---------------------------------------------------------------------------

test_that("group.by = NULL produces 'all_cells' group in metrics", {
  skip_if_no_data()

  m <- res_ungrouped$metrics

  expect_true(all(m$group == "all_cells"),
    label = "group column is 'all_cells' when group.by = NULL")
})

# ---------------------------------------------------------------------------
# 7. group.by = "Group" produces n_groups × n_pathways rows in metrics
# ---------------------------------------------------------------------------

test_that("grouped metrics has n_groups × n_retained_pathways rows", {
  skip_if_no_data()

  m <- res_grouped$metrics

  n_groups   <- length(unique(seurat_obj$Group))
  # tiny_set is filtered out (3 genes < min_genes=5), so 2 pathways remain
  n_pathways <- 2L

  expect_equal(nrow(m), n_groups * n_pathways,
    label = "metrics has n_groups × n_pathways rows when group.by is set")
})

# ---------------------------------------------------------------------------
# 8. Gene sets below min_genes are filtered out
# ---------------------------------------------------------------------------

test_that("gene sets with fewer than min_genes genes are dropped from metrics", {
  skip_if_no_data()

  m <- res_grouped$metrics

  expect_false("tiny_set" %in% m$pathway,
    label = "tiny_set (3 genes < min_genes=5) is absent from metrics")
  expect_true("hub_set"    %in% m$pathway, label = "hub_set is present in metrics")
  expect_true("nonmod_set" %in% m$pathway, label = "nonmod_set is present in metrics")
})

# ---------------------------------------------------------------------------
# 9. No gene sets survive filtering → informative error
# ---------------------------------------------------------------------------

test_that("AlertSystemScore errors when no gene sets pass min_genes filter", {
  skip_if_no_data()

  all_tiny <- list(
    s1 = head(rownames(seurat_obj), 2),
    s2 = head(rownames(seurat_obj), 3)
  )

  expect_error(
    AlertSystemScore(
      seurat_obj        = seurat_obj,
      gene_sets         = all_tiny,
      perturbation_name = "red_up",
      base_tag          = "red",
      min_genes         = 5,
      verbose           = FALSE
    ),
    regexp = "No gene sets retained",
    label  = "informative error when no gene sets pass min_genes filter"
  )
})

# ---------------------------------------------------------------------------
# 10. compute_concordance = FALSE → concordance column is all NA
# ---------------------------------------------------------------------------

test_that("compute_concordance = FALSE sets concordance to NA", {
  skip_if_no_data()

  res <- AlertSystemScore(
    seurat_obj          = seurat_obj,
    gene_sets           = gene_sets_ok,
    perturbation_name   = "red_up",
    base_tag            = "red",
    min_genes           = 5,
    compute_concordance = FALSE,
    return_details      = TRUE,
    verbose             = FALSE
  )

  expect_true(all(is.na(res$metrics$concordance)),
    label = "concordance is NA when compute_concordance = FALSE")
})

# ---------------------------------------------------------------------------
# 11. robust_effect_thresh = 0 → robust_gt_thresh column is all NA
# ---------------------------------------------------------------------------

test_that("robust_effect_thresh = 0 sets robust_gt_thresh to NA", {
  skip_if_no_data()

  res <- AlertSystemScore(
    seurat_obj           = seurat_obj,
    gene_sets            = gene_sets_ok,
    perturbation_name    = "red_up",
    base_tag             = "red",
    min_genes            = 5,
    robust_effect_thresh = 0,
    return_details       = TRUE,
    verbose              = FALSE
  )

  expect_true(all(is.na(res$metrics$robust_gt_thresh)),
    label = "robust_gt_thresh is NA when robust_effect_thresh = 0")
})

# ---------------------------------------------------------------------------
# 12. Biological direction: hub_set shows positive effect for red_up
#     Hub genes are directly perturbed 2x → UCell score higher in red_up
#     than in RNA baseline → mean log2FC > 0
# ---------------------------------------------------------------------------

test_that("hub gene set shows positive mean effect for UP perturbation", {
  skip_if_no_data()

  hub_row <- subset(res_ungrouped$metrics, pathway == "hub_set")
  cat(sprintf("\n  hub_set effect (all_cells): %.4f\n", hub_row$effect))

  expect_gt(hub_row$effect, 0,
    label = "hub gene set has positive mean log2FC in red_up vs RNA baseline")
})

# ---------------------------------------------------------------------------
# 13. concordance values are in [0, 1] when computed
# ---------------------------------------------------------------------------

test_that("concordance values are in the valid range [0, 1]", {
  skip_if_no_data()

  conc <- res_grouped$metrics$concordance
  conc <- conc[!is.na(conc)]

  expect_true(all(conc >= 0 & conc <= 1),
    label = "all concordance values are between 0 and 1")
})

# ---------------------------------------------------------------------------
# 14. out_dir creates RDS and CSV output files
# ---------------------------------------------------------------------------

test_that("out_dir causes RDS and CSV files to be written", {
  skip_if_no_data()

  expect_false(is.null(res_grouped$info_file),    label = "info_file path is not NULL")
  expect_false(is.null(res_grouped$metrics_file), label = "metrics_file path is not NULL")
  expect_true(file.exists(res_grouped$info_file),    label = "AlertSystem_info RDS file exists")
  expect_true(file.exists(res_grouped$metrics_file), label = "metrics CSV file exists")
})

# ---------------------------------------------------------------------------
# 15. return_details = TRUE returns a list with the expected elements
# ---------------------------------------------------------------------------

test_that("return_details = TRUE returns a list with seurat, per_cell, metrics, info_file, metrics_file", {
  skip_if_no_data()

  expect_true(is.list(res_grouped),                        label = "return_details result is a list")
  expect_true("seurat"       %in% names(res_grouped),      label = "list has 'seurat' element")
  expect_true("per_cell"     %in% names(res_grouped),      label = "list has 'per_cell' element")
  expect_true("metrics"      %in% names(res_grouped),      label = "list has 'metrics' element")
  expect_true("info_file"    %in% names(res_grouped),      label = "list has 'info_file' element")
  expect_true("metrics_file" %in% names(res_grouped),      label = "list has 'metrics_file' element")
  expect_true(inherits(res_grouped$seurat, "Seurat"),      label = "seurat element is a Seurat object")
  expect_true(is.data.frame(res_grouped$per_cell),         label = "per_cell element is a data frame")
  expect_true(is.data.frame(res_grouped$metrics),          label = "metrics element is a data frame")
})

# ---------------------------------------------------------------------------
# 16. SaveAlertSystem writes an RDS file
# ---------------------------------------------------------------------------

test_that("SaveAlertSystem writes an RDS file that exists on disk", {
  skip_if_no_data()

  save_path <- tempfile(fileext = ".rds")
  SaveAlertSystem(res_grouped$seurat, file = save_path)

  expect_true(file.exists(save_path),
    label = "SaveAlertSystem creates an RDS file")
})

# ---------------------------------------------------------------------------
# 17. LoadAlertSystem round-trip restores misc$AlertSystem
# ---------------------------------------------------------------------------

test_that("LoadAlertSystem round-trip restores misc$AlertSystem correctly", {
  skip_if_no_data()

  save_path <- tempfile(fileext = ".rds")
  SaveAlertSystem(res_grouped$seurat, file = save_path)

  # strip misc$AlertSystem then reload
  stripped <- seurat_obj
  stripped@misc[["AlertSystem"]] <- NULL

  reloaded <- LoadAlertSystem(stripped, file = save_path, verbose = FALSE)

  expect_true(!is.null(reloaded@misc$AlertSystem),
    label = "misc$AlertSystem is present after LoadAlertSystem")
  expect_equal(
    names(reloaded@misc$AlertSystem),
    names(res_grouped$seurat@misc$AlertSystem),
    label = "reloaded AlertSystem has the same named elements"
  )
})
