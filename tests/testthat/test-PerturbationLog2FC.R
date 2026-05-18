library(testthat)
library(compact)
library(Seurat)

# ---------------------------------------------------------------------------
# Test fixtures — loaded once at file scope
# ---------------------------------------------------------------------------

data_path <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/simulation_branch.rds"

skip_if_no_data <- function() {
  skip_if(!file.exists(data_path), "Test dataset not found")
}

if (file.exists(data_path)) {
  seurat_obj <- readRDS(data_path)
  # red_up (knock-in) and red_down (knock-down) are pre-stored in the object
  mods      <- hdWGCNA::GetModules(seurat_obj, "simulation")
  red_genes <- subset(mods, module == "red")$gene_name
  hub_df    <- hdWGCNA::GetHubGenes(seurat_obj, n_hubs = 5, wgcna_name = "simulation")
  hub_genes <- subset(hub_df, module == "red")$gene_name
}

# ---------------------------------------------------------------------------
# 1. Returns a data.frame with expected columns
# ---------------------------------------------------------------------------

test_that("returns a data.frame with gene_name, mean_delta, avg_exp, log2fc columns", {
  skip_if_no_data()

  result <- PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                               wgcna_name = "simulation")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("gene_name", "mean_delta", "avg_exp", "log2fc") %in% colnames(result)),
    label = "all four core columns present")
})

# ---------------------------------------------------------------------------
# 2. Without module: all perturbation assay genes are returned
# ---------------------------------------------------------------------------

test_that("without module argument, all perturbed assay genes are returned", {
  skip_if_no_data()

  result <- PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                               wgcna_name = "simulation")

  expect_equal(sort(result$gene_name), sort(rownames(seurat_obj[["red_up"]])),
    label = "gene list matches perturbation assay rows")
})

# ---------------------------------------------------------------------------
# 3. With module: only that module's genes are returned (scoping bug regression)
# ---------------------------------------------------------------------------

test_that("with module, returns only that module's genes (not all modules)", {
  skip_if_no_data()

  result <- PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                               module = "red", wgcna_name = "simulation")

  expect_equal(nrow(result), length(red_genes),
    label = "row count equals red module size, not all-module gene count")
  expect_true(all(result$gene_name %in% red_genes),
    label = "all returned genes belong to the red module")
})

# ---------------------------------------------------------------------------
# 4. With module: hub column is present
# ---------------------------------------------------------------------------

test_that("with module, output contains a 'hub' column", {
  skip_if_no_data()

  result <- PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                               module = "red", wgcna_name = "simulation")

  expect_true("hub" %in% colnames(result),
    label = "hub column present when module is specified")
})

# ---------------------------------------------------------------------------
# 5. Exactly n_hubs genes flagged as 'hub'
# ---------------------------------------------------------------------------

test_that("exactly n_hubs genes are flagged as 'hub' in the hub column", {
  skip_if_no_data()

  n <- 5
  result <- PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                               module = "red", n_hubs = n, wgcna_name = "simulation")

  expect_equal(sum(result$hub == "hub"), n,
    label = "hub count equals n_hubs")
})

# ---------------------------------------------------------------------------
# 6. Knock-in: mean log2fc of hub genes is positive
# ---------------------------------------------------------------------------

test_that("knock-in (red_up) produces positive mean log2fc for hub genes", {
  skip_if_no_data()

  result <- PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                               module = "red", n_hubs = 5, wgcna_name = "simulation")

  hub_l2fc <- result$log2fc[result$hub == "hub"]
  expect_true(mean(hub_l2fc) > 0,
    label = "mean hub log2fc is positive after knock-in")
})

# ---------------------------------------------------------------------------
# 7. Knock-down: mean log2fc of hub genes is non-positive
# ---------------------------------------------------------------------------

test_that("knock-down (red_down) produces non-positive mean log2fc for hub genes", {
  skip_if_no_data()

  result <- PerturbationLog2FC(seurat_obj, perturbation_name = "red_down",
                               module = "red", n_hubs = 5, wgcna_name = "simulation")

  hub_l2fc <- result$log2fc[result$hub == "hub"]
  expect_true(mean(hub_l2fc) <= 0,
    label = "mean hub log2fc is non-positive after knock-down")
})

# ---------------------------------------------------------------------------
# 8. Knock-in: mean_delta is positive for hub genes
# ---------------------------------------------------------------------------

test_that("knock-in (red_up) produces positive mean_delta for hub genes", {
  skip_if_no_data()

  result <- PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                               module = "red", n_hubs = 5, wgcna_name = "simulation")

  hub_delta <- result$mean_delta[result$hub == "hub"]
  expect_true(mean(hub_delta) > 0,
    label = "mean_delta is positive for hub genes after knock-in")
})

# ---------------------------------------------------------------------------
# 9. Non-negative avg_exp values (normalized expression is >= 0)
# ---------------------------------------------------------------------------

test_that("avg_exp values are all non-negative", {
  skip_if_no_data()

  result <- PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                               wgcna_name = "simulation")

  expect_true(all(result$avg_exp >= 0),
    label = "baseline average expression is non-negative")
})

# ---------------------------------------------------------------------------
# 10. Non-existent perturbation assay → informative error
# ---------------------------------------------------------------------------

test_that("non-existent perturbation assay gives an informative error", {
  skip_if_no_data()

  expect_error(
    PerturbationLog2FC(seurat_obj, perturbation_name = "nonexistent_assay"),
    regexp = "not found",
    label = "bad perturbation assay name gives error"
  )
})

# ---------------------------------------------------------------------------
# 11. Non-existent baseline assay → informative error
# ---------------------------------------------------------------------------

test_that("non-existent baseline assay gives an informative error", {
  skip_if_no_data()

  expect_error(
    PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                       assay_obs = "nonexistent_assay"),
    regexp = "not found",
    label = "bad baseline assay name gives error"
  )
})

# ---------------------------------------------------------------------------
# 12. Non-existent module → informative error
# ---------------------------------------------------------------------------

test_that("non-existent module gives an informative error", {
  skip_if_no_data()

  expect_error(
    PerturbationLog2FC(seurat_obj, perturbation_name = "red_up",
                       module = "nonexistent_module", wgcna_name = "simulation"),
    regexp = "not found",
    label = "bad module name gives error"
  )
})
