library(testthat)
library(compact)
library(Seurat)
library(Matrix)

# ---------------------------------------------------------------------------
# Test fixtures — loaded once at file scope
# ---------------------------------------------------------------------------

data_path <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/simulation_branch.rds"
tom_path  <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/TOM/sim_TOM.rda"

skip_if_no_data <- function() {
  skip_if(!file.exists(data_path) || !file.exists(tom_path),
          "Test dataset or TOM not found")
}

if (file.exists(data_path) && file.exists(tom_path)) {
  seurat_obj <- readRDS(data_path)
  exp        <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

  mods      <- hdWGCNA::GetModules(seurat_obj, "simulation")
  red_genes <- subset(mods, module == "red")$gene_name   # 93 genes

  # load TOM directly to avoid working-directory dependency of GetTOM
  tom_env <- new.env()
  load(tom_path, envir = tom_env)
  TOM <- as.matrix(tom_env$consTomDS)
  rownames(TOM) <- mods$gene_name
  colnames(TOM) <- mods$gene_name
  cur_net <- TOM[red_genes, red_genes]

  # baseline and perturbed expression for the red module
  exp_red  <- exp[red_genes, ]
  hub_df   <- hdWGCNA::GetHubGenes(seurat_obj, n = 5, wgcna_name = "simulation")
  hub_genes <- subset(hub_df, module == "red")$gene_name

  # simulated knock-in: double hub gene expression
  exp_per_ki <- exp_red
  exp_per_ki[hub_genes, ] <- exp_per_ki[hub_genes, ] * 2

  # simulated knock-down: halve hub gene expression (floor to 0)
  exp_per_kd <- exp_red
  exp_per_kd[hub_genes, ] <- Matrix::Matrix(
    pmax(as.matrix(exp_per_kd[hub_genes, ]) / 2, 0), sparse = TRUE
  )
}

# ---------------------------------------------------------------------------
# 1. Output dimensions match input
# ---------------------------------------------------------------------------

test_that("output has same dimensions as the input expression matrices", {
  skip_if_no_data()

  result <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 3
  )

  expect_equal(dim(result), dim(exp_red),
    label = "output dimensions match input")
})

# ---------------------------------------------------------------------------
# 2. Row and column names are preserved
# ---------------------------------------------------------------------------

test_that("output preserves gene (row) and cell (column) names", {
  skip_if_no_data()

  result <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 3
  )

  expect_equal(rownames(result), rownames(exp_red),
    label = "rownames preserved")
  expect_equal(colnames(result), colnames(exp_red),
    label = "colnames preserved")
})

# ---------------------------------------------------------------------------
# 3. No negative values in output
# ---------------------------------------------------------------------------

test_that("output contains no negative values", {
  skip_if_no_data()

  result <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 3
  )

  expect_true(min(result) >= 0, label = "all output values are non-negative")
})

# ---------------------------------------------------------------------------
# 4. Output is integer-valued (due to round() at the end)
# ---------------------------------------------------------------------------

test_that("output values are integers after internal rounding", {
  skip_if_no_data()

  result <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 3
  )

  vals <- as.numeric(result)
  expect_true(all(vals == floor(vals)),
    label = "all output values are integer-valued")
})

# ---------------------------------------------------------------------------
# 5. n_iters < 1 throws a clear error (not a silent wrong result)
# ---------------------------------------------------------------------------

test_that("n_iters = 0 throws an informative error", {
  skip_if_no_data()

  expect_error(
    ApplyPropagation(
      seurat_obj, exp_red, exp_per_ki, cur_net,
      perturb_dir = 1, n_iters = 0
    ),
    regexp = "n_iters must be a positive integer",
    label = "n_iters=0 gives a clear error"
  )
})

test_that("n_iters = -1 throws an informative error", {
  skip_if_no_data()

  expect_error(
    ApplyPropagation(
      seurat_obj, exp_red, exp_per_ki, cur_net,
      perturb_dir = 1, n_iters = -1
    ),
    regexp = "n_iters must be a positive integer",
    label = "n_iters=-1 gives a clear error"
  )
})

# ---------------------------------------------------------------------------
# 6. Signal propagation direction: knock-in of hub genes raises expression
#    of non-hub module genes (which are connected to the hubs via the network)
# ---------------------------------------------------------------------------

test_that("knock-in propagates: non-hub module genes increase in mean expression", {
  skip_if_no_data()

  result <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 3, delta_scale = 0.2
  )

  non_hub_genes <- setdiff(red_genes, hub_genes)
  mean_before   <- mean(as.numeric(exp_red[non_hub_genes, ]))
  mean_after    <- mean(as.numeric(result[non_hub_genes, ]))

  expect_gt(mean_after, mean_before,
    label = "non-hub gene mean expression increases after knock-in propagation")
})

# ---------------------------------------------------------------------------
# 7. apply_ceiling caps all gene expression within the observed baseline max
# ---------------------------------------------------------------------------

test_that("apply_ceiling=TRUE keeps every gene within its observed maximum", {
  skip_if_no_data()

  result_ceil <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 3,
    apply_ceiling = TRUE, ceiling_multiplier = 1.0
  )

  max_obs    <- apply(exp_red, 1, max)
  result_max <- apply(result_ceil, 1, max)

  expect_true(all(result_max <= max_obs + 1e-9),
    label = "no gene exceeds its observed baseline maximum when apply_ceiling=TRUE")
})

test_that("without ceiling, expression can exceed the observed maximum", {
  skip_if_no_data()

  result_noceil <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 3, apply_ceiling = FALSE
  )

  max_obs    <- apply(exp_red, 1, max)
  result_max <- apply(result_noceil, 1, max)

  expect_true(any(result_max > max_obs),
    label = "at least one gene exceeds observed max without ceiling")
})

# ---------------------------------------------------------------------------
# 8. prune_network=TRUE runs without error and returns valid output
# ---------------------------------------------------------------------------

test_that("prune_network=TRUE completes without error and returns valid output", {
  skip_if_no_data()

  result <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 3,
    prune_network = TRUE, prune_percentile = 0.95
  )

  expect_equal(dim(result), dim(exp_red),  label = "pruned: dimensions correct")
  expect_true(min(result) >= 0,            label = "pruned: no negative values")
  expect_true(all(as.numeric(result) == floor(as.numeric(result))),
              label = "pruned: integer-valued output")
})

# ---------------------------------------------------------------------------
# 9. row_normalize=TRUE runs without error and returns valid output
# ---------------------------------------------------------------------------

test_that("row_normalize=TRUE completes without error and returns valid output", {
  skip_if_no_data()

  result <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 3, row_normalize = TRUE
  )

  expect_equal(dim(result), dim(exp_red),  label = "row_norm: dimensions correct")
  expect_true(min(result) >= 0,            label = "row_norm: no negative values")
})

# ---------------------------------------------------------------------------
# 10. More iterations produce greater propagation (larger mean shift)
# ---------------------------------------------------------------------------

test_that("more iterations produce a larger propagation effect than fewer", {
  skip_if_no_data()

  non_hub_genes <- setdiff(red_genes, hub_genes)
  mean_baseline <- mean(as.numeric(exp_red[non_hub_genes, ]))

  result_1 <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 1, delta_scale = 0.2
  )
  result_5 <- ApplyPropagation(
    seurat_obj, exp_red, exp_per_ki, cur_net,
    perturb_dir = 1, n_iters = 5, delta_scale = 0.2
  )

  shift_1 <- mean(as.numeric(result_1[non_hub_genes, ])) - mean_baseline
  shift_5 <- mean(as.numeric(result_5[non_hub_genes, ])) - mean_baseline

  expect_gt(shift_5, shift_1,
    label = "5 iterations propagate signal further than 1 iteration")
})
