library(testthat)
library(compact)
library(hdWGCNA)
library(Seurat)
library(Matrix)

# ---------------------------------------------------------------------------
# Test fixtures — loaded once at file scope
#
# ApplyPropagation works in log-normalized expression space.
# Inputs: log_obs_mod (log baseline), delta_log (log delta, hub rows non-zero),
#         network (gene-gene).
# Output: max(0, log_obs_mod + propagated_delta_log)  — floored at zero.
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
  col_sums   <- Matrix::colSums(exp)

  mods      <- hdWGCNA::GetModules(seurat_obj, "simulation")
  red_genes <- subset(mods, module == "red")$gene_name

  tom_env <- new.env()
  load(tom_path, envir = tom_env)
  TOM <- as.matrix(tom_env$consTomDS)
  rownames(TOM) <- mods$gene_name
  colnames(TOM) <- mods$gene_name
  cur_net <- TOM[red_genes, red_genes]

  hub_df    <- hdWGCNA::GetHubGenes(seurat_obj, n = 5, wgcna_name = "simulation")
  hub_genes <- subset(hub_df, module == "red")$gene_name
  non_hub_genes <- setdiff(red_genes, hub_genes)

  # log-normalized baseline for module genes (using package helper)
  log_obs_mod <- compact:::log_normalize(exp[red_genes, ], col_sums)

  # knock-in (2x): only hub gene rows differ from baseline
  exp_ki <- exp
  exp_ki[hub_genes, ] <- round(exp[hub_genes, ] * 2)
  log_per_ki  <- compact:::log_normalize(exp_ki[red_genes, ], col_sums)
  delta_log_ki <- Matrix::Matrix(log_per_ki - log_obs_mod, sparse = TRUE)

  # knock-down (0.5x): only hub gene rows differ from baseline
  exp_kd <- exp
  exp_kd[hub_genes, ] <- round(exp[hub_genes, ] * 0.5)
  log_per_kd  <- compact:::log_normalize(exp_kd[red_genes, ], col_sums)
  delta_log_kd <- Matrix::Matrix(log_per_kd - log_obs_mod, sparse = TRUE)
}

# ---------------------------------------------------------------------------
# 1. Output dimensions match input log_obs_mod
# ---------------------------------------------------------------------------

test_that("output has same dimensions as log_obs_mod", {
  skip_if_no_data()

  result <- ApplyPropagation(log_obs_mod, delta_log_ki, cur_net, n_iters = 3)

  expect_equal(dim(result), dim(log_obs_mod),
    label = "output dimensions match input")
})

# ---------------------------------------------------------------------------
# 2. Row and column names are preserved
# ---------------------------------------------------------------------------

test_that("output preserves gene (row) and cell (column) names", {
  skip_if_no_data()

  result <- ApplyPropagation(log_obs_mod, delta_log_ki, cur_net, n_iters = 3)

  expect_equal(rownames(result), rownames(log_obs_mod), label = "rownames preserved")
  expect_equal(colnames(result), colnames(log_obs_mod), label = "colnames preserved")
})

# ---------------------------------------------------------------------------
# 3. n_iters < 1 throws a clear error
# ---------------------------------------------------------------------------

test_that("n_iters = 0 throws an informative error", {
  skip_if_no_data()

  expect_error(
    ApplyPropagation(log_obs_mod, delta_log_ki, cur_net, n_iters = 0),
    regexp = "n_iters must be a positive integer",
    label  = "n_iters=0 gives a clear error"
  )
})

test_that("n_iters = -1 throws an informative error", {
  skip_if_no_data()

  expect_error(
    ApplyPropagation(log_obs_mod, delta_log_ki, cur_net, n_iters = -1),
    regexp = "n_iters must be a positive integer",
    label  = "n_iters=-1 gives a clear error"
  )
})

# ---------------------------------------------------------------------------
# 4. Knock-in propagates: non-hub log-delta is positive after propagation
# ---------------------------------------------------------------------------

test_that("knock-in propagates: non-hub genes have positive mean log-delta", {
  skip_if_no_data()

  result <- ApplyPropagation(
    log_obs_mod, delta_log_ki, cur_net,
    n_iters = 3, delta_scale = 0.2, row_normalize = FALSE
  )

  delta_nonhub <- as.numeric(result[non_hub_genes, ] - log_obs_mod[non_hub_genes, ])
  expect_gt(mean(delta_nonhub), 0,
    label = "non-hub mean log-delta is positive after knock-in propagation")
})

# ---------------------------------------------------------------------------
# 5. Knock-down propagates: non-hub log-delta is negative after propagation
#    delta = result - log_obs_mod; the floor clips result to 0 for
#    zero-expressing cells, so delta can be 0 or negative — mean should be < 0
# ---------------------------------------------------------------------------

test_that("knock-down propagates: non-hub genes have negative mean log-delta", {
  skip_if_no_data()

  result <- ApplyPropagation(
    log_obs_mod, delta_log_kd, cur_net,
    n_iters = 3, delta_scale = 0.2, row_normalize = FALSE
  )

  delta_nonhub <- as.numeric(result[non_hub_genes, ] - log_obs_mod[non_hub_genes, ])
  expect_lt(mean(delta_nonhub), 0,
    label = "non-hub mean log-delta (result - baseline) is negative after knock-down")
})

# ---------------------------------------------------------------------------
# 6. Non-hub log-deltas for UP and DOWN are negatively correlated
#    (symmetry — the central property of log-space propagation)
# ---------------------------------------------------------------------------

test_that("non-hub log-deltas for UP and DOWN are negatively correlated", {
  skip_if_no_data()

  result_ki <- ApplyPropagation(
    log_obs_mod, delta_log_ki, cur_net,
    n_iters = 3, delta_scale = 0.2, row_normalize = FALSE
  )
  result_kd <- ApplyPropagation(
    log_obs_mod, delta_log_kd, cur_net,
    n_iters = 3, delta_scale = 0.2, row_normalize = FALSE
  )

  delta_up <- as.numeric(result_ki[non_hub_genes, ] - log_obs_mod[non_hub_genes, ])
  delta_dn <- as.numeric(result_kd[non_hub_genes, ] - log_obs_mod[non_hub_genes, ])

  r <- cor(delta_up, delta_dn)
  cat(sprintf("\n  Non-hub log-delta cor (UP vs DOWN): %.4f\n", r))

  expect_lt(r, 0,
    label = "non-hub log-deltas for UP and DOWN are negatively correlated")
})

# ---------------------------------------------------------------------------
# 7. Non-hub DOWN delta coverage: reaches a meaningful fraction of gene-cell entries
#    count-space propagation produced ~18% non-zero; log-space with floor ~40%+.
#    Zero-expressing cells (log_obs=0) whose clipped result stays at 0 contribute
#    a delta of 0, so 100% coverage is no longer expected, but coverage should
#    still be well above the count-space baseline.
# ---------------------------------------------------------------------------

test_that("knock-down propagates non-zero delta to a meaningful fraction of non-hub entries", {
  skip_if_no_data()

  result_kd <- ApplyPropagation(
    log_obs_mod, delta_log_kd, cur_net,
    n_iters = 3, delta_scale = 0.2, row_normalize = FALSE
  )

  delta_dn    <- result_kd[non_hub_genes, ] - log_obs_mod[non_hub_genes, ]
  pct_nonzero <- mean(as.numeric(delta_dn) != 0)
  cat(sprintf("\n  Pct non-zero non-hub DOWN log-delta: %.1f%%\n", pct_nonzero * 100))

  expect_gt(pct_nonzero, 0.25,
    label = "at least 25% of non-hub DOWN log-delta entries are non-zero (well above count-space ~18%)")
})

# ---------------------------------------------------------------------------
# 8. prune_network=TRUE completes without error and returns valid dimensions
# ---------------------------------------------------------------------------

test_that("prune_network=TRUE completes without error and returns correct dimensions", {
  skip_if_no_data()

  result <- ApplyPropagation(
    log_obs_mod, delta_log_ki, cur_net,
    n_iters = 3, prune_network = TRUE, prune_percentile = 0.95
  )

  expect_equal(dim(result), dim(log_obs_mod), label = "pruned: dimensions correct")
})

# ---------------------------------------------------------------------------
# 9. row_normalize=TRUE completes without error and returns correct dimensions
# ---------------------------------------------------------------------------

test_that("row_normalize=TRUE completes without error and returns correct dimensions", {
  skip_if_no_data()

  result <- ApplyPropagation(
    log_obs_mod, delta_log_ki, cur_net,
    n_iters = 3, row_normalize = TRUE
  )

  expect_equal(dim(result), dim(log_obs_mod), label = "row_norm: dimensions correct")
})

# ---------------------------------------------------------------------------
# 10. More iterations produce a larger propagation effect than fewer
# ---------------------------------------------------------------------------

test_that("more iterations produce a larger propagation effect than fewer", {
  skip_if_no_data()

  result_1 <- ApplyPropagation(
    log_obs_mod, delta_log_ki, cur_net,
    n_iters = 1, delta_scale = 0.2, row_normalize = FALSE
  )
  result_5 <- ApplyPropagation(
    log_obs_mod, delta_log_ki, cur_net,
    n_iters = 5, delta_scale = 0.2, row_normalize = FALSE
  )

  shift_1 <- mean(abs(as.numeric(result_1[non_hub_genes, ] - log_obs_mod[non_hub_genes, ])))
  shift_5 <- mean(abs(as.numeric(result_5[non_hub_genes, ] - log_obs_mod[non_hub_genes, ])))

  cat(sprintf("\n  Mean |delta| at n_iters=1: %.5f, n_iters=5: %.5f\n", shift_1, shift_5))

  expect_gt(shift_5, shift_1,
    label = "5 iterations propagate signal further than 1 iteration")
})

# ---------------------------------------------------------------------------
# 11. Fix 1: output is always >= 0 (no negative simulated log-expression)
# ---------------------------------------------------------------------------

test_that("ApplyPropagation output is non-negative for all parameters", {
  skip_if_no_data()

  # knock-down with parameters that trigger amplification (large delta_scale)
  result_kd <- ApplyPropagation(
    log_obs_mod, delta_log_kd, cur_net,
    n_iters = 3, delta_scale = 0.5, row_normalize = FALSE
  )

  expect_true(min(result_kd) >= 0,
    label = "no negative simulated log-expression after knock-down with delta_scale=0.5")

  # knock-in should also be non-negative
  result_ki <- ApplyPropagation(
    log_obs_mod, delta_log_ki, cur_net,
    n_iters = 3, delta_scale = 0.5, row_normalize = FALSE
  )

  expect_true(min(result_ki) >= 0,
    label = "no negative simulated log-expression after knock-in with delta_scale=0.5")
})

# ---------------------------------------------------------------------------
# 12. Fix 5: warning issued when mean non-hub row sum * delta_scale > 1
#    (hub rows are excluded because they are restored each iteration and
#    cannot amplify; only the non-hub subgraph determines amplification)
# ---------------------------------------------------------------------------

test_that("amplification warning fires when mean non-hub row sum x delta_scale > 1", {
  skip_if_no_data()

  # red module non-hub mean row sum is ~15; delta_scale = 0.5 gives ~7.5 > 1
  expect_warning(
    ApplyPropagation(
      log_obs_mod, delta_log_kd, cur_net,
      n_iters = 3, delta_scale = 0.5, row_normalize = FALSE
    ),
    regexp = "Signal amplification detected",
    label  = "amplification warning is issued when non-hub effective multiplier > 1"
  )
})

test_that("no amplification warning when mean non-hub row sum x delta_scale <= 1", {
  skip_if_no_data()

  # red module non-hub mean row sum is ~15; delta_scale = 0.05 gives ~0.75 < 1
  # row_normalize = FALSE so CheckSignalDecay decay-warning also does not fire
  warnings_seen <- character(0)
  withCallingHandlers(
    ApplyPropagation(
      log_obs_mod, delta_log_kd, cur_net,
      n_iters = 3, delta_scale = 0.05, row_normalize = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(
    any(grepl("Signal amplification", warnings_seen)),
    label = "no amplification warning when effective multiplier <= 1"
  )
})

# ---------------------------------------------------------------------------
# 13. Hub preservation: hub gene expression reflects the primary perturbation
#    exactly, regardless of n_iters or delta_scale. Hub genes are held as a
#    persistent source and restored to delta_initial after every iteration.
# ---------------------------------------------------------------------------

test_that("hub gene expression equals baseline + primary perturbation delta (n_iters=3)", {
  skip_if_no_data()

  result_ki <- ApplyPropagation(
    log_obs_mod, delta_log_ki, cur_net,
    n_iters = 3, delta_scale = 0.2, row_normalize = FALSE
  )

  # for knock-in, hub deltas are positive so the floor does not clip them;
  # result[hub,] must exactly equal log_obs_mod[hub,] + delta_log_ki[hub,]
  expected <- as.matrix(log_obs_mod[hub_genes, ] + delta_log_ki[hub_genes, ])
  actual   <- as.matrix(result_ki[hub_genes, ])

  expect_equal(actual, expected, tolerance = 1e-10,
    label = "hub gene result equals baseline + primary delta after n_iters=3")
})

test_that("hub gene expression is invariant to n_iters (persistent source property)", {
  skip_if_no_data()

  hub_expected <- as.matrix(log_obs_mod[hub_genes, ] + delta_log_ki[hub_genes, ])

  for (ni in c(1L, 3L, 5L)) {
    result <- ApplyPropagation(
      log_obs_mod, delta_log_ki, cur_net,
      n_iters = ni, delta_scale = 0.2, row_normalize = FALSE
    )
    expect_equal(as.matrix(result[hub_genes, ]), hub_expected, tolerance = 1e-10,
      label = paste0("hub result invariant at n_iters=", ni))
  }
})
