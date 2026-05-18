library(testthat)
library(compact)
library(Seurat)
library(Matrix)

# ---------------------------------------------------------------------------
# Test fixtures — loaded once at file scope
# ---------------------------------------------------------------------------

data_path <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/simulation_branch.rds"

skip_if_no_data <- function() {
  skip_if(!file.exists(data_path), "Test dataset not found")
}

seurat_obj <- if (file.exists(data_path)) {
  readRDS(data_path)
} else {
  NULL
}

if (!is.null(seurat_obj)) {
  exp <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

  # five hub genes from the red module (small module, 93 genes, ~74% non-zero)
  hub_df  <- hdWGCNA::GetHubGenes(seurat_obj, n = 5, wgcna_name = "simulation")
  hub_genes <- subset(hub_df, module == "red")$gene_name   # length 5
  one_hub   <- hub_genes[1]                                 # single-gene case
}

# ---------------------------------------------------------------------------
# 1. Knockout (perturb_dir = 0)
# ---------------------------------------------------------------------------

test_that("knockout zeroes targeted features and leaves all other genes unchanged", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = 0,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch"
  )

  # perturbed genes must be all zero
  expect_true(all(result[hub_genes, ] == 0),
    label = "perturbed features are zero after knockout")

  # all other genes must be identical to baseline
  other_genes <- setdiff(rownames(exp), hub_genes)
  expect_equal(
    as.matrix(result[other_genes, ]),
    as.matrix(exp[other_genes, ]),
    label = "non-perturbed genes are unchanged after knockout"
  )
})

# ---------------------------------------------------------------------------
# 2. Output dimensions
# ---------------------------------------------------------------------------

test_that("output has the same dimensions as the input expression matrix", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = 1,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch"
  )

  expect_equal(dim(result), dim(exp),
    label = "output dimensions match input dimensions")
})

# ---------------------------------------------------------------------------
# 3. Row and column names are preserved
# ---------------------------------------------------------------------------

test_that("output preserves gene (row) and cell (column) names from the input", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = 1,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch"
  )

  expect_equal(rownames(result), rownames(exp),
    label = "rownames preserved")
  expect_equal(colnames(result), colnames(seurat_obj),
    label = "colnames preserved")
})

# ---------------------------------------------------------------------------
# 4. No negative counts after knock-down
# ---------------------------------------------------------------------------

test_that("no negative counts exist in output after knock-down", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = -1,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch"
  )

  expect_true(min(result) >= 0,
    label = "all counts are non-negative after knock-down")
})

# ---------------------------------------------------------------------------
# 5. Knock-in increases mean expression of perturbed genes
# ---------------------------------------------------------------------------

test_that("knock-in raises mean expression of the targeted genes", {
  skip_if_no_data()
  set.seed(42)

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = 2,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch"
  )

  mean_before <- mean(as.numeric(exp[hub_genes, ]))
  mean_after  <- mean(as.numeric(result[hub_genes, ]))

  expect_gt(mean_after, mean_before,
    label = "mean expression of perturbed genes increases after knock-in")
})

# ---------------------------------------------------------------------------
# 6. Knock-down decreases mean expression of perturbed genes
# ---------------------------------------------------------------------------

test_that("knock-down lowers mean expression of the targeted genes", {
  skip_if_no_data()
  set.seed(42)

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = -1,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch"
  )

  mean_before <- mean(as.numeric(exp[hub_genes, ]))
  mean_after  <- mean(as.numeric(result[hub_genes, ]))

  expect_lt(mean_after, mean_before,
    label = "mean expression of perturbed genes decreases after knock-down")
})

# ---------------------------------------------------------------------------
# 7. Single-gene perturbation — dimensions and names
# ---------------------------------------------------------------------------

test_that("single-gene perturbation returns correct dimensions and names", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features    = one_hub,
    perturb_dir = 1,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch"
  )

  expect_equal(dim(result), dim(exp),
    label = "single-gene: output dimensions match input")
  expect_equal(rownames(result), rownames(exp),
    label = "single-gene: rownames preserved")
  expect_equal(colnames(result), colnames(seurat_obj),
    label = "single-gene: colnames preserved")
  expect_true(min(result) >= 0,
    label = "single-gene: no negative counts")
})

# ---------------------------------------------------------------------------
# 8. Single-gene knockout zeroes only the targeted gene
# ---------------------------------------------------------------------------

test_that("single-gene knockout zeroes only the one targeted gene", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features    = one_hub,
    perturb_dir = 0,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch"
  )

  expect_true(all(result[one_hub, ] == 0),
    label = "single-gene knockout: targeted gene is zero")

  other_genes <- setdiff(rownames(exp), one_hub)
  expect_equal(
    as.matrix(result[other_genes, ]),
    as.matrix(exp[other_genes, ]),
    label = "single-gene knockout: all other genes unchanged"
  )
})

# ---------------------------------------------------------------------------
# 9. NULL group.by defaults gracefully (no error)
# ---------------------------------------------------------------------------

test_that("NULL group.by runs without error, treating all cells as one group", {
  skip_if_no_data()

  expect_no_error(
    ApplyPerturbation(
      seurat_obj, exp,
      features    = hub_genes,
      perturb_dir = 0,
      cells_use   = colnames(seurat_obj),
      group.by    = NULL
    )
  )
})

# ---------------------------------------------------------------------------
# 10. NULL cells_use defaults to all cells
# ---------------------------------------------------------------------------

test_that("NULL cells_use defaults to all cells and produces correct dimensions", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = 0,
    cells_use   = NULL,
    group.by    = "branch"
  )

  expect_equal(ncol(result), ncol(seurat_obj),
    label = "NULL cells_use: output has all cells")
  expect_equal(colnames(result), colnames(seurat_obj),
    label = "NULL cells_use: cell order matches seurat_obj")
})

# ---------------------------------------------------------------------------
# 11. n_workers > 1: parallel path produces same dimensions and non-negative
#     counts as the serial path (correctness, not exact value equality since
#     ZINB sampling is stochastic)
# ---------------------------------------------------------------------------

test_that("n_workers=2 parallel path returns same dimensions and valid counts as serial", {
  skip_if_no_data()
  set.seed(1)

  result_par <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = 1,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch",
    n_workers   = 2
  )

  expect_equal(dim(result_par), dim(exp),
    label = "parallel: output dimensions match input")
  expect_equal(rownames(result_par), rownames(exp),
    label = "parallel: rownames preserved")
  expect_equal(colnames(result_par), colnames(seurat_obj),
    label = "parallel: colnames preserved")
  expect_true(min(result_par) >= 0,
    label = "parallel: no negative counts")
})

# ---------------------------------------------------------------------------
# 12. n_workers > 1 knockout: same deterministic result as serial knockout
# ---------------------------------------------------------------------------

test_that("n_workers=2 knockout gives identical result to serial knockout", {
  skip_if_no_data()

  result_serial <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = 0,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch",
    n_workers   = 1
  )

  result_par <- ApplyPerturbation(
    seurat_obj, exp,
    features    = hub_genes,
    perturb_dir = 0,
    cells_use   = colnames(seurat_obj),
    group.by    = "branch",
    n_workers   = 2
  )

  expect_equal(result_par, result_serial,
    label = "parallel knockout is bit-for-bit identical to serial knockout")
})

# ===========================================================================
# Multiplicative mode tests
# ===========================================================================

# ---------------------------------------------------------------------------
# 13. Multiplicative UP: dimensions and names preserved
# ---------------------------------------------------------------------------

test_that("multiplicative UP preserves output dimensions and names", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features     = hub_genes,
    perturb_dir  = 2,
    perturb_mode = "multiplicative"
  )

  expect_equal(dim(result), dim(exp),
    label = "multiplicative UP: dimensions match input")
  expect_equal(rownames(result), rownames(exp),
    label = "multiplicative UP: rownames preserved")
  expect_equal(colnames(result), colnames(seurat_obj),
    label = "multiplicative UP: colnames preserved")
})

# ---------------------------------------------------------------------------
# 14. Multiplicative UP: no negative values, integer-valued output
# ---------------------------------------------------------------------------

test_that("multiplicative UP produces non-negative integer-valued counts", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features     = hub_genes,
    perturb_dir  = 2,
    perturb_mode = "multiplicative"
  )

  vals <- as.numeric(result)
  expect_true(min(vals) >= 0,
    label = "multiplicative UP: all values are non-negative")
  expect_true(all(vals == floor(vals)),
    label = "multiplicative UP: all values are integer-valued")
})

# ---------------------------------------------------------------------------
# 15. Multiplicative UP: fold change is exactly applied to hub genes
# ---------------------------------------------------------------------------

test_that("multiplicative UP applies fold change exactly to hub genes", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features     = hub_genes,
    perturb_dir  = 2,
    perturb_mode = "multiplicative"
  )

  expected <- round(as.matrix(exp[hub_genes, ]) * 2)
  observed <- as.matrix(result[hub_genes, ])

  expect_equal(observed, expected,
    label = "multiplicative UP: hub gene counts equal round(baseline * fold_change)")
})

# ---------------------------------------------------------------------------
# 16. Multiplicative UP: non-hub genes are exactly unchanged
# ---------------------------------------------------------------------------

test_that("multiplicative UP leaves non-hub genes exactly unchanged", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features     = hub_genes,
    perturb_dir  = 2,
    perturb_mode = "multiplicative"
  )

  other_genes <- setdiff(rownames(exp), hub_genes)
  expect_equal(
    as.matrix(result[other_genes, ]),
    as.matrix(exp[other_genes, ]),
    label = "multiplicative UP: non-hub genes are unchanged"
  )
})

# ---------------------------------------------------------------------------
# 17. Multiplicative UP: cell-specificity — cells with zero baseline remain zero
#     (unlike ZINB, which can assign counts to zero-expressing cells)
# ---------------------------------------------------------------------------

test_that("multiplicative UP keeps zero-expressing cells at zero for each hub gene", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features     = hub_genes,
    perturb_dir  = 2,
    perturb_mode = "multiplicative"
  )

  for(g in hub_genes){
    zero_cells <- which(as.numeric(exp[g, ]) == 0)
    if(length(zero_cells) > 0){
      expect_true(
        all(as.numeric(result[g, zero_cells]) == 0),
        label = paste0("multiplicative UP: zero-expression cells remain zero for ", g)
      )
    }
  }
})

# ---------------------------------------------------------------------------
# 18. Multiplicative UP: mean expression of hub genes increases
# ---------------------------------------------------------------------------

test_that("multiplicative UP raises mean expression of hub genes", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features     = hub_genes,
    perturb_dir  = 2,
    perturb_mode = "multiplicative"
  )

  mean_before <- mean(as.numeric(exp[hub_genes, ]))
  mean_after  <- mean(as.numeric(result[hub_genes, ]))

  expect_gt(mean_after, mean_before,
    label = "multiplicative UP: mean hub gene expression increases")
})

# ---------------------------------------------------------------------------
# 19. Multiplicative DOWN (0.5): mean expression of hub genes decreases,
#     no negatives, non-hub genes unchanged
# ---------------------------------------------------------------------------

test_that("multiplicative DOWN (0.5) lowers hub expression and leaves other genes intact", {
  skip_if_no_data()

  result <- ApplyPerturbation(
    seurat_obj, exp,
    features     = hub_genes,
    perturb_dir  = 0.5,
    perturb_mode = "multiplicative"
  )

  mean_before <- mean(as.numeric(exp[hub_genes, ]))
  mean_after  <- mean(as.numeric(result[hub_genes, ]))

  expect_lt(mean_after, mean_before,
    label = "multiplicative DOWN: mean hub expression decreases")
  expect_true(min(as.numeric(result)) >= 0,
    label = "multiplicative DOWN: no negative values")

  other_genes <- setdiff(rownames(exp), hub_genes)
  expect_equal(
    as.matrix(result[other_genes, ]),
    as.matrix(exp[other_genes, ]),
    label = "multiplicative DOWN: non-hub genes unchanged"
  )
})

# ---------------------------------------------------------------------------
# 20. Multiplicative mode: negative perturb_dir stops with an informative error
# ---------------------------------------------------------------------------

test_that("multiplicative mode with negative perturb_dir stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ApplyPerturbation(
      seurat_obj, exp,
      features     = hub_genes,
      perturb_dir  = -1,
      perturb_mode = "multiplicative"
    ),
    regexp = "positive fold change",
    label  = "negative perturb_dir in multiplicative mode gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 21. Invalid perturb_mode stops with a match.arg error
# ---------------------------------------------------------------------------

test_that("invalid perturb_mode stops with a match.arg error", {
  skip_if_no_data()

  expect_error(
    ApplyPerturbation(
      seurat_obj, exp,
      features     = hub_genes,
      perturb_dir  = 2,
      perturb_mode = "additive"
    ),
    label = "invalid perturb_mode triggers match.arg error"
  )
})
