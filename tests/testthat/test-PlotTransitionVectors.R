library(testthat)
library(compact)
library(hdWGCNA)
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
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE)

  tom_env <- new.env()
  load(tom_path, envir = tom_env)
  mods <- GetModules(seurat_obj, "simulation")
  TOM  <- as.matrix(tom_env$consTomDS)
  rownames(TOM) <- mods$gene_name
  colnames(TOM) <- mods$gene_name

  seurat_obj <- ModulePerturbation(
    seurat_obj,
    mod               = "red",
    perturb_dir       = 1,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    n_hubs            = 5,
    n_iters           = 3,
    use_velocyto      = FALSE,
    custom_network    = TOM,
    custom_modules    = mods
  )

  n_cells <- ncol(seurat_obj)
}

# ===========================================================================
# PerturbationVectors
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. Returns a list with the expected components
# ---------------------------------------------------------------------------

test_that("PerturbationVectors returns a named list with ars and arsd", {
  skip_if_no_data()

  result <- PerturbationVectors(
    seurat_obj,
    perturbation_name = "red_up",
    reduction         = "umap",
    use_velocyto      = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("ars", "arsd"), ignore.order = FALSE)
})

# ---------------------------------------------------------------------------
# 2. ars has the right shape and column names
# ---------------------------------------------------------------------------

test_that("PerturbationVectors ars has correct dimensions and column names", {
  skip_if_no_data()

  result <- PerturbationVectors(
    seurat_obj,
    perturbation_name = "red_up",
    reduction         = "umap",
    use_velocyto      = FALSE
  )

  ars <- result$ars
  expect_equal(nrow(ars), n_cells, label = "ars row count matches cell count")
  expect_equal(ncol(ars), 4,       label = "ars has 4 columns")
  expect_equal(colnames(ars), c("x0", "y0", "x1", "y1"),
               label = "ars column names are x0 y0 x1 y1")
})

# ---------------------------------------------------------------------------
# 3. arsd has the right shape and column names
# ---------------------------------------------------------------------------

test_that("PerturbationVectors arsd has correct dimensions and column names", {
  skip_if_no_data()

  result <- PerturbationVectors(
    seurat_obj,
    perturbation_name = "red_up",
    reduction         = "umap",
    use_velocyto      = FALSE
  )

  arsd <- result$arsd
  expect_equal(nrow(arsd), n_cells, label = "arsd row count matches cell count")
  expect_equal(ncol(arsd), 2,       label = "arsd has 2 columns")
  expect_equal(colnames(arsd), c("xd", "yd"),
               label = "arsd column names are xd yd")
})

# ---------------------------------------------------------------------------
# 4. arsd values are capped in [-1, 1] after max_pct normalization
# ---------------------------------------------------------------------------

test_that("PerturbationVectors arsd values are in [-1, 1]", {
  skip_if_no_data()

  result <- PerturbationVectors(
    seurat_obj,
    perturbation_name = "red_up",
    reduction         = "umap",
    use_velocyto      = FALSE
  )

  arsd <- result$arsd
  expect_true(all(arsd >= -1 - 1e-9 & arsd <= 1 + 1e-9),
              label = "all arsd values are in [-1, 1]")
})

# ---------------------------------------------------------------------------
# 5. Output contains no NaN or Inf values
# ---------------------------------------------------------------------------

test_that("PerturbationVectors output contains no NaN or Inf", {
  skip_if_no_data()

  result <- PerturbationVectors(
    seurat_obj,
    perturbation_name = "red_up",
    reduction         = "umap",
    use_velocyto      = FALSE
  )

  ars_vals  <- unlist(result$ars)
  arsd_vals <- unlist(result$arsd)
  expect_false(any(is.nan(ars_vals)),  label = "no NaN in ars")
  expect_false(any(is.nan(arsd_vals)), label = "no NaN in arsd")
  expect_false(any(is.infinite(ars_vals)),  label = "no Inf in ars")
  expect_false(any(is.infinite(arsd_vals)), label = "no Inf in arsd")
})

# ---------------------------------------------------------------------------
# 6. arrow_scale parameter changes vector magnitudes
# ---------------------------------------------------------------------------

test_that("arrow_scale changes vector magnitudes", {
  skip_if_no_data()

  r1 <- PerturbationVectors(seurat_obj, "red_up", use_velocyto = FALSE,
                             arrow_scale = 1)
  r2 <- PerturbationVectors(seurat_obj, "red_up", use_velocyto = FALSE,
                             arrow_scale = 2)

  mag1 <- sqrt(r1$arsd$xd^2 + r1$arsd$yd^2)
  mag2 <- sqrt(r2$arsd$xd^2 + r2$arsd$yd^2)

  # after normalization magnitudes are scaled relative to max_pct, so
  # scale=2 produces vectors at least as large on average
  expect_gte(mean(mag2), mean(mag1) - 1e-9,
             label = "larger arrow_scale produces equal or larger mean magnitude")
})

# ---------------------------------------------------------------------------
# 7. Invalid reduction → informative error
# ---------------------------------------------------------------------------

test_that("PerturbationVectors stops on an invalid reduction", {
  skip_if_no_data()

  expect_error(
    PerturbationVectors(
      seurat_obj,
      perturbation_name = "red_up",
      reduction         = "nonexistent_reduction",
      use_velocyto      = FALSE
    ),
    label = "bad reduction raises an error"
  )
})

# ---------------------------------------------------------------------------
# 8. Missing transition-probability graph → informative error
# ---------------------------------------------------------------------------

test_that("PerturbationVectors stops when the TP graph is absent", {
  skip_if_no_data()

  expect_error(
    PerturbationVectors(
      seurat_obj,
      perturbation_name = "nonexistent_perturbation",
      reduction         = "umap",
      use_velocyto      = FALSE
    ),
    label = "missing TP graph raises an error"
  )
})

# ===========================================================================
# grid_vectors
# ===========================================================================

# ---------------------------------------------------------------------------
# 9. Returns a data.frame with start/end columns
# ---------------------------------------------------------------------------

test_that("grid_vectors returns a data.frame with expected columns", {
  skip_if_no_data()

  result <- PerturbationVectors(
    seurat_obj, "red_up", use_velocyto = FALSE
  )
  emb  <- Seurat::Reductions(seurat_obj, "umap")@cell.embeddings[, 1:2]
  colnames(emb) <- paste0("emb_", 1:2)
  arsd <- result$arsd

  gv <- grid_vectors(emb[rownames(arsd), ], arsd, resolution = 20)

  expect_s3_class(gv, "data.frame")
  expect_true(all(c("start.emb_1", "start.emb_2", "end.xd", "end.yd") %in% colnames(gv)),
              label = "grid_vectors has start and end coordinate columns")
})

# ---------------------------------------------------------------------------
# 10. Higher resolution → at least as many grid cells
# ---------------------------------------------------------------------------

test_that("higher resolution produces more grid arrows", {
  skip_if_no_data()

  result <- PerturbationVectors(
    seurat_obj, "red_up", use_velocyto = FALSE
  )
  emb  <- Seurat::Reductions(seurat_obj, "umap")@cell.embeddings[, 1:2]
  colnames(emb) <- paste0("emb_", 1:2)
  arsd <- result$arsd

  gv_lo <- grid_vectors(emb[rownames(arsd), ], arsd, resolution = 10)
  gv_hi <- grid_vectors(emb[rownames(arsd), ], arsd, resolution = 40)

  expect_gte(nrow(gv_hi), nrow(gv_lo),
             label = "finer resolution yields at least as many arrows")
})

# ---------------------------------------------------------------------------
# 11. group_min filter removes sparse grid cells
# ---------------------------------------------------------------------------

test_that("grid_vectors group_min filters low-count grid cells", {
  skip_if_no_data()

  result <- PerturbationVectors(
    seurat_obj, "red_up", use_velocyto = FALSE
  )
  emb  <- Seurat::Reductions(seurat_obj, "umap")@cell.embeddings[, 1:2]
  colnames(emb) <- paste0("emb_", 1:2)
  arsd <- result$arsd

  gv_strict <- grid_vectors(emb[rownames(arsd), ], arsd,
                             resolution = 25, group_min = 20)
  gv_loose  <- grid_vectors(emb[rownames(arsd), ], arsd,
                             resolution = 25, group_min = 1)

  expect_lte(nrow(gv_strict), nrow(gv_loose),
             label = "stricter group_min yields fewer or equal arrows")
})

# ===========================================================================
# PlotTransitionVectors
# ===========================================================================

# ---------------------------------------------------------------------------
# 12. Returns a ggplot object
# ---------------------------------------------------------------------------

test_that("PlotTransitionVectors returns a ggplot object", {
  skip_if_no_data()

  p <- PlotTransitionVectors(
    seurat_obj,
    perturbation_name = "red_up",
    color.by          = "branch",
    reduction         = "umap",
    use_velocyto      = FALSE
  )

  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# 13. arrow_alpha = FALSE branch also returns a ggplot
# ---------------------------------------------------------------------------

test_that("PlotTransitionVectors with arrow_alpha = FALSE returns a ggplot", {
  skip_if_no_data()

  p <- PlotTransitionVectors(
    seurat_obj,
    perturbation_name = "red_up",
    color.by          = "branch",
    reduction         = "umap",
    arrow_alpha       = FALSE,
    use_velocyto      = FALSE
  )

  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# 14. Invalid reduction → error
# ---------------------------------------------------------------------------

test_that("PlotTransitionVectors stops on an invalid reduction", {
  skip_if_no_data()

  expect_error(
    PlotTransitionVectors(
      seurat_obj,
      perturbation_name = "red_up",
      color.by          = "branch",
      reduction         = "nonexistent_reduction",
      use_velocyto      = FALSE
    ),
    label = "invalid reduction raises an error"
  )
})

# ---------------------------------------------------------------------------
# 15. Invalid color.by → error
# ---------------------------------------------------------------------------

test_that("PlotTransitionVectors stops on an invalid color.by column", {
  skip_if_no_data()

  expect_error(
    PlotTransitionVectors(
      seurat_obj,
      perturbation_name = "red_up",
      color.by          = "nonexistent_column",
      reduction         = "umap",
      use_velocyto      = FALSE
    ),
    label = "invalid color.by raises an error"
  )
})

# ===========================================================================
# VectorFieldCoherence
# ===========================================================================

# ---------------------------------------------------------------------------
# 16. Returns a numeric vector with length equal to cell count
# ---------------------------------------------------------------------------

test_that("VectorFieldCoherence returns a numeric vector of length n_cells", {
  skip_if_no_data()

  coh <- VectorFieldCoherence(
    seurat_obj,
    perturbation_name = "red_up",
    reduction         = "umap",
    graph             = "RNA_nn",
    use_velocyto      = FALSE
  )

  expect_type(coh, "double")
  expect_length(coh, n_cells)
})

# ---------------------------------------------------------------------------
# 17. Values are in the valid cosine similarity range [-1, 1]
# ---------------------------------------------------------------------------

test_that("VectorFieldCoherence values are in [-1, 1]", {
  skip_if_no_data()

  coh <- VectorFieldCoherence(
    seurat_obj,
    perturbation_name = "red_up",
    reduction         = "umap",
    graph             = "RNA_nn",
    use_velocyto      = FALSE
  )

  expect_true(all(coh >= -1 - 1e-9 & coh <= 1 + 1e-9),
              label = "all coherence values are in [-1, 1]")
})

# ---------------------------------------------------------------------------
# 18. weighted = TRUE also returns valid coherence values
# ---------------------------------------------------------------------------

test_that("VectorFieldCoherence with weighted = TRUE returns valid values", {
  skip_if_no_data()

  coh <- VectorFieldCoherence(
    seurat_obj,
    perturbation_name = "red_up",
    reduction         = "umap",
    graph             = "RNA_snn",
    weighted          = TRUE,
    use_velocyto      = FALSE
  )

  expect_type(coh, "double")
  expect_length(coh, n_cells)
  expect_true(all(coh >= -1 - 1e-9 & coh <= 1 + 1e-9),
              label = "weighted coherence values are in [-1, 1]")
})

# ---------------------------------------------------------------------------
# 19. Invalid reduction → informative error
# ---------------------------------------------------------------------------

test_that("VectorFieldCoherence stops on an invalid reduction", {
  skip_if_no_data()

  expect_error(
    VectorFieldCoherence(
      seurat_obj,
      perturbation_name = "red_up",
      reduction         = "nonexistent_reduction",
      graph             = "RNA_nn",
      use_velocyto      = FALSE
    ),
    regexp = "not found in Seurat object",
    label  = "invalid reduction gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 20. Missing TP graph → informative error (tests the fixed validation bug)
# ---------------------------------------------------------------------------

test_that("VectorFieldCoherence stops when the TP graph is absent", {
  skip_if_no_data()

  expect_error(
    VectorFieldCoherence(
      seurat_obj,
      perturbation_name = "nonexistent_perturbation",
      reduction         = "umap",
      graph             = "RNA_nn",
      use_velocyto      = FALSE
    ),
    regexp = "not found in Seurat object graphs",
    label  = "missing TP graph gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 21. Invalid neighborhood graph → informative error
# ---------------------------------------------------------------------------

test_that("VectorFieldCoherence stops on an invalid neighborhood graph", {
  skip_if_no_data()

  expect_error(
    VectorFieldCoherence(
      seurat_obj,
      perturbation_name = "red_up",
      reduction         = "umap",
      graph             = "nonexistent_graph",
      use_velocyto      = FALSE
    ),
    regexp = "not found in Seurat object graphs",
    label  = "invalid neighborhood graph gives informative error"
  )
})
