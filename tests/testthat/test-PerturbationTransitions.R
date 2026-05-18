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

  tom_env <- new.env()
  load(tom_path, envir = tom_env)
  mods <- GetModules(seurat_obj, "simulation")
  TOM  <- as.matrix(tom_env$consTomDS)
  rownames(TOM) <- mods$gene_name
  colnames(TOM) <- mods$gene_name

  # run a full ModulePerturbation first so the perturbation assay exists
  seurat_obj <- ModulePerturbation(
    seurat_obj,
    mod               = "red",
    perturb_dir       = 1,
    perturbation_name = "red_up_tp_test",
    graph             = "RNA_nn",
    n_hubs            = 5,
    n_iters           = 3,
    use_velocyto      = FALSE,
    custom_network    = TOM,
    custom_modules    = mods
  )

  red_genes <- subset(mods, module == "red")$gene_name
  n_cells   <- ncol(seurat_obj)
}

# ---------------------------------------------------------------------------
# 1. Function runs without error and returns a Seurat object
# ---------------------------------------------------------------------------

test_that("PerturbationTransitions returns a Seurat object", {
  skip_if_no_data()

  result <- PerturbationTransitions(
    seurat_obj,
    perturbation_name = "red_up_tp_test",
    features          = red_genes,
    graph             = "RNA_nn",
    use_velocyto      = FALSE
  )

  expect_s4_class(result, "Seurat")
})

# ---------------------------------------------------------------------------
# 2. Transition probability graph is stored with correct name
# ---------------------------------------------------------------------------

test_that("TP graph is stored as {perturbation_name}_tp", {
  skip_if_no_data()

  result <- PerturbationTransitions(
    seurat_obj,
    perturbation_name = "red_up_tp_test",
    features          = red_genes,
    graph             = "RNA_nn",
    use_velocyto      = FALSE
  )

  expect_true("red_up_tp_test_tp" %in% names(result@graphs),
    label = "TP graph present with expected name")
})

# ---------------------------------------------------------------------------
# 3. TP graph is square cells x cells
# ---------------------------------------------------------------------------

test_that("TP graph is a square cells x cells matrix", {
  skip_if_no_data()

  result <- PerturbationTransitions(
    seurat_obj,
    perturbation_name = "red_up_tp_test",
    features          = red_genes,
    graph             = "RNA_nn",
    use_velocyto      = FALSE
  )

  tp <- result@graphs[["red_up_tp_test_tp"]]
  expect_equal(nrow(tp), n_cells, label = "TP row count matches cell count")
  expect_equal(ncol(tp), n_cells, label = "TP col count matches cell count")
})

# ---------------------------------------------------------------------------
# 4. TP graph has no NaN values (column-normalization guard)
# ---------------------------------------------------------------------------

test_that("TP graph contains no NaN values", {
  skip_if_no_data()

  result <- PerturbationTransitions(
    seurat_obj,
    perturbation_name = "red_up_tp_test",
    features          = red_genes,
    graph             = "RNA_nn",
    use_velocyto      = FALSE
  )

  tp_vals <- as.numeric(result@graphs[["red_up_tp_test_tp"]])
  expect_false(any(is.nan(tp_vals)), label = "no NaN values in TP graph")
})

# ---------------------------------------------------------------------------
# 5. TP graph column sums are in [0, 1] (valid probability distribution)
# ---------------------------------------------------------------------------

test_that("TP graph column sums are between 0 and 1", {
  skip_if_no_data()

  result <- PerturbationTransitions(
    seurat_obj,
    perturbation_name = "red_up_tp_test",
    features          = red_genes,
    graph             = "RNA_nn",
    use_velocyto      = FALSE
  )

  tp      <- result@graphs[["red_up_tp_test_tp"]]
  cs      <- Matrix::colSums(tp)
  expect_true(all(cs >= -1e-9 & cs <= 1 + 1e-9),
    label = "all column sums are in [0, 1]")
})

# ---------------------------------------------------------------------------
# 6. TP graph has no negative values
# ---------------------------------------------------------------------------

test_that("TP graph contains no negative values", {
  skip_if_no_data()

  result <- PerturbationTransitions(
    seurat_obj,
    perturbation_name = "red_up_tp_test",
    features          = red_genes,
    graph             = "RNA_nn",
    use_velocyto      = FALSE
  )

  tp <- result@graphs[["red_up_tp_test_tp"]]
  expect_true(min(tp) >= 0, label = "all TP values are non-negative")
})

# ---------------------------------------------------------------------------
# 7. Invalid graph → informative error
# ---------------------------------------------------------------------------

test_that("invalid graph name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    PerturbationTransitions(
      seurat_obj,
      perturbation_name = "red_up_tp_test",
      features          = red_genes,
      graph             = "nonexistent_graph",
      use_velocyto      = FALSE
    ),
    regexp = "not found in seurat_obj@graphs",
    label  = "bad graph name gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 8. Invalid observed assay → informative error
# ---------------------------------------------------------------------------

test_that("invalid observed assay name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    PerturbationTransitions(
      seurat_obj,
      perturbation_name = "red_up_tp_test",
      features          = red_genes,
      graph             = "RNA_nn",
      assay             = "nonexistent_assay",
      use_velocyto      = FALSE
    ),
    regexp = "not found in seurat_obj@assays",
    label  = "bad assay name gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 9. Non-existent perturbation assay → informative error
# ---------------------------------------------------------------------------

test_that("non-existent perturbation assay stops with an informative error", {
  skip_if_no_data()

  expect_error(
    PerturbationTransitions(
      seurat_obj,
      perturbation_name = "nonexistent_perturbation",
      features          = red_genes,
      graph             = "RNA_nn",
      use_velocyto      = FALSE
    ),
    regexp = "not found in seurat_obj@assays",
    label  = "missing perturbation assay gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 10. Feature not in observed assay → informative error
# ---------------------------------------------------------------------------

test_that("feature not in observed assay stops with an informative error", {
  skip_if_no_data()

  expect_error(
    PerturbationTransitions(
      seurat_obj,
      perturbation_name = "red_up_tp_test",
      features          = c(red_genes, "FAKEGENE_XYZ"),
      graph             = "RNA_nn",
      use_velocyto      = FALSE
    ),
    regexp = "not found in assay",
    label  = "missing feature gives informative error"
  )
})
