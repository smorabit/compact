library(testthat)
library(compact)
library(hdWGCNA)
library(Seurat)
library(Matrix)

# ---------------------------------------------------------------------------
# Test fixtures — run once at file scope
# ---------------------------------------------------------------------------

data_path <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/simulation_branch.rds"
tom_path  <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/TOM/sim_TOM.rda"

skip_if_no_data <- function() {
  skip_if(!file.exists(data_path) || !file.exists(tom_path),
          "Test dataset or TOM not found")
}

if (file.exists(data_path) && file.exists(tom_path)) {
  seurat_obj <- readRDS(data_path)

  # build KNN graph — not stored in the saved object, needed for PerturbationTransitions
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)

  # load TOM directly to avoid working-directory dependency of GetTOM
  tom_env <- new.env()
  load(tom_path, envir = tom_env)
  mods <- GetModules(seurat_obj, "simulation")
  TOM  <- as.matrix(tom_env$consTomDS)
  rownames(TOM) <- mods$gene_name
  colnames(TOM) <- mods$gene_name

  exp <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

  # run full pipeline twice (knock-in + knockout) — used by most assertions below
  result_up <- ModulePerturbation(
    seurat_obj,
    mod               = "red",
    perturb_dir       = 1,
    perturbation_name = "red_up_test",
    graph             = "RNA_nn",
    n_hubs            = 5,
    n_iters           = 3,
    use_velocyto      = FALSE,
    custom_network    = TOM,
    custom_modules    = mods
  )

  result_ko <- ModulePerturbation(
    seurat_obj,
    mod               = "red",
    perturb_dir       = 0,
    perturbation_name = "red_ko_test",
    graph             = "RNA_nn",
    n_hubs            = 5,
    n_iters           = 3,
    use_velocyto      = FALSE,
    custom_network    = TOM,
    custom_modules    = mods
  )

  hub_df    <- GetHubGenes(seurat_obj, n = 5, wgcna_name = "simulation")
  hub_genes <- subset(hub_df, module == "red")$gene_name
}

# ---------------------------------------------------------------------------
# 1. New assay is created with correct name
# ---------------------------------------------------------------------------

test_that("ModulePerturbation creates a new assay with the given perturbation_name", {
  skip_if_no_data()

  expect_true("red_up_test" %in% names(result_up@assays),
    label = "perturbation assay is present in Seurat object")
})

# ---------------------------------------------------------------------------
# 2. New assay has the same dimensions as the original RNA assay
# ---------------------------------------------------------------------------

test_that("perturbation assay has same dimensions as the source assay", {
  skip_if_no_data()

  expect_equal(
    dim(result_up[["red_up_test"]]),
    dim(seurat_obj[["RNA"]]),
    label = "perturbation assay dimensions match RNA assay"
  )
})

# ---------------------------------------------------------------------------
# 3. Transition probability graph is stored with correct name
# ---------------------------------------------------------------------------

test_that("transition probability graph is stored as {perturbation_name}_tp", {
  skip_if_no_data()

  expect_true("red_up_test_tp" %in% names(result_up@graphs),
    label = "TP graph present with expected name")
})

# ---------------------------------------------------------------------------
# 4. TP graph is a square cells x cells matrix
# ---------------------------------------------------------------------------

test_that("transition probability graph is a square cells x cells matrix", {
  skip_if_no_data()

  tp <- result_up@graphs[["red_up_test_tp"]]
  n_cells <- ncol(seurat_obj)

  expect_equal(nrow(tp), n_cells, label = "TP graph has correct number of rows")
  expect_equal(ncol(tp), n_cells, label = "TP graph has correct number of columns")
})

# ---------------------------------------------------------------------------
# 5. Both counts and data layers are present in the perturbation assay
# ---------------------------------------------------------------------------

test_that("perturbation assay contains both counts and data layers", {
  skip_if_no_data()

  layers <- SeuratObject::Layers(result_up[["red_up_test"]])

  expect_true("counts" %in% layers, label = "counts layer present")
  expect_true("data"   %in% layers, label = "data (log-normalized) layer present")
})

# ---------------------------------------------------------------------------
# 6. Counts layer contains no negative values
# ---------------------------------------------------------------------------

test_that("counts layer of the perturbation assay has no negative values", {
  skip_if_no_data()

  counts <- Seurat::GetAssayData(result_up, assay = "red_up_test", layer = "counts")
  expect_true(min(counts) >= 0, label = "all perturbed counts are non-negative")
})

# ---------------------------------------------------------------------------
# 7. Knock-out sets hub gene counts to zero in the output assay
# ---------------------------------------------------------------------------

test_that("knock-out (perturb_dir=0) sets hub gene counts to zero", {
  skip_if_no_data()

  ko_counts <- Seurat::GetAssayData(result_ko, assay = "red_ko_test", layer = "counts")
  expect_true(all(ko_counts[hub_genes, ] == 0),
    label = "hub gene counts are zero after knockout")
})

# ---------------------------------------------------------------------------
# 8. Knock-in raises mean hub gene expression relative to baseline
# ---------------------------------------------------------------------------

test_that("knock-in (perturb_dir=1) raises mean hub gene counts vs baseline", {
  skip_if_no_data()

  baseline_mean <- mean(as.numeric(exp[hub_genes, ]))
  ki_counts     <- Seurat::GetAssayData(result_up, assay = "red_up_test", layer = "counts")
  ki_mean       <- mean(as.numeric(ki_counts[hub_genes, ]))

  expect_gt(ki_mean, baseline_mean,
    label = "mean hub expression increases after knock-in")
})

# ---------------------------------------------------------------------------
# 9. Non-hub genes are unchanged in the knockout (perturbation limited to module)
# ---------------------------------------------------------------------------

test_that("genes outside the red module are unchanged in the knockout assay", {
  skip_if_no_data()

  red_genes     <- subset(mods, module == "red")$gene_name
  non_mod_genes <- setdiff(rownames(seurat_obj), red_genes)

  # sample 200 non-module genes for efficiency
  set.seed(1)
  check_genes <- sample(non_mod_genes, 200)

  ko_counts  <- Seurat::GetAssayData(result_ko, assay = "red_ko_test",  layer = "counts")
  baseline   <- exp[check_genes, ]

  expect_equal(
    as.matrix(ko_counts[check_genes, ]),
    as.matrix(baseline),
    label = "non-module gene counts unchanged in knockout"
  )
})

# ---------------------------------------------------------------------------
# 10. Invalid module name → informative error
# ---------------------------------------------------------------------------

test_that("invalid module name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ModulePerturbation(
      seurat_obj,
      mod               = "nonexistent_module",
      perturb_dir       = 1,
      perturbation_name = "test_err",
      graph             = "RNA_nn",
      use_velocyto      = FALSE,
      custom_network    = TOM,
      custom_modules    = mods
    ),
    regexp = "not found",
    label  = "bad module name gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 11. Invalid assay → informative error
# ---------------------------------------------------------------------------

test_that("invalid assay name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ModulePerturbation(
      seurat_obj,
      mod               = "red",
      perturb_dir       = 1,
      perturbation_name = "test_err",
      graph             = "RNA_nn",
      assay             = "nonexistent_assay",
      use_velocyto      = FALSE,
      custom_network    = TOM,
      custom_modules    = mods
    ),
    regexp = "Invalid assay",
    label  = "bad assay name gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 12. Invalid graph → informative error (new check)
# ---------------------------------------------------------------------------

test_that("invalid graph name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ModulePerturbation(
      seurat_obj,
      mod               = "red",
      perturb_dir       = 1,
      perturbation_name = "test_err",
      graph             = "nonexistent_graph",
      use_velocyto      = FALSE,
      custom_network    = TOM,
      custom_modules    = mods
    ),
    regexp = "not found",
    label  = "bad graph name gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 13. Non-numeric perturb_dir → informative error
# ---------------------------------------------------------------------------

test_that("non-numeric perturb_dir stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ModulePerturbation(
      seurat_obj,
      mod               = "red",
      perturb_dir       = "up",
      perturbation_name = "test_err",
      graph             = "RNA_nn",
      use_velocyto      = FALSE,
      custom_network    = TOM,
      custom_modules    = mods
    ),
    regexp = "Invalid choice for perturb_dir",
    label  = "non-numeric perturb_dir gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 14. n_hubs < 1 → informative error (new check)
# ---------------------------------------------------------------------------

test_that("n_hubs < 1 stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ModulePerturbation(
      seurat_obj,
      mod               = "red",
      perturb_dir       = 1,
      perturbation_name = "test_err",
      graph             = "RNA_nn",
      n_hubs            = 0,
      use_velocyto      = FALSE,
      custom_network    = TOM,
      custom_modules    = mods
    ),
    regexp = "n_hubs must be a positive integer",
    label  = "n_hubs=0 gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 15. Invalid group.by column → informative error (new check)
# ---------------------------------------------------------------------------

test_that("non-existent group.by column stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ModulePerturbation(
      seurat_obj,
      mod               = "red",
      perturb_dir       = 1,
      perturbation_name = "test_err",
      graph             = "RNA_nn",
      group.by          = "no_such_column",
      use_velocyto      = FALSE,
      custom_network    = TOM,
      custom_modules    = mods
    ),
    regexp = "not found in seurat_obj@meta.data",
    label  = "invalid group.by column gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 16. group_name → warning that it has no effect
# ---------------------------------------------------------------------------

test_that("providing group_name emits a warning that it is not yet implemented", {
  skip_if_no_data()

  expect_warning(
    ModulePerturbation(
      seurat_obj,
      mod               = "red",
      perturb_dir       = 0,
      perturbation_name = "test_gn_warn",
      graph             = "RNA_nn",
      group_name        = "Branch 1",
      use_velocyto      = FALSE,
      custom_network    = TOM,
      custom_modules    = mods
    ),
    regexp = "not yet implemented",
    label  = "group_name triggers a warning"
  )
})
