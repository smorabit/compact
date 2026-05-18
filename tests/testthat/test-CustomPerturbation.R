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

  # build KNN graph — not stored in the saved object
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)

  # load TOM directly to avoid working-directory dependency of GetTOM
  tom_env <- new.env()
  load(tom_path, envir = tom_env)
  mods <- GetModules(seurat_obj, "simulation")
  TOM  <- as.matrix(tom_env$consTomDS)
  rownames(TOM) <- mods$gene_name
  colnames(TOM) <- mods$gene_name

  exp <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

  # selected features: top 5 red-module hub genes
  hub_df    <- GetHubGenes(seurat_obj, n = 5, wgcna_name = "simulation")
  sel_genes <- subset(hub_df, module == "red")$gene_name

  n_conn_test <- 20L

  # precompute which genes the connectivity path will select as the neighborhood
  conn_scores      <- colSums(TOM[sel_genes, ])
  candidates_test  <- setdiff(names(sort(conn_scores, decreasing = TRUE)), sel_genes)
  non_hub_test     <- candidates_test[seq_len(min(n_conn_test, length(candidates_test)))]
  module_genes_test <- unique(c(sel_genes, non_hub_test))

  # run knock-in and knock-out once each at file scope
  result_ki <- CustomPerturbation(
    seurat_obj,
    selected_features = sel_genes,
    perturb_dir       = 1,
    perturbation_name = "custom_ki_test",
    graph             = "RNA_nn",
    n_connections     = n_conn_test,
    n_iters           = 3,
    use_velocyto      = FALSE,
    custom_network    = TOM,
    custom_modules    = mods
  )

  result_ko <- CustomPerturbation(
    seurat_obj,
    selected_features = sel_genes,
    perturb_dir       = 0,
    perturbation_name = "custom_ko_test",
    graph             = "RNA_nn",
    n_connections     = n_conn_test,
    n_iters           = 3,
    use_velocyto      = FALSE,
    custom_network    = TOM,
    custom_modules    = mods
  )
}

# ---------------------------------------------------------------------------
# 1. New assay is created with the correct name
# ---------------------------------------------------------------------------

test_that("CustomPerturbation creates a new assay with the given perturbation_name", {
  skip_if_no_data()

  expect_true("custom_ki_test" %in% names(result_ki@assays),
    label = "perturbation assay is present in Seurat object")
})

# ---------------------------------------------------------------------------
# 2. Perturbation assay has the same dimensions as the source assay
# ---------------------------------------------------------------------------

test_that("perturbation assay has same dimensions as the source assay", {
  skip_if_no_data()

  expect_equal(
    dim(result_ki[["custom_ki_test"]]),
    dim(seurat_obj[["RNA"]]),
    label = "perturbation assay dimensions match RNA assay"
  )
})

# ---------------------------------------------------------------------------
# 3. Transition probability graph stored with correct name
# ---------------------------------------------------------------------------

test_that("transition probability graph is stored as {perturbation_name}_tp", {
  skip_if_no_data()

  expect_true("custom_ki_test_tp" %in% names(result_ki@graphs),
    label = "TP graph present with expected name")
})

# ---------------------------------------------------------------------------
# 4. TP graph is square cells x cells
# ---------------------------------------------------------------------------

test_that("transition probability graph is a square cells x cells matrix", {
  skip_if_no_data()

  tp     <- result_ki@graphs[["custom_ki_test_tp"]]
  n_cells <- ncol(seurat_obj)

  expect_equal(nrow(tp), n_cells, label = "TP graph row count matches cell count")
  expect_equal(ncol(tp), n_cells, label = "TP graph col count matches cell count")
})

# ---------------------------------------------------------------------------
# 5. Both counts and data layers are present
# ---------------------------------------------------------------------------

test_that("perturbation assay contains both counts and data layers", {
  skip_if_no_data()

  layers <- SeuratObject::Layers(result_ki[["custom_ki_test"]])

  expect_true("counts" %in% layers, label = "counts layer present")
  expect_true("data"   %in% layers, label = "data (log-normalized) layer present")
})

# ---------------------------------------------------------------------------
# 6. Counts layer has no negative values
# ---------------------------------------------------------------------------

test_that("counts layer of the perturbation assay has no negative values", {
  skip_if_no_data()

  counts <- Seurat::GetAssayData(result_ki, assay = "custom_ki_test", layer = "counts")
  expect_true(min(counts) >= 0, label = "all perturbed counts are non-negative")
})

# ---------------------------------------------------------------------------
# 7. Knock-out sets selected_features counts to zero
# ---------------------------------------------------------------------------

test_that("knock-out (perturb_dir=0) sets selected_features counts to zero", {
  skip_if_no_data()

  ko_counts <- Seurat::GetAssayData(result_ko, assay = "custom_ko_test", layer = "counts")
  expect_true(all(ko_counts[sel_genes, ] == 0),
    label = "selected gene counts are zero after knockout")
})

# ---------------------------------------------------------------------------
# 8. Knock-in raises mean selected_features expression vs baseline
# ---------------------------------------------------------------------------

test_that("knock-in (perturb_dir=1) raises mean selected_features counts vs baseline", {
  skip_if_no_data()

  baseline_mean <- mean(as.numeric(exp[sel_genes, ]))
  ki_counts     <- Seurat::GetAssayData(result_ki, assay = "custom_ki_test", layer = "counts")
  ki_mean       <- mean(as.numeric(ki_counts[sel_genes, ]))

  expect_gt(ki_mean, baseline_mean,
    label = "mean selected gene expression increases after knock-in")
})

# ---------------------------------------------------------------------------
# 9. Genes outside the co-expression neighborhood are unchanged in the knockout
# ---------------------------------------------------------------------------

test_that("genes outside the co-expression neighborhood are unchanged in the knockout", {
  skip_if_no_data()

  non_module_genes <- setdiff(rownames(seurat_obj), module_genes_test)
  set.seed(1)
  check_genes <- sample(non_module_genes, min(200L, length(non_module_genes)))

  ko_counts <- Seurat::GetAssayData(result_ko, assay = "custom_ko_test", layer = "counts")

  expect_equal(
    as.matrix(ko_counts[check_genes, ]),
    as.matrix(exp[check_genes, ]),
    label = "non-neighborhood gene counts unchanged in knockout"
  )
})

# ---------------------------------------------------------------------------
# 10. random_connections path runs and returns correct structure
# ---------------------------------------------------------------------------

test_that("random_connections=TRUE path runs and returns a correctly structured result", {
  skip_if_no_data()

  set.seed(42)
  result_rand <- CustomPerturbation(
    seurat_obj,
    selected_features  = sel_genes,
    perturb_dir        = 0,
    perturbation_name  = "custom_rand_test",
    graph              = "RNA_nn",
    n_connections      = 10L,
    random_connections = TRUE,
    n_iters            = 3,
    use_velocyto       = FALSE,
    custom_network     = TOM,
    custom_modules     = mods
  )

  expect_true("custom_rand_test" %in% names(result_rand@assays),
    label = "random_connections perturbation assay present")
  expect_equal(dim(result_rand[["custom_rand_test"]]), dim(seurat_obj[["RNA"]]),
    label = "random_connections assay has correct dimensions")
})

# ---------------------------------------------------------------------------
# 11. custom_weights path runs and returns correct structure
# ---------------------------------------------------------------------------

test_that("custom_weights path selects neighborhood genes by the supplied column", {
  skip_if_no_data()

  # attach a synthetic weight column so the test doesn't depend on kME names
  mods_w <- mods
  set.seed(7)
  mods_w$test_weight <- runif(nrow(mods_w))

  result_cw <- CustomPerturbation(
    seurat_obj,
    selected_features = sel_genes,
    perturb_dir       = 0,
    perturbation_name = "custom_cw_test",
    graph             = "RNA_nn",
    n_connections     = 10L,
    n_iters           = 3,
    use_velocyto      = FALSE,
    custom_network    = TOM,
    custom_modules    = mods_w,
    custom_weights    = "test_weight"
  )

  expect_true("custom_cw_test" %in% names(result_cw@assays),
    label = "custom_weights perturbation assay present")
  expect_equal(dim(result_cw[["custom_cw_test"]]), dim(seurat_obj[["RNA"]]),
    label = "custom_weights assay has correct dimensions")
})

# ---------------------------------------------------------------------------
# 12. Invalid assay → informative error
# ---------------------------------------------------------------------------

test_that("invalid assay name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    CustomPerturbation(
      seurat_obj,
      selected_features = sel_genes,
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
# 13. Invalid graph → informative error
# ---------------------------------------------------------------------------

test_that("invalid graph name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    CustomPerturbation(
      seurat_obj,
      selected_features = sel_genes,
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
# 14. Non-numeric perturb_dir → informative error
# ---------------------------------------------------------------------------

test_that("non-numeric perturb_dir stops with an informative error", {
  skip_if_no_data()

  expect_error(
    CustomPerturbation(
      seurat_obj,
      selected_features = sel_genes,
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
# 15. Invalid group.by column → informative error
# ---------------------------------------------------------------------------

test_that("invalid group.by column stops with an informative error", {
  skip_if_no_data()

  expect_error(
    CustomPerturbation(
      seurat_obj,
      selected_features = sel_genes,
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
# 16. selected_features not in network → informative error
# ---------------------------------------------------------------------------

test_that("selected_features not in network stops with an informative error", {
  skip_if_no_data()

  expect_error(
    CustomPerturbation(
      seurat_obj,
      selected_features = c(sel_genes, "FAKEGENE_NOTEXIST"),
      perturb_dir       = 1,
      perturbation_name = "test_err",
      graph             = "RNA_nn",
      use_velocyto      = FALSE,
      custom_network    = TOM,
      custom_modules    = mods
    ),
    regexp = "not found in the network",
    label  = "missing feature gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 17. custom_modules missing required columns → informative error
# ---------------------------------------------------------------------------

test_that("custom_modules without required columns stops with an informative error", {
  skip_if_no_data()

  bad_modules <- data.frame(gene = mods$gene_name, color = mods$color)

  expect_error(
    CustomPerturbation(
      seurat_obj,
      selected_features = sel_genes,
      perturb_dir       = 1,
      perturbation_name = "test_err",
      graph             = "RNA_nn",
      use_velocyto      = FALSE,
      custom_network    = TOM,
      custom_modules    = bad_modules
    ),
    regexp = "custom_modules must contain columns",
    label  = "malformed custom_modules gives informative error"
  )
})
