library(testthat)
library(compact)
library(hdWGCNA)
library(Seurat)

# ---------------------------------------------------------------------------
# Test fixtures — loaded once at file scope
# ---------------------------------------------------------------------------

data_path <- "/home/groups/singlecell/smorabito/analysis/COMPACT/immunotherapy/NSCLC/data/NSCLC_CD8_hdWGCNA_TFs.rds"

skip_if_no_data <- function() {
  skip_if(!file.exists(data_path), "NSCLC dataset not found")
}

if (file.exists(data_path)) {
  seurat_full <- readRDS(data_path)

  # subset to a single sample for speed; WGCNA/TF data lives in @misc and is preserved
  seurat_obj <- seurat_full[, seurat_full$SampleID == "CNAG_128"]

  # build KNN graph — not stored in the saved object
  seurat_obj <- Seurat::FindNeighbors(
    seurat_obj,
    reduction = "HARMONY",
    dims      = 1:20,
    k.param   = 30,
    annoy.metric = "cosine",
    verbose   = FALSE
  )

  # ETV6 has 100 depth=1 target genes — small enough to be fast, large enough
  # for SparseColDeltaCor (ETV5's 1-gene network causes a C++ out-of-bounds error)
  test_tf <- "ETV6"

  result_ko <- TFPerturbation(
    seurat_obj,
    selected_tf       = test_tf,
    perturb_dir       = 0,
    perturbation_name = "ETV6_ko_test",
    graph             = "originalexp_snn",
    assay             = "originalexp",
    depth             = 1,
    n_iters           = 1,
    use_velocyto      = FALSE
  )

  result_kd <- TFPerturbation(
    seurat_obj,
    selected_tf       = test_tf,
    perturb_dir       = -1,
    perturbation_name = "ETV6_kd_test",
    graph             = "originalexp_snn",
    assay             = "originalexp",
    depth             = 1,
    n_iters           = 1,
    use_velocyto      = FALSE
  )

  exp_base <- Seurat::GetAssayData(seurat_obj, assay = "originalexp", layer = "counts")
}

# ---------------------------------------------------------------------------
# 1. New assay is created with the correct name
# ---------------------------------------------------------------------------

test_that("TFPerturbation creates a new assay with the given perturbation_name", {
  skip_if_no_data()

  expect_true("ETV6_ko_test" %in% names(result_ko@assays),
    label = "perturbation assay present in Seurat object")
})

# ---------------------------------------------------------------------------
# 2. Perturbation assay has the same dimensions as the source assay
# ---------------------------------------------------------------------------

test_that("perturbation assay has same dimensions as the source assay", {
  skip_if_no_data()

  expect_equal(
    dim(result_ko[["ETV6_ko_test"]]),
    dim(seurat_obj[["originalexp"]]),
    label = "perturbation assay dimensions match source assay"
  )
})

# ---------------------------------------------------------------------------
# 3. Transition probability graph stored with correct name
# ---------------------------------------------------------------------------

test_that("transition probability graph is stored as {perturbation_name}_tp", {
  skip_if_no_data()

  expect_true("ETV6_ko_test_tp" %in% names(result_ko@graphs),
    label = "TP graph present with expected name")
})

# ---------------------------------------------------------------------------
# 4. TP graph is a square cells x cells matrix
# ---------------------------------------------------------------------------

test_that("transition probability graph is square cells x cells", {
  skip_if_no_data()

  tp <- result_ko@graphs[["ETV6_ko_test_tp"]]
  n  <- ncol(seurat_obj)

  expect_equal(nrow(tp), n, label = "TP graph row count matches cell count")
  expect_equal(ncol(tp), n, label = "TP graph col count matches cell count")
})

# ---------------------------------------------------------------------------
# 5. Both counts and data layers are present
# ---------------------------------------------------------------------------

test_that("perturbation assay contains both counts and data layers", {
  skip_if_no_data()

  layers <- SeuratObject::Layers(result_ko[["ETV6_ko_test"]])

  expect_true("counts" %in% layers, label = "counts layer present")
  expect_true("data"   %in% layers, label = "data layer present")
})

# ---------------------------------------------------------------------------
# 6. Counts layer has no negative values
# ---------------------------------------------------------------------------

test_that("counts layer of the perturbation assay has no negative values", {
  skip_if_no_data()

  counts <- Seurat::GetAssayData(result_ko, assay = "ETV6_ko_test", layer = "counts")
  expect_true(min(counts) >= 0, label = "all perturbed counts are non-negative")
})

# ---------------------------------------------------------------------------
# 7. Knock-out sets TF counts to zero
# ---------------------------------------------------------------------------

test_that("knock-out (perturb_dir=0) sets TF gene counts to zero", {
  skip_if_no_data()

  ko_counts <- Seurat::GetAssayData(result_ko, assay = "ETV6_ko_test", layer = "counts")
  expect_true(all(ko_counts[test_tf, ] == 0),
    label = "TF counts are zero after knockout")
})

# ---------------------------------------------------------------------------
# 8. Knock-down lowers mean TF expression vs baseline
# ---------------------------------------------------------------------------

test_that("knock-down (perturb_dir=-1) lowers mean TF counts vs baseline", {
  skip_if_no_data()

  baseline_mean <- mean(as.numeric(exp_base[test_tf, ]))
  kd_counts     <- Seurat::GetAssayData(result_kd, assay = "ETV6_kd_test", layer = "counts")
  kd_mean       <- mean(as.numeric(kd_counts[test_tf, ]))

  expect_lt(kd_mean, baseline_mean,
    label = "mean TF expression decreases after knock-down")
})

# ---------------------------------------------------------------------------
# 9. Genes outside the TF network are unchanged in the knockout
# ---------------------------------------------------------------------------

test_that("genes outside the TF network are unchanged in the knockout", {
  skip_if_no_data()

  tf_net   <- GetTFTargetGenes(seurat_obj, selected_tfs = test_tf,
                               depth = 1, target_type = "both",
                               use_regulons = TRUE)
  net_genes <- unique(c(test_tf, tf_net$gene))

  # sample up to 200 non-network genes for efficiency
  non_net_genes <- setdiff(rownames(seurat_obj), net_genes)
  set.seed(1)
  check_genes <- sample(non_net_genes, min(200, length(non_net_genes)))

  ko_counts <- Seurat::GetAssayData(result_ko, assay = "ETV6_ko_test", layer = "counts")

  expect_equal(
    as.matrix(ko_counts[check_genes, ]),
    as.matrix(exp_base[check_genes, ]),
    label = "non-network gene counts unchanged in knockout"
  )
})

# ---------------------------------------------------------------------------
# 10. Invalid TF name → informative error
# ---------------------------------------------------------------------------

test_that("TF not in regulons stops with an informative error", {
  skip_if_no_data()

  expect_error(
    TFPerturbation(
      seurat_obj,
      selected_tf       = "NONEXISTENT_TF",
      perturb_dir       = -1,
      perturbation_name = "test_err",
      graph             = "originalexp_snn",
      assay             = "originalexp",
      use_velocyto      = FALSE
    ),
    regexp = "not found in TF regulons",
    label  = "bad TF name gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 11. Invalid assay → informative error
# ---------------------------------------------------------------------------

test_that("invalid assay name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    TFPerturbation(
      seurat_obj,
      selected_tf       = test_tf,
      perturb_dir       = -1,
      perturbation_name = "test_err",
      graph             = "originalexp_snn",
      assay             = "nonexistent_assay",
      use_velocyto      = FALSE
    ),
    regexp = "Invalid assay",
    label  = "bad assay name gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 12. Invalid graph → informative error
# ---------------------------------------------------------------------------

test_that("invalid graph name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    TFPerturbation(
      seurat_obj,
      selected_tf       = test_tf,
      perturb_dir       = -1,
      perturbation_name = "test_err",
      graph             = "nonexistent_graph",
      assay             = "originalexp",
      use_velocyto      = FALSE
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
    TFPerturbation(
      seurat_obj,
      selected_tf       = test_tf,
      perturb_dir       = "down",
      perturbation_name = "test_err",
      graph             = "originalexp_snn",
      assay             = "originalexp",
      use_velocyto      = FALSE
    ),
    regexp = "Invalid choice for perturb_dir",
    label  = "non-numeric perturb_dir gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 14. Invalid group.by column → informative error
# ---------------------------------------------------------------------------

test_that("invalid group.by column stops with an informative error", {
  skip_if_no_data()

  expect_error(
    TFPerturbation(
      seurat_obj,
      selected_tf       = test_tf,
      perturb_dir       = -1,
      perturbation_name = "test_err",
      graph             = "originalexp_snn",
      assay             = "originalexp",
      group.by          = "no_such_column",
      use_velocyto      = FALSE
    ),
    regexp = "not found in seurat_obj@meta.data",
    label  = "invalid group.by column gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 15. Single-gene TF network → informative error (not a C++ crash)
# ---------------------------------------------------------------------------

test_that("TF with a 1-gene network stops with an informative error", {
  skip_if_no_data()

  # ETV5 has exactly 1 target gene at depth=1; without the size check this
  # would crash SparseColDeltaCor with a C++ out-of-bounds error
  expect_error(
    TFPerturbation(
      seurat_obj,
      selected_tf       = "ETV5",
      perturb_dir       = -1,
      perturbation_name = "test_err",
      graph             = "originalexp_snn",
      assay             = "originalexp",
      depth             = 1,
      use_velocyto      = FALSE
    ),
    regexp = "minimum of 2 target genes",
    label  = "1-gene network gives informative error instead of C++ crash"
  )
})
