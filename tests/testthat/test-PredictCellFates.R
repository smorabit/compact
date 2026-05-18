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

  n_cells    <- ncol(seurat_obj)
  branch1    <- colnames(seurat_obj)[seurat_obj$branch == "Branch 1"]
  branch3    <- colnames(seurat_obj)[seurat_obj$branch == "Branch 3"]
}

# ===========================================================================
# PredictFates
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. Returns a Seurat object with the expected metadata column
# ---------------------------------------------------------------------------

test_that("PredictFates returns a Seurat object and adds metadata column", {
  skip_if_no_data()

  result <- PredictFates(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    output_name       = "test_fate",
    t_steps           = 20,
    rank_transform    = FALSE,
    return_seurat     = TRUE
  )

  expect_s4_class(result, "Seurat")
  expect_true("test_fate" %in% colnames(result@meta.data),
              label = "fate output column present in metadata")
})

# ---------------------------------------------------------------------------
# 2. Fate scores are in [0, 1] when rank_transform = FALSE
# ---------------------------------------------------------------------------

test_that("PredictFates raw fate scores are in [0, 1]", {
  skip_if_no_data()

  result <- PredictFates(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    output_name       = "test_fate",
    t_steps           = 20,
    rank_transform    = FALSE,
    return_seurat     = FALSE
  )

  scores <- result$fate_score
  expect_true(all(scores >= -1e-9 & scores <= 1 + 1e-9),
              label = "all fate scores in [0, 1]")
})

# ---------------------------------------------------------------------------
# 3. Probability mass is conserved across diffusion steps
# ---------------------------------------------------------------------------

test_that("PredictFates conserves total probability mass", {
  skip_if_no_data()

  result <- PredictFates(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    output_name       = "test_fate",
    t_steps           = 20,
    rank_transform    = FALSE,
    return_seurat     = FALSE
  )

  # for a row-stochastic P, p %*% P preserves the sum exactly
  expect_equal(sum(result$fate_score), 1, tolerance = 1e-6,
               label = "total probability mass sums to 1")
})

# ---------------------------------------------------------------------------
# 4. rank_transform produces ECDF-rank scores in [0, 1]
# ---------------------------------------------------------------------------

test_that("PredictFates rank-transformed scores are in [0, 1]", {
  skip_if_no_data()

  result <- PredictFates(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    output_name       = "test_fate",
    t_steps           = 20,
    rank_transform    = TRUE,
    return_seurat     = FALSE
  )

  scores <- result$fate_score
  expect_true(all(scores >= -1e-9 & scores <= 1 + 1e-9),
              label = "rank-transformed scores in [0, 1]")
})

# ---------------------------------------------------------------------------
# 5. return_seurat = FALSE returns the expected list structure
# ---------------------------------------------------------------------------

test_that("PredictFates with return_seurat = FALSE returns a list", {
  skip_if_no_data()

  result <- PredictFates(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    t_steps           = 10,
    return_seurat     = FALSE
  )

  expect_type(result, "list")
  expect_true("fate_score" %in% names(result))
  expect_true("source_cells" %in% names(result))
  expect_length(result$fate_score, n_cells)
})

# ---------------------------------------------------------------------------
# 6. group.by / group_name interface works
# ---------------------------------------------------------------------------

test_that("PredictFates accepts group.by and group_name instead of source_cells", {
  skip_if_no_data()

  result <- PredictFates(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    group.by          = "branch",
    group_name        = "Branch 1",
    t_steps           = 10,
    return_seurat     = FALSE
  )

  expect_length(result$fate_score, n_cells)
})

# ---------------------------------------------------------------------------
# 7. Invalid graph → informative error
# ---------------------------------------------------------------------------

test_that("PredictFates stops on an invalid KNN graph name", {
  skip_if_no_data()

  expect_error(
    PredictFates(
      seurat_obj,
      perturbation_name = "red_up",
      graph             = "nonexistent_graph",
      source_cells      = branch1
    ),
    regexp = "not found in seurat_obj@graphs"
  )
})

# ---------------------------------------------------------------------------
# 8. Missing perturbation graph → informative error
# ---------------------------------------------------------------------------

test_that("PredictFates stops when the TP graph is absent", {
  skip_if_no_data()

  expect_error(
    PredictFates(
      seurat_obj,
      perturbation_name = "nonexistent_perturbation",
      graph             = "RNA_nn",
      source_cells      = branch1
    ),
    regexp = "not found"
  )
})

# ---------------------------------------------------------------------------
# 9. No source_cells and no group specification → informative error
# ---------------------------------------------------------------------------

test_that("PredictFates stops when neither source_cells nor group info is given", {
  skip_if_no_data()

  expect_error(
    PredictFates(
      seurat_obj,
      perturbation_name = "red_up",
      graph             = "RNA_nn"
    ),
    regexp = "source_cells"
  )
})

# ---------------------------------------------------------------------------
# 10. group_name not in metadata column → informative error
# ---------------------------------------------------------------------------

test_that("PredictFates stops when group_name yields no cells", {
  skip_if_no_data()

  expect_error(
    PredictFates(
      seurat_obj,
      perturbation_name = "red_up",
      graph             = "RNA_nn",
      group.by          = "branch",
      group_name        = "NonexistentBranch"
    ),
    regexp = "No cells found"
  )
})

# ===========================================================================
# PredictPerturbationTime
# ===========================================================================

# ---------------------------------------------------------------------------
# 11. Returns a Seurat object and adds the pseudotime metadata column
# ---------------------------------------------------------------------------

test_that("PredictPerturbationTime returns a Seurat object with metadata", {
  skip_if_no_data()

  result <- PredictPerturbationTime(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    sink_cells        = branch3,
    output_name       = "test_ptime",
    max_iter          = 200,
    return_seurat     = TRUE
  )

  expect_s4_class(result, "Seurat")
  expect_true("test_ptime" %in% colnames(result@meta.data))
})

# ---------------------------------------------------------------------------
# 12. Sink cells are assigned time = 0
# ---------------------------------------------------------------------------

test_that("PredictPerturbationTime assigns sink cells a time of 0", {
  skip_if_no_data()

  result <- PredictPerturbationTime(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    sink_cells        = branch3,
    output_name       = "test_ptime",
    max_iter          = 200,
    return_seurat     = FALSE
  )

  expect_true(all(result[branch3] == 0),
              label = "sink cells have pseudotime = 0")
})

# ---------------------------------------------------------------------------
# 13. All pseudotime values are non-negative
# ---------------------------------------------------------------------------

test_that("PredictPerturbationTime pseudotime values are non-negative", {
  skip_if_no_data()

  result <- PredictPerturbationTime(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    sink_cells        = branch3,
    output_name       = "test_ptime",
    max_iter          = 200,
    return_seurat     = FALSE
  )

  expect_true(all(result >= 0), label = "all pseudotime values are non-negative")
})

# ---------------------------------------------------------------------------
# 14. group.by / group_name interface
# ---------------------------------------------------------------------------

test_that("PredictPerturbationTime accepts group.by and group_name", {
  skip_if_no_data()

  result <- PredictPerturbationTime(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    group.by          = "branch",
    group_name        = "Branch 3",
    max_iter          = 200,
    return_seurat     = FALSE
  )

  expect_length(result, n_cells)
})

# ---------------------------------------------------------------------------
# 15. All cells are sinks → error
# ---------------------------------------------------------------------------

test_that("PredictPerturbationTime stops when all cells are sinks", {
  skip_if_no_data()

  expect_error(
    PredictPerturbationTime(
      seurat_obj,
      perturbation_name = "red_up",
      graph             = "RNA_nn",
      sink_cells        = colnames(seurat_obj)
    ),
    regexp = "All cells"
  )
})

# ---------------------------------------------------------------------------
# 16. Invalid graph → informative error
# ---------------------------------------------------------------------------

test_that("PredictPerturbationTime stops on an invalid KNN graph name", {
  skip_if_no_data()

  expect_error(
    PredictPerturbationTime(
      seurat_obj,
      perturbation_name = "red_up",
      graph             = "nonexistent_graph",
      sink_cells        = branch3
    ),
    regexp = "not found in seurat_obj@graphs"
  )
})

# ---------------------------------------------------------------------------
# 17. Missing perturbation graph → informative error
# ---------------------------------------------------------------------------

test_that("PredictPerturbationTime stops when the TP graph is absent", {
  skip_if_no_data()

  expect_error(
    PredictPerturbationTime(
      seurat_obj,
      perturbation_name = "nonexistent_perturbation",
      graph             = "RNA_nn",
      sink_cells        = branch3
    ),
    regexp = "not found"
  )
})

# ===========================================================================
# PredictAttractors
# ===========================================================================

# ---------------------------------------------------------------------------
# 18. Returns a Seurat object with the attractor score column
# ---------------------------------------------------------------------------

test_that("PredictAttractors returns a Seurat object with metadata", {
  skip_if_no_data()

  result <- PredictAttractors(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    output_name       = "test_attractor",
    return_seurat     = TRUE
  )

  expect_s4_class(result, "Seurat")
  expect_true("test_attractor" %in% colnames(result@meta.data))
})

# ---------------------------------------------------------------------------
# 19. Attractor scores are all non-negative (validates the abs() fix)
# ---------------------------------------------------------------------------

test_that("PredictAttractors scores are all non-negative", {
  skip_if_no_data()

  result <- PredictAttractors(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    output_name       = "test_attractor",
    return_seurat     = FALSE
  )

  scores <- result$attractor_score
  expect_true(all(scores >= -1e-12),
              label = "all attractor scores are non-negative (no negative eigenvector artifacts)")
})

# ---------------------------------------------------------------------------
# 20. Attractor scores form a valid probability distribution (sum to 1)
# ---------------------------------------------------------------------------

test_that("PredictAttractors scores sum to 1", {
  skip_if_no_data()

  result <- PredictAttractors(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    return_seurat     = FALSE
  )

  expect_equal(sum(result$attractor_score), 1, tolerance = 1e-6,
               label = "attractor scores sum to 1")
})

# ---------------------------------------------------------------------------
# 21. return_seurat = FALSE returns list with attractor_score and sink_cells
# ---------------------------------------------------------------------------

test_that("PredictAttractors with return_seurat = FALSE returns correct structure", {
  skip_if_no_data()

  result <- PredictAttractors(
    seurat_obj,
    perturbation_name    = "red_up",
    graph                = "RNA_nn",
    quantile_threshold   = 0.95,
    return_seurat        = FALSE
  )

  expect_type(result, "list")
  expect_true("attractor_score" %in% names(result))
  expect_true("sink_cells" %in% names(result))
  expect_length(result$attractor_score, n_cells)
  expect_true(length(result$sink_cells) > 0, label = "some sink cells identified")
})

# ---------------------------------------------------------------------------
# 22. Invalid graph → informative error
# ---------------------------------------------------------------------------

test_that("PredictAttractors stops on an invalid KNN graph name", {
  skip_if_no_data()

  expect_error(
    PredictAttractors(
      seurat_obj,
      perturbation_name = "red_up",
      graph             = "nonexistent_graph"
    ),
    regexp = "not found in seurat_obj@graphs"
  )
})

# ---------------------------------------------------------------------------
# 23. Missing perturbation graph → informative error
# ---------------------------------------------------------------------------

test_that("PredictAttractors stops when the TP graph is absent", {
  skip_if_no_data()

  expect_error(
    PredictAttractors(
      seurat_obj,
      perturbation_name = "nonexistent_perturbation",
      graph             = "RNA_nn"
    ),
    regexp = "not found"
  )
})

# ===========================================================================
# PredictCommitment
# ===========================================================================

# ---------------------------------------------------------------------------
# 24. Returns a Seurat object with the commitment score column
# ---------------------------------------------------------------------------

test_that("PredictCommitment returns a Seurat object with metadata", {
  skip_if_no_data()

  result <- PredictCommitment(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    sink_cells        = branch3,
    output_name       = "test_commit",
    max_iter          = 200,
    return_seurat     = TRUE
  )

  expect_s4_class(result, "Seurat")
  expect_true("test_commit" %in% colnames(result@meta.data))
})

# ---------------------------------------------------------------------------
# 25. Sink cells have commitment score = 1 (boundary condition)
# ---------------------------------------------------------------------------

test_that("PredictCommitment assigns sink cells a score of 1", {
  skip_if_no_data()

  result <- PredictCommitment(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    sink_cells        = branch3,
    output_name       = "test_commit",
    max_iter          = 200,
    return_seurat     = TRUE
  )

  sink_scores <- result@meta.data[branch3, "test_commit"]
  expect_true(all(sink_scores == 1),
              label = "sink cells have commitment score = 1")
})

# ---------------------------------------------------------------------------
# 26. Source cells have commitment score = 0 (boundary condition)
# ---------------------------------------------------------------------------

test_that("PredictCommitment assigns source cells a score of 0", {
  skip_if_no_data()

  result <- PredictCommitment(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    sink_cells        = branch3,
    output_name       = "test_commit",
    max_iter          = 200,
    return_seurat     = TRUE
  )

  source_scores <- result@meta.data[branch1, "test_commit"]
  expect_true(all(source_scores == 0),
              label = "source cells have commitment score = 0")
})

# ---------------------------------------------------------------------------
# 27. Transient cells have scores strictly in [0, 1]
# ---------------------------------------------------------------------------

test_that("PredictCommitment transient cell scores are in [0, 1]", {
  skip_if_no_data()

  result <- PredictCommitment(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    sink_cells        = branch3,
    max_iter          = 200,
    return_seurat     = FALSE
  )

  transient_cells <- setdiff(colnames(seurat_obj), c(branch1, branch3))
  transient_scores <- result[transient_cells]
  expect_true(all(transient_scores >= -1e-9 & transient_scores <= 1 + 1e-9),
              label = "transient cell scores are in [0, 1]")
})

# ---------------------------------------------------------------------------
# 28. group.by / source_group / sink_group interface
# ---------------------------------------------------------------------------

test_that("PredictCommitment accepts group-based source and sink specification", {
  skip_if_no_data()

  result <- PredictCommitment(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    group.by          = "branch",
    source_group      = "Branch 1",
    sink_group        = "Branch 3",
    max_iter          = 200,
    return_seurat     = FALSE
  )

  expect_length(result, n_cells)
})

# ---------------------------------------------------------------------------
# 29. Overlapping source and sink cells → informative error
# ---------------------------------------------------------------------------

test_that("PredictCommitment stops when source and sink cells overlap", {
  skip_if_no_data()

  shared <- branch1[1:10]

  expect_error(
    PredictCommitment(
      seurat_obj,
      perturbation_name = "red_up",
      graph             = "RNA_nn",
      source_cells      = shared,
      sink_cells        = shared
    ),
    regexp = "overlap"
  )
})

# ---------------------------------------------------------------------------
# 30. Invalid graph → informative error
# ---------------------------------------------------------------------------

test_that("PredictCommitment stops on an invalid KNN graph name", {
  skip_if_no_data()

  expect_error(
    PredictCommitment(
      seurat_obj,
      perturbation_name = "red_up",
      graph             = "nonexistent_graph",
      source_cells      = branch1,
      sink_cells        = branch3
    ),
    regexp = "not found in seurat_obj@graphs"
  )
})

# ---------------------------------------------------------------------------
# 31. Missing perturbation graph → informative error
# ---------------------------------------------------------------------------

test_that("PredictCommitment stops when the TP graph is absent", {
  skip_if_no_data()

  expect_error(
    PredictCommitment(
      seurat_obj,
      perturbation_name = "nonexistent_perturbation",
      graph             = "RNA_nn",
      source_cells      = branch1,
      sink_cells        = branch3
    ),
    regexp = "not found"
  )
})

# ===========================================================================
# Mathematical consistency cross-checks
# ===========================================================================

# ---------------------------------------------------------------------------
# 32. PredictFates direction check: source-branch cells receive higher
#     fate scores than the opposite-branch cells on average
# ---------------------------------------------------------------------------

test_that("PredictFates source cells have higher scores than distant cells", {
  skip_if_no_data()

  result <- PredictFates(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    t_steps           = 50,
    rank_transform    = FALSE,
    return_seurat     = FALSE
  )

  scores <- result$fate_score
  mean_source <- mean(scores[branch1])
  mean_distant <- mean(scores[branch3])

  # cells near the source should retain more probability mass than
  # cells at the far end of a different branch
  expect_gt(mean_source, mean_distant,
            label = "source-branch cells have higher fate scores than distant-branch cells")
})

# ---------------------------------------------------------------------------
# 33. PredictCommitment monotonicity: mean commitment score should be
#     strictly between 0 (source) and 1 (sink) for transient cells
# ---------------------------------------------------------------------------

test_that("PredictCommitment transient mean score is strictly between 0 and 1", {
  skip_if_no_data()

  result <- PredictCommitment(
    seurat_obj,
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    source_cells      = branch1,
    sink_cells        = branch3,
    max_iter          = 500,
    return_seurat     = FALSE
  )

  transient_cells  <- setdiff(colnames(seurat_obj), c(branch1, branch3))
  mean_transient   <- mean(result[transient_cells])

  expect_gt(mean_transient, 0,   label = "mean transient score > 0")
  expect_lt(mean_transient, 1,   label = "mean transient score < 1")
})
