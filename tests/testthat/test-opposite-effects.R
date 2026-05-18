library(testthat)
library(compact)
library(hdWGCNA)
library(Seurat)
library(Matrix)

# ---------------------------------------------------------------------------
# Fixtures — run once at file scope
#
# Runs ModulePerturbation with perturb_mode = "multiplicative" for both
# UP (perturb_dir = 2, i.e. 2x) and DOWN (perturb_dir = 0.5, i.e. 0.5x)
# on the red module. Tests below verify that the two perturbations push
# gene expression in opposite directions.
#
# Note on the TP correlation metric: even when log-normalized deltas are
# correctly anticorrelated (UP vs DOWN), the SparseColDeltaCor-based TP
# computation can still yield positively correlated TP matrices because
# the softmax normalization compresses cell-to-cell differences in the
# transition probabilities. The log-normalized delta correlation is therefore
# a more direct and reliable metric for opposite effects — it measures
# what PerturbationTransitions actually receives as input.
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
  mods <- hdWGCNA::GetModules(seurat_obj, "simulation")
  TOM  <- as.matrix(tom_env$consTomDS)
  rownames(TOM) <- mods$gene_name
  colnames(TOM) <- mods$gene_name

  # row_normalize = FALSE preserves propagation signal in downstream genes;
  # n_iters = 1 + delta_scale = 0.05 keeps the signal small enough to avoid
  # floor saturation for DOWN (which would make UP and DOWN look similar)
  seurat_obj <- ModulePerturbation(
    seurat_obj,
    mod               = "red",
    perturb_dir       = 2,
    perturb_mode      = "multiplicative",
    perturbation_name = "red_up",
    graph             = "RNA_nn",
    n_hubs            = 5,
    n_iters           = 1,
    delta_scale       = 0.05,
    row_normalize     = FALSE,
    use_velocyto      = FALSE,
    custom_network    = TOM,
    custom_modules    = mods
  )

  seurat_obj <- ModulePerturbation(
    seurat_obj,
    mod               = "red",
    perturb_dir       = 0.5,
    perturb_mode      = "multiplicative",
    perturbation_name = "red_dn",
    graph             = "RNA_nn",
    n_hubs            = 5,
    n_iters           = 1,
    delta_scale       = 0.05,
    row_normalize     = FALSE,
    use_velocyto      = FALSE,
    custom_network    = TOM,
    custom_modules    = mods
  )
}

# helper: Pearson correlation between the log-normalized expression deltas
# of UP and DOWN across all module gene-cell entries. This is the direct
# input to PerturbationTransitions (the data layer delta), so a negative
# value confirms the two perturbations push expression in opposite directions
# in the space used for transition probability computation.
log_delta_correlation <- function(seurat_obj, name_up, name_dn, features,
                                  assay_obs = "RNA") {
  obs_log <- Seurat::GetAssayData(seurat_obj, assay = assay_obs,  layer = "data")[features, ]
  up_log  <- Seurat::GetAssayData(seurat_obj, assay = name_up,    layer = "data")[features, ]
  dn_log  <- Seurat::GetAssayData(seurat_obj, assay = name_dn,    layer = "data")[features, ]
  delta_up <- as.numeric(as.matrix(up_log - obs_log))
  delta_dn <- as.numeric(as.matrix(dn_log - obs_log))
  cor(delta_up, delta_dn)
}

# ---------------------------------------------------------------------------
# 1. Log-normalized expression deltas for UP and DOWN are negatively correlated
#    across all module genes and cells — confirming opposite directions in the
#    space that PerturbationTransitions uses to compute transition probabilities
# ---------------------------------------------------------------------------

test_that("log-normalized expression deltas for multiplicative UP and DOWN are negatively correlated", {
  skip_if_no_data()

  red_genes <- subset(mods, module == "red")$gene_name
  log_cor   <- log_delta_correlation(seurat_obj, "red_up", "red_dn", red_genes)
  cat(sprintf("\n  Log-delta correlation (multiplicative UP vs DOWN, all module genes): %.4f\n", log_cor))

  expect_lt(log_cor, 0,
    label = "log-normalized deltas for UP and DOWN are negatively correlated (opposite directions)")
})

# ---------------------------------------------------------------------------
# 2. Hub gene expression moves in opposite directions for UP vs DOWN
#    (sanity check that the primary perturbation is correctly applied)
# ---------------------------------------------------------------------------

test_that("hub genes increase in UP and decrease in DOWN relative to baseline", {
  skip_if_no_data()

  hub_df    <- hdWGCNA::GetHubGenes(seurat_obj, n = 5, wgcna_name = "simulation")
  hub_genes <- subset(hub_df, module == "red")$gene_name

  exp_obs <- Seurat::GetAssayData(seurat_obj, assay = "RNA",    layer = "counts")
  exp_up  <- Seurat::GetAssayData(seurat_obj, assay = "red_up", layer = "counts")
  exp_dn  <- Seurat::GetAssayData(seurat_obj, assay = "red_dn", layer = "counts")

  mean_obs <- mean(as.numeric(exp_obs[hub_genes, ]))
  mean_up  <- mean(as.numeric(exp_up[hub_genes, ]))
  mean_dn  <- mean(as.numeric(exp_dn[hub_genes, ]))

  expect_gt(mean_up, mean_obs,
    label = "multiplicative UP: hub gene mean expression is higher than baseline")
  expect_lt(mean_dn, mean_obs,
    label = "multiplicative DOWN: hub gene mean expression is lower than baseline")
})

# ---------------------------------------------------------------------------
# 3. Hub gene deltas for UP and DOWN are negatively correlated (cell-level)
#    — confirms the multiplicative mode preserves the anticorrelation structure
#    that was lost in the ZINB additive model
# ---------------------------------------------------------------------------

test_that("hub gene count deltas for UP and DOWN are negatively correlated across cells", {
  skip_if_no_data()

  hub_df    <- hdWGCNA::GetHubGenes(seurat_obj, n = 5, wgcna_name = "simulation")
  hub_genes <- subset(hub_df, module == "red")$gene_name

  exp_obs <- Seurat::GetAssayData(seurat_obj, assay = "RNA",    layer = "counts")
  exp_up  <- Seurat::GetAssayData(seurat_obj, assay = "red_up", layer = "counts")
  exp_dn  <- Seurat::GetAssayData(seurat_obj, assay = "red_dn", layer = "counts")

  delta_up <- as.numeric(as.matrix(exp_up[hub_genes, ] - exp_obs[hub_genes, ]))
  delta_dn <- as.numeric(as.matrix(exp_dn[hub_genes, ] - exp_obs[hub_genes, ]))

  hub_delta_cor <- cor(delta_up, delta_dn)
  cat(sprintf("\n  Hub delta correlation (UP vs DOWN): %.4f\n", hub_delta_cor))

  expect_lt(hub_delta_cor, 0,
    label = "hub gene deltas for UP and DOWN are negatively correlated (cell-specific structure preserved)")
})

# ---------------------------------------------------------------------------
# 4. Transition probability matrices for UP and DOWN are weakly correlated
#    (the softmax normalization in SparseColDeltaCor compresses cell-to-cell
#    differences, so the TP matrices won't be perfectly anticorrelated, but
#    with log-space propagation the correlation should be substantially lower
#    than the ~0.9 seen with count-space propagation)
# ---------------------------------------------------------------------------

test_that("TP matrices for UP and DOWN have low positive correlation (< 0.5)", {
  skip_if_no_data()

  tp_up <- Matrix::Matrix(Graphs(seurat_obj, "red_up_tp"))
  tp_dn <- Matrix::Matrix(Graphs(seurat_obj, "red_dn_tp"))

  cell_graph <- Matrix::Matrix(Graphs(seurat_obj, "RNA_nn"))
  diag(cell_graph) <- 1
  edges <- summary(methods::as(cell_graph, "TsparseMatrix"))

  vals_up <- tp_up[cbind(edges$i, edges$j)]
  vals_dn <- tp_dn[cbind(edges$i, edges$j)]
  tp_cor  <- cor(vals_up, vals_dn, use = "complete.obs")
  cat(sprintf("\n  TP correlation (UP vs DOWN): %.4f\n", tp_cor))

  expect_lt(tp_cor, 0.5,
    label = "TP matrices for UP and DOWN are not highly correlated (log-space propagation reduces similarity)")
})

# ---------------------------------------------------------------------------
# 5. Vector field cosine similarity: majority of cells have vectors pointing
#    in opposite directions between UP and DOWN perturbations.
#    Hub preservation keeps hub genes fixed at their primary perturbation
#    level throughout propagation, which strongly enforces the antisymmetry
#    between UP and DOWN. The floor (result >= 0) introduces residual
#    asymmetry for zero-expressing non-hub cells but is secondary.
# ---------------------------------------------------------------------------

test_that("vector fields for UP and DOWN point in opposite directions for most cells", {
  skip_if_no_data()

  cos_sim_rowwise <- function(a, b) {
    norm_a <- sqrt(rowSums(a^2))
    norm_b <- sqrt(rowSums(b^2))
    ok     <- norm_a > 0 & norm_b > 0
    result <- rep(NA_real_, nrow(a))
    result[ok] <- rowSums((a[ok, ] / norm_a[ok]) * (b[ok, ] / norm_b[ok]))
    result
  }

  vf_up <- PerturbationVectors(
    seurat_obj, "red_up",
    reduction    = "pca",
    arrow_scale  = 2,
    max_pct      = 0.5,
    use_velocyto = FALSE
  )$arsd

  vf_dn <- PerturbationVectors(
    seurat_obj, "red_dn",
    reduction    = "pca",
    arrow_scale  = 2,
    max_pct      = 0.5,
    use_velocyto = FALSE
  )$arsd

  sim        <- cos_sim_rowwise(as.matrix(vf_up), as.matrix(vf_dn))
  med_sim    <- median(sim, na.rm = TRUE)
  pct_neg    <- mean(sim < 0, na.rm = TRUE)

  cat(sprintf("\n  Vector field cosine similarity (UP vs DOWN):\n"))
  cat(sprintf("    Median: %.4f\n", med_sim))
  cat(sprintf("    Pct cells with opposite direction: %.1f%%\n", pct_neg * 100))

  expect_lt(med_sim, 0,
    label = "median vector field cosine similarity is negative (UP and DOWN point in opposite directions)")
  expect_gt(pct_neg, 0.6,
    label = "at least 60% of cells have opposite-direction vectors for UP vs DOWN")
})
