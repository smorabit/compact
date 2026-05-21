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

    branch1 <- colnames(seurat_obj)[seurat_obj$branch == "Branch 1"]
    branch3 <- colnames(seurat_obj)[seurat_obj$branch == "Branch 3"]

    # run compute functions so we have scores in metadata
    seurat_obj <- PredictAttractors(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        output_name       = "attractor_score",
        return_seurat     = TRUE
    )

    seurat_obj <- PredictPerturbationTime(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        sink_cells        = branch3,
        output_name       = "perturbation_pseudotime",
        max_iter          = 200,
        return_seurat     = TRUE
    )

    seurat_obj <- PredictCommitment(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        source_cells      = branch1,
        sink_cells        = branch3,
        output_name       = "commitment_score",
        max_iter          = 200,
        return_seurat     = TRUE
    )
}

# ===========================================================================
# PlotMarkovEmbedding
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. Basic usage — returns a ggplot object
# ---------------------------------------------------------------------------

test_that("PlotMarkovEmbedding returns a ggplot for basic usage", {
    skip_if_no_data()

    p <- PlotMarkovEmbedding(seurat_obj, feature = "attractor_score", reduction = "pca")

    expect_true(inherits(p, "ggplot"),
                label = "returns a ggplot object")
})

# ---------------------------------------------------------------------------
# 2. source_cells and sink_cells provided — still returns a ggplot
# ---------------------------------------------------------------------------

test_that("PlotMarkovEmbedding returns a ggplot when source and sink cells are provided", {
    skip_if_no_data()

    p <- PlotMarkovEmbedding(
        seurat_obj,
        feature      = "commitment_score",
        reduction    = "pca",
        sink_cells   = branch3,
        source_cells = branch1
    )

    expect_true(inherits(p, "ggplot"),
                label = "returns a ggplot with source/sink overlays")
})

# ---------------------------------------------------------------------------
# 3. unconverged_cells provided — still returns a ggplot
# ---------------------------------------------------------------------------

test_that("PlotMarkovEmbedding returns a ggplot when unconverged_cells is provided", {
    skip_if_no_data()

    max_iter_val <- 200
    unconverged  <- colnames(seurat_obj)[
        seurat_obj@meta.data[["perturbation_pseudotime"]] >= max_iter_val
    ]

    p <- PlotMarkovEmbedding(
        seurat_obj,
        feature           = "perturbation_pseudotime",
        reduction         = "pca",
        unconverged_cells = unconverged
    )

    expect_true(inherits(p, "ggplot"),
                label = "returns a ggplot with unconverged cells greyed out")
})

# ---------------------------------------------------------------------------
# 4. Missing feature → informative error
# ---------------------------------------------------------------------------

test_that("PlotMarkovEmbedding stops when feature is not in @meta.data", {
    skip_if_no_data()

    expect_error(
        PlotMarkovEmbedding(seurat_obj, feature = "nonexistent_score", reduction = "pca"),
        regexp = "not found in seurat_obj@meta.data"
    )
})

# ---------------------------------------------------------------------------
# 5. Missing reduction → informative error
# ---------------------------------------------------------------------------

test_that("PlotMarkovEmbedding stops when reduction is not in @reductions", {
    skip_if_no_data()

    expect_error(
        PlotMarkovEmbedding(
            seurat_obj,
            feature   = "attractor_score",
            reduction = "nonexistent_reduction"
        ),
        regexp = "not found in seurat_obj@reductions"
    )
})

# ===========================================================================
# PlotMarkovScatter
# ===========================================================================

# ---------------------------------------------------------------------------
# 6. Basic usage (no color.by) — returns a ggplot
# ---------------------------------------------------------------------------

test_that("PlotMarkovScatter returns a ggplot for basic usage without color.by", {
    skip_if_no_data()

    p <- PlotMarkovScatter(
        seurat_obj,
        x_feature = "attractor_score",
        y_feature = "perturbation_pseudotime"
    )

    expect_true(inherits(p, "ggplot"),
                label = "returns a ggplot with no color.by")
})

# ---------------------------------------------------------------------------
# 7. Categorical color.by — returns a ggplot
# ---------------------------------------------------------------------------

test_that("PlotMarkovScatter returns a ggplot with a categorical color.by", {
    skip_if_no_data()

    p <- PlotMarkovScatter(
        seurat_obj,
        x_feature = "attractor_score",
        y_feature = "perturbation_pseudotime",
        color.by  = "branch"
    )

    expect_true(inherits(p, "ggplot"),
                label = "returns a ggplot with categorical color.by")
})

# ---------------------------------------------------------------------------
# 8. Numeric color.by — returns a ggplot
# ---------------------------------------------------------------------------

test_that("PlotMarkovScatter returns a ggplot with a numeric color.by", {
    skip_if_no_data()

    p <- PlotMarkovScatter(
        seurat_obj,
        x_feature = "attractor_score",
        y_feature = "perturbation_pseudotime",
        color.by  = "commitment_score"
    )

    expect_true(inherits(p, "ggplot"),
                label = "returns a ggplot with numeric color.by")
})

# ---------------------------------------------------------------------------
# 9. add_smooth = TRUE — returns a ggplot
# ---------------------------------------------------------------------------

test_that("PlotMarkovScatter returns a ggplot with add_smooth = TRUE", {
    skip_if_no_data()

    p <- PlotMarkovScatter(
        seurat_obj,
        x_feature  = "attractor_score",
        y_feature  = "perturbation_pseudotime",
        color.by   = "branch",
        add_smooth = TRUE
    )

    expect_true(inherits(p, "ggplot"),
                label = "returns a ggplot with smooth line")
})

# ---------------------------------------------------------------------------
# 10. hline and vline provided — returns a ggplot
# ---------------------------------------------------------------------------

test_that("PlotMarkovScatter returns a ggplot when hline and vline are provided", {
    skip_if_no_data()

    p <- PlotMarkovScatter(
        seurat_obj,
        x_feature = "attractor_score",
        y_feature = "perturbation_pseudotime",
        hline     = 0.5,
        vline     = 0.5
    )

    expect_true(inherits(p, "ggplot"),
                label = "returns a ggplot with hline and vline")
})

# ---------------------------------------------------------------------------
# 11. Missing x_feature → informative error
# ---------------------------------------------------------------------------

test_that("PlotMarkovScatter stops when x_feature is not in @meta.data", {
    skip_if_no_data()

    expect_error(
        PlotMarkovScatter(
            seurat_obj,
            x_feature = "nonexistent_feature",
            y_feature = "attractor_score"
        ),
        regexp = "not found in seurat_obj@meta.data"
    )
})

# ---------------------------------------------------------------------------
# 12. Missing y_feature → informative error
# ---------------------------------------------------------------------------

test_that("PlotMarkovScatter stops when y_feature is not in @meta.data", {
    skip_if_no_data()

    expect_error(
        PlotMarkovScatter(
            seurat_obj,
            x_feature = "attractor_score",
            y_feature = "nonexistent_feature"
        ),
        regexp = "not found in seurat_obj@meta.data"
    )
})
