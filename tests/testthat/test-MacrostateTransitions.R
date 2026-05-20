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

    # run a perturbation so we have a TP graph to work with
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

    # also run a DOWN perturbation for the directional cross-check test
    seurat_obj <- ModulePerturbation(
        seurat_obj,
        mod               = "red",
        perturb_dir       = -1,
        perturbation_name = "red_down",
        graph             = "RNA_nn",
        n_hubs            = 5,
        n_iters           = 3,
        use_velocyto      = FALSE,
        custom_network    = TOM,
        custom_modules    = mods
    )

    n_cells   <- ncol(seurat_obj)
    n_groups  <- length(unique(seurat_obj@meta.data[["branch"]]))
}

# ===========================================================================
# MacrostateTransitions — computation
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. Returns a Seurat object with result stored in @misc
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions returns a Seurat object and stores result in @misc", {
    skip_if_no_data()

    result_obj <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = TRUE,
        verbose           = FALSE
    )

    expect_s4_class(result_obj, "Seurat")
    expect_true("MacrostateTransitions" %in% names(result_obj@misc),
                label = "MacrostateTransitions entry created in @misc")
    expect_true("red_up" %in% names(result_obj@misc$MacrostateTransitions),
                label = "perturbation_name key exists in @misc$MacrostateTransitions")
})

# ---------------------------------------------------------------------------
# 2. Q is row-stochastic (rows sum to 1)
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions: Q is row-stochastic", {
    skip_if_no_data()

    result <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = FALSE,
        verbose           = FALSE
    )

    row_sums <- rowSums(result$Q)
    expect_true(all(abs(row_sums - 1) < 1e-10),
                label = "all rows of Q sum to 1")
})

# ---------------------------------------------------------------------------
# 3. Stability index equals the diagonal of Q
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions: stability index matches diagonal of Q", {
    skip_if_no_data()

    result <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = FALSE,
        verbose           = FALSE
    )

    expect_equal(result$stability, diag(result$Q), tolerance = 1e-12,
                 label = "stability == diag(Q)")
})

# ---------------------------------------------------------------------------
# 4. Q dimensions match the number of groups in the metadata column
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions: Q dimensions equal the number of groups", {
    skip_if_no_data()

    result <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = FALSE,
        verbose           = FALSE
    )

    expect_equal(nrow(result$Q), n_groups,
                 label = "Q has correct number of rows")
    expect_equal(ncol(result$Q), n_groups,
                 label = "Q has correct number of columns")
    expect_equal(nrow(result$Q), ncol(result$Q),
                 label = "Q is square")
})

# ---------------------------------------------------------------------------
# 5. group_sizes match actual cell counts per group
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions: group_sizes match table() of metadata column", {
    skip_if_no_data()

    result <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = FALSE,
        verbose           = FALSE
    )

    expected <- table(seurat_obj@meta.data[["branch"]])
    # result$group_sizes should match expected counts for the same group names
    for (g in names(result$group_sizes)) {
        expect_equal(
            result$group_sizes[[g]],
            as.integer(expected[[g]]),
            label = paste0("group size for '", g, "' matches table()")
        )
    }
})

# ---------------------------------------------------------------------------
# 6. min_group_size filter excludes small groups
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions: min_group_size filter excludes small groups", {
    skip_if_no_data()

    # use a large threshold that will exclude some groups
    large_threshold <- floor(n_cells / n_groups) * 2  # bigger than any single group if unbalanced

    # guard: only run if the threshold would actually filter something
    group_tab <- table(seurat_obj@meta.data[["branch"]])
    small_groups_exist <- any(group_tab < large_threshold)
    skip_if(!small_groups_exist, "No small groups to filter with chosen threshold")

    expect_warning(
        result <- MacrostateTransitions(
            seurat_obj,
            perturbation_name = "red_up",
            graph             = "RNA_nn",
            group.by          = "branch",
            min_group_size    = large_threshold,
            store_result      = FALSE,
            verbose           = FALSE
        )
    )

    # all retained groups must have >= threshold cells
    expect_true(all(result$group_sizes >= large_threshold),
                label = "all retained groups have >= min_group_size cells")
})

# ---------------------------------------------------------------------------
# 7. Stability values are in [0, 1]
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions: stability index values are in [0, 1]", {
    skip_if_no_data()

    result <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = FALSE,
        verbose           = FALSE
    )

    expect_true(all(result$stability >= -1e-10 & result$stability <= 1 + 1e-10),
                label = "all stability values in [0, 1]")
})

# ---------------------------------------------------------------------------
# 8. Missing TP graph → informative error
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions stops when the TP graph is absent", {
    skip_if_no_data()

    expect_error(
        MacrostateTransitions(
            seurat_obj,
            perturbation_name = "nonexistent_perturbation",
            graph             = "RNA_nn",
            group.by          = "branch"
        ),
        regexp = "not found in seurat_obj@graphs"
    )
})

# ---------------------------------------------------------------------------
# 9. Missing group.by column → informative error
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions stops when group.by column is absent", {
    skip_if_no_data()

    expect_error(
        MacrostateTransitions(
            seurat_obj,
            perturbation_name = "red_up",
            graph             = "RNA_nn",
            group.by          = "not_a_real_column"
        ),
        regexp = "not found in seurat_obj@meta.data"
    )
})

# ---------------------------------------------------------------------------
# 10. Missing KNN graph → informative error
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions stops when the KNN graph is absent", {
    skip_if_no_data()

    expect_error(
        MacrostateTransitions(
            seurat_obj,
            perturbation_name = "red_up",
            graph             = "nonexistent_knn_graph",
            group.by          = "branch"
        ),
        regexp = "not found in seurat_obj@graphs"
    )
})

# ---------------------------------------------------------------------------
# 11. Single-group edge case: Q is 1×1 with Q[1,1] = 1
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions single-group edge case: Q is 1x1 with value 1", {
    skip_if_no_data()

    seurat_obj_tmp <- seurat_obj
    seurat_obj_tmp$all_one <- "group_A"

    result <- MacrostateTransitions(
        seurat_obj_tmp,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "all_one",
        store_result      = FALSE,
        verbose           = FALSE
    )

    expect_equal(dim(result$Q), c(1L, 1L),
                 label = "Q is 1×1 for single-group input")
    expect_equal(result$Q[1, 1], 1, tolerance = 1e-10,
                 label = "Q[1,1] = 1 for single-group input")
    expect_equal(result$stability[[1]], 1, tolerance = 1e-10,
                 label = "stability = 1 for single group")
})

# ---------------------------------------------------------------------------
# 12. store_result = FALSE returns the list directly (not a Seurat object)
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions store_result = FALSE returns the result list", {
    skip_if_no_data()

    result <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = FALSE,
        verbose           = FALSE
    )

    expect_type(result, "list")
    expect_true(all(c("Q", "stability", "group_sizes", "group.by", "perturbation_name")
                    %in% names(result)),
                label = "result list contains all expected fields")
})

# ---------------------------------------------------------------------------
# 13. result_name parameter controls the key in @misc
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions stores result under custom result_name", {
    skip_if_no_data()

    result_obj <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = TRUE,
        result_name       = "my_custom_key",
        verbose           = FALSE
    )

    expect_true("my_custom_key" %in% names(result_obj@misc$MacrostateTransitions),
                label = "custom result_name used as @misc key")
})

# ===========================================================================
# PlotMacrostateTransitions — visualization
# ===========================================================================

# ---------------------------------------------------------------------------
# 14. PlotMacrostateTransitions returns a ggplot object
# ---------------------------------------------------------------------------

test_that("PlotMacrostateTransitions returns a ggplot object", {
    skip_if_no_data()

    result_obj <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = TRUE,
        verbose           = FALSE
    )

    p <- PlotMacrostateTransitions(
        result_obj,
        perturbation_name = "red_up"
    )

    expect_true(inherits(p, "ggplot"),
                label = "PlotMacrostateTransitions returns a ggplot object")
})

# ---------------------------------------------------------------------------
# 15. PlotMacrostateTransitions errors when no MacrostateTransitions in @misc
# ---------------------------------------------------------------------------

test_that("PlotMacrostateTransitions stops when @misc$MacrostateTransitions is absent", {
    skip_if_no_data()

    # use a fresh seurat_obj that hasn't had MacrostateTransitions run
    seurat_obj_fresh <- seurat_obj
    seurat_obj_fresh@misc$MacrostateTransitions <- NULL

    expect_error(
        PlotMacrostateTransitions(
            seurat_obj_fresh,
            perturbation_name = "red_up"
        ),
        regexp = "No MacrostateTransitions results found"
    )
})

# ---------------------------------------------------------------------------
# 16. PlotMacrostateTransitions errors when result_name key is absent
# ---------------------------------------------------------------------------

test_that("PlotMacrostateTransitions stops when result_name key is missing", {
    skip_if_no_data()

    result_obj <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = TRUE,
        verbose           = FALSE
    )

    expect_error(
        PlotMacrostateTransitions(
            result_obj,
            perturbation_name = "red_up",
            result_name       = "nonexistent_key"
        ),
        regexp = "not found in seurat_obj@misc\\$MacrostateTransitions"
    )
})

# ===========================================================================
# Mathematical consistency cross-checks
# ===========================================================================

# ---------------------------------------------------------------------------
# 17. Q all-cell self-transition check: diagonal entries are non-negative
#     and do not exceed 1 even for small datasets
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions: all Q entries are in [0, 1]", {
    skip_if_no_data()

    result <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = FALSE,
        verbose           = FALSE
    )

    expect_true(all(result$Q >= -1e-10),
                label = "all Q entries are non-negative")
    expect_true(all(result$Q <= 1 + 1e-10),
                label = "all Q entries are at most 1")
})

# ---------------------------------------------------------------------------
# 18. UP vs DOWN perturbation: stability should differ
#     (the two perturbations should produce different Q matrices)
# ---------------------------------------------------------------------------

test_that("MacrostateTransitions: UP and DOWN perturbations yield different Q matrices", {
    skip_if_no_data()

    result_up <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_up",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = FALSE,
        verbose           = FALSE
    )

    result_down <- MacrostateTransitions(
        seurat_obj,
        perturbation_name = "red_down",
        graph             = "RNA_nn",
        group.by          = "branch",
        store_result      = FALSE,
        verbose           = FALSE
    )

    # the two Q matrices should not be identical
    expect_false(
        isTRUE(all.equal(result_up$Q, result_down$Q, tolerance = 1e-6)),
        label = "Q_up and Q_down are not identical matrices"
    )
})
