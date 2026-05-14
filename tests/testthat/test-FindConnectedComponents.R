library(testthat)
library(compact)
library(Seurat)
library(Matrix)

# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------

data_path <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/simulation_branch.rds"

skip_if_no_data <- function() {
    skip_if(!file.exists(data_path), "Test dataset not found")
}

# build a minimal toy Seurat object with a known two-component graph
make_toy_seurat <- function() {
    set.seed(42)
    n_genes <- 20
    n_cells <- 10
    counts <- matrix(
        rpois(n_genes * n_cells, lambda = 5),
        nrow = n_genes, ncol = n_cells,
        dimnames = list(
            paste0("Gene", seq_len(n_genes)),
            paste0("Cell", seq_len(n_cells))
        )
    )
    obj <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)

    # build a block-diagonal SNN: cells 1-6 form one component, cells 7-10 another
    adj <- Matrix(0, nrow = n_cells, ncol = n_cells,
                  dimnames = list(colnames(obj), colnames(obj)), sparse = TRUE)
    # component 1: cells 1-6
    for (i in 1:6) for (j in 1:6) if (i != j) adj[i, j] <- 0.5
    # component 2: cells 7-10
    for (i in 7:10) for (j in 7:10) if (i != j) adj[i, j] <- 0.5

    snn <- as.Graph(adj)
    obj@graphs[["RNA_snn"]] <- snn

    # also add a second graph for auto-detect tests
    adj2 <- adj
    snn2 <- as.Graph(adj2)
    obj@graphs[["RNA_nn"]] <- snn2

    obj
}

toy <- make_toy_seurat()

# ---------------------------------------------------------------------------
# 1. Returns a Seurat object
# ---------------------------------------------------------------------------

test_that("FindConnectedComponents returns a Seurat object", {
    result <- FindConnectedComponents(toy, graph = "RNA_snn", verbose = FALSE)
    expect_s4_class(result, "Seurat")
})

# ---------------------------------------------------------------------------
# 2. meta.data column is created with the default name
# ---------------------------------------------------------------------------

test_that("default meta_data_name column is created", {
    result <- FindConnectedComponents(toy, graph = "RNA_snn", verbose = FALSE)
    expect_true("connected_component" %in% colnames(result@meta.data))
})

# ---------------------------------------------------------------------------
# 3. Custom meta_data_name is respected
# ---------------------------------------------------------------------------

test_that("custom meta_data_name is written correctly", {
    result <- FindConnectedComponents(toy, graph = "RNA_snn",
                                      meta_data_name = "my_components",
                                      verbose = FALSE)
    expect_true("my_components" %in% colnames(result@meta.data))
})

# ---------------------------------------------------------------------------
# 4. Correct number of components detected (toy graph has 2 components)
# ---------------------------------------------------------------------------

test_that("correct number of components detected in toy graph", {
    result <- FindConnectedComponents(toy, graph = "RNA_snn", verbose = FALSE)
    n_comp <- nlevels(result@meta.data[["connected_component"]])
    expect_equal(n_comp, 2)
})

# ---------------------------------------------------------------------------
# 5. Component labels are a factor
# ---------------------------------------------------------------------------

test_that("component labels are stored as a factor", {
    result <- FindConnectedComponents(toy, graph = "RNA_snn", verbose = FALSE)
    expect_true(is.factor(result@meta.data[["connected_component"]]))
})

# ---------------------------------------------------------------------------
# 6. Every cell receives a label (no NAs)
# ---------------------------------------------------------------------------

test_that("no NAs in component labels", {
    result <- FindConnectedComponents(toy, graph = "RNA_snn", verbose = FALSE)
    expect_false(any(is.na(result@meta.data[["connected_component"]])))
})

# ---------------------------------------------------------------------------
# 7. Largest component gets label "1"
# ---------------------------------------------------------------------------

test_that("largest component gets label 1", {
    result <- FindConnectedComponents(toy, graph = "RNA_snn", verbose = FALSE)
    md <- result@meta.data
    # cells 1-6 form the larger component (6 cells) and should be labelled "1"
    large_cells <- paste0("Cell", 1:6)
    small_cells <- paste0("Cell", 7:10)
    expect_true(all(md[large_cells, "connected_component"] == "1"))
    expect_true(all(md[small_cells, "connected_component"] == "2"))
})

# ---------------------------------------------------------------------------
# 8. Auto-detection picks the _snn graph
# ---------------------------------------------------------------------------

test_that("auto-detect picks the _snn graph when available", {
    # toy has both RNA_snn and RNA_nn; auto-detect should prefer RNA_snn
    result <- NULL
    expect_message(
        { result <- FindConnectedComponents(toy, verbose = TRUE) },
        regexp = "RNA_snn"
    )
    expect_s4_class(result, "Seurat")
})

# ---------------------------------------------------------------------------
# 9. Fully connected graph yields one component
# ---------------------------------------------------------------------------

test_that("fully connected graph yields one component", {
    n <- 8
    full_adj <- matrix(0.5, n, n)
    diag(full_adj) <- 0
    cells <- paste0("C", seq_len(n))
    dimnames(full_adj) <- list(cells, cells)

    counts <- matrix(rpois(20 * n, 5), 20, n,
                     dimnames = list(paste0("G", 1:20), cells))
    obj2 <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
    obj2@graphs[["RNA_snn"]] <- as.Graph(Matrix(full_adj, sparse = TRUE))

    result <- FindConnectedComponents(obj2, graph = "RNA_snn", verbose = FALSE)
    n_comp <- nlevels(result@meta.data[["connected_component"]])
    expect_equal(n_comp, 1)
})

# ---------------------------------------------------------------------------
# 10. Invalid graph name → informative error
# ---------------------------------------------------------------------------

test_that("invalid graph name gives an informative error", {
    expect_error(
        FindConnectedComponents(toy, graph = "nonexistent_graph", verbose = FALSE),
        regexp = "not found in seurat_obj@graphs"
    )
})

# ---------------------------------------------------------------------------
# 11. No graphs in seurat object → informative error
# ---------------------------------------------------------------------------

test_that("no graphs in object gives an informative error", {
    obj_empty <- toy
    obj_empty@graphs <- list()
    expect_error(
        FindConnectedComponents(obj_empty, verbose = FALSE),
        regexp = "No graphs found"
    )
})

# ---------------------------------------------------------------------------
# 12. Real data: simulation_branch dataset with varying FindNeighbors params
# ---------------------------------------------------------------------------

test_that("component count decreases as k.param increases on real data", {
    skip_if_no_data()

    seurat_obj <- readRDS(data_path)

    k_params  <- c(5, 15, 30)
    n_comps   <- integer(length(k_params))

    for (i in seq_along(k_params)) {
        obj_k <- Seurat::FindNeighbors(
            seurat_obj,
            dims      = 1:20,
            k.param   = k_params[i],
            verbose   = FALSE
        )
        obj_k    <- FindConnectedComponents(obj_k, verbose = FALSE)
        comp_col <- obj_k@meta.data[["connected_component"]]
        n_comps[i] <- nlevels(comp_col)
    }

    # higher k should produce fewer or equal components
    expect_true(n_comps[1] >= n_comps[2],
        label = "k=5 has >= components than k=15")
    expect_true(n_comps[2] >= n_comps[3],
        label = "k=15 has >= components than k=30")
})

test_that("component count increases as prune.SNN threshold increases on real data", {
    skip_if_no_data()

    seurat_obj <- readRDS(data_path)

    prune_vals <- c(0, 1/10, 1/5)
    n_comps    <- integer(length(prune_vals))

    for (i in seq_along(prune_vals)) {
        obj_p <- Seurat::FindNeighbors(
            seurat_obj,
            dims      = 1:20,
            k.param   = 20,
            prune.SNN = prune_vals[i],
            verbose   = FALSE
        )
        obj_p    <- FindConnectedComponents(obj_p, verbose = FALSE)
        comp_col <- obj_p@meta.data[["connected_component"]]
        n_comps[i] <- nlevels(comp_col)
    }

    # stronger pruning removes more edges → more or equal components
    expect_true(n_comps[1] <= n_comps[2],
        label = "prune=0 has <= components than prune=1/10")
    expect_true(n_comps[2] <= n_comps[3],
        label = "prune=1/10 has <= components than prune=1/5")
})
