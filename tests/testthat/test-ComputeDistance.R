library(testthat)
library(compact)
library(Seurat)

# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------

data_path <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/simulation_branch.rds"

skip_if_no_data <- function() {
  skip_if(!file.exists(data_path), "Test dataset not found")
}

if (file.exists(data_path)) {
  seurat_obj <- readRDS(data_path)
  groups     <- sort(unique(as.character(seurat_obj@meta.data[["branch"]])))
  n_groups   <- length(groups)  # 3
}

# ---------------------------------------------------------------------------
# 1. Euclidean: returns a data.frame
# ---------------------------------------------------------------------------

test_that("euclidean method returns a data.frame", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "euclidean",
                            verbose = FALSE)

  expect_true(is.data.frame(result), label = "euclidean result is a data.frame")
})

# ---------------------------------------------------------------------------
# 2. Euclidean: correct dimensions
# ---------------------------------------------------------------------------

test_that("euclidean result has n_groups x n_groups dimensions", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "euclidean",
                            verbose = FALSE)

  expect_equal(nrow(result), n_groups, label = "euclidean: correct number of rows")
  expect_equal(ncol(result), n_groups, label = "euclidean: correct number of columns")
})

# ---------------------------------------------------------------------------
# 3. Euclidean: row and column names match group names
# ---------------------------------------------------------------------------

test_that("euclidean result row and column names match group labels", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "euclidean",
                            verbose = FALSE)

  expect_equal(sort(rownames(result)), sort(groups),
    label = "euclidean: rownames match group labels")
  expect_equal(sort(colnames(result)), sort(groups),
    label = "euclidean: colnames match group labels")
})

# ---------------------------------------------------------------------------
# 4. Euclidean: symmetric matrix
# ---------------------------------------------------------------------------

test_that("euclidean result is a symmetric matrix", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "euclidean",
                            verbose = FALSE)

  mat <- as.matrix(result)
  expect_equal(mat, t(mat), label = "euclidean: matrix is symmetric")
})

# ---------------------------------------------------------------------------
# 5. Euclidean: diagonal is zero
# ---------------------------------------------------------------------------

test_that("euclidean result has zero diagonal (distance from group to itself)", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "euclidean",
                            verbose = FALSE)

  mat <- as.matrix(result)
  expect_true(all(diag(mat) == 0), label = "euclidean: diagonal entries are zero")
})

# ---------------------------------------------------------------------------
# 6. Euclidean: all off-diagonal values are non-negative and finite
# ---------------------------------------------------------------------------

test_that("euclidean off-diagonal values are non-negative and finite", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "euclidean",
                            verbose = FALSE)

  mat <- as.matrix(result)
  off_diag <- mat[row(mat) != col(mat)]

  expect_true(all(is.finite(off_diag)), label = "euclidean: all off-diagonal values are finite")
  expect_true(all(off_diag >= 0),       label = "euclidean: all off-diagonal values are non-negative")
})

# ---------------------------------------------------------------------------
# 7. Spearman: returns a data.frame with correct dimensions and names
# ---------------------------------------------------------------------------

test_that("spearman method returns a correctly shaped data.frame", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "spearman",
                            verbose = FALSE)

  expect_true(is.data.frame(result), label = "spearman result is a data.frame")
  expect_equal(nrow(result), n_groups,  label = "spearman: correct row count")
  expect_equal(ncol(result), n_groups,  label = "spearman: correct column count")
  expect_equal(sort(rownames(result)), sort(groups), label = "spearman: rownames match groups")
  expect_equal(sort(colnames(result)), sort(groups), label = "spearman: colnames match groups")
})

# ---------------------------------------------------------------------------
# 8. Spearman: symmetric with zero diagonal
# ---------------------------------------------------------------------------

test_that("spearman result is symmetric with zero diagonal", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "spearman",
                            verbose = FALSE)

  mat <- as.matrix(result)
  expect_equal(mat, t(mat),            label = "spearman: matrix is symmetric")
  expect_true(all(diag(mat) == 0),     label = "spearman: diagonal entries are zero")
})

# ---------------------------------------------------------------------------
# 9. Spearman: off-diagonal values in [0, 2] (distance = 1 - rho)
# ---------------------------------------------------------------------------

test_that("spearman off-diagonal values are in [0, 2]", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "spearman",
                            verbose = FALSE)

  mat      <- as.matrix(result)
  off_diag <- mat[row(mat) != col(mat)]

  expect_true(all(off_diag >= 0), label = "spearman: distances are non-negative (rho <= 1)")
  expect_true(all(off_diag <= 2), label = "spearman: distances at most 2 (rho >= -1)")
})

# ---------------------------------------------------------------------------
# 10. Energy distance: returns a data.frame with correct dimensions and names
# ---------------------------------------------------------------------------

test_that("edist method returns a correctly shaped data.frame", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "edist",
                            verbose = FALSE)

  expect_true(is.data.frame(result), label = "edist result is a data.frame")
  expect_equal(nrow(result), n_groups,  label = "edist: correct row count")
  expect_equal(ncol(result), n_groups,  label = "edist: correct column count")
  expect_equal(sort(rownames(result)), sort(groups), label = "edist: rownames match groups")
  expect_equal(sort(colnames(result)), sort(groups), label = "edist: colnames match groups")
})

# ---------------------------------------------------------------------------
# 11. Energy distance: symmetric with zero diagonal
# ---------------------------------------------------------------------------

test_that("edist result is symmetric with zero diagonal", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "edist",
                            verbose = FALSE)

  mat <- as.matrix(result)
  expect_equal(mat, t(mat),        label = "edist: matrix is symmetric")
  expect_true(all(diag(mat) == 0), label = "edist: diagonal entries are zero")
})

# ---------------------------------------------------------------------------
# 12. Energy distance: all values are finite (no NA or NaN)
# ---------------------------------------------------------------------------

test_that("edist result contains no NA or NaN values", {
  skip_if_no_data()

  result <- ComputeDistance(seurat_obj, groupby = "branch",
                            reduction = "pca", method = "edist",
                            verbose = FALSE)

  expect_true(all(is.finite(as.matrix(result))),
    label = "edist: all entries are finite")
})

# ---------------------------------------------------------------------------
# 13. Invalid method → match.arg error
# ---------------------------------------------------------------------------

test_that("invalid method argument stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ComputeDistance(seurat_obj, groupby = "branch",
                    reduction = "pca", method = "cosine"),
    label = "invalid method triggers match.arg error"
  )
})

# ---------------------------------------------------------------------------
# 14. Invalid reduction → informative stop
# ---------------------------------------------------------------------------

test_that("invalid reduction name stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ComputeDistance(seurat_obj, groupby = "branch",
                    reduction = "nonexistent_reduction", method = "euclidean"),
    label = "invalid reduction triggers stop"
  )
  expect_error(
    ComputeDistance(seurat_obj, groupby = "branch",
                    reduction = "nonexistent_reduction", method = "spearman"),
    label = "invalid reduction triggers stop (spearman)"
  )
  expect_error(
    ComputeDistance(seurat_obj, groupby = "branch",
                    reduction = "nonexistent_reduction", method = "edist"),
    label = "invalid reduction triggers stop (edist)"
  )
})

# ---------------------------------------------------------------------------
# 15. Invalid groupby column → informative stop
# ---------------------------------------------------------------------------

test_that("invalid groupby column stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ComputeDistance(seurat_obj, groupby = "nonexistent_column",
                    reduction = "pca", method = "euclidean"),
    label = "invalid groupby triggers stop (euclidean)"
  )
  expect_error(
    ComputeDistance(seurat_obj, groupby = "nonexistent_column",
                    reduction = "pca", method = "edist"),
    label = "invalid groupby triggers stop (edist)"
  )
})

# ---------------------------------------------------------------------------
# 16. Non-Seurat input → error
# ---------------------------------------------------------------------------

test_that("non-Seurat first argument stops with an error", {
  skip_if_no_data()

  expect_error(
    ComputeDistance(list(), groupby = "branch",
                    reduction = "pca", method = "euclidean"),
    label = "list input triggers error (euclidean)"
  )
  expect_error(
    ComputeDistance(list(), groupby = "branch",
                    reduction = "pca", method = "edist"),
    label = "list input triggers error (edist)"
  )
})

# ---------------------------------------------------------------------------
# 17. verbose = FALSE runs without error for all methods
# ---------------------------------------------------------------------------

test_that("verbose=FALSE runs without error for all methods", {
  skip_if_no_data()

  expect_no_error(
    ComputeDistance(seurat_obj, groupby = "branch",
                    reduction = "pca", method = "euclidean", verbose = FALSE)
  )
  expect_no_error(
    ComputeDistance(seurat_obj, groupby = "branch",
                    reduction = "pca", method = "spearman", verbose = FALSE)
  )
  expect_no_error(
    ComputeDistance(seurat_obj, groupby = "branch",
                    reduction = "pca", method = "edist", verbose = FALSE)
  )
})
