library(testthat)
library(compact)
library(Seurat)
library(Matrix)

# ---------------------------------------------------------------------------
# Test fixtures — loaded once at file scope
# ---------------------------------------------------------------------------

data_path <- "/home/groups/singlecell/smorabito/analysis/COMPACT/data/simulation_branch.rds"

skip_if_no_data <- function() {
  skip_if(!file.exists(data_path), "Test dataset not found")
}

if (file.exists(data_path)) {
  seurat_obj <- readRDS(data_path)
  exp        <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

  # pick a gene known to be expressed (non-zero in most cells)
  gene_counts   <- Matrix::rowSums(exp)
  expressed_gene <- names(which(gene_counts > 0))[1]

  # build a zero-expression gene by temporarily zeroing out one row;
  # do this as a modified Seurat object so ModelZINB can use FetchData
  seurat_zero           <- seurat_obj
  exp_zero              <- exp
  exp_zero[expressed_gene, ] <- 0
  seurat_zero[["RNA"]]  <- Seurat::SetAssayData(
    seurat_zero[["RNA"]], layer = "counts", new.data = exp_zero
  )
}

###############################################################################
# ModelZINB tests
###############################################################################

# ---------------------------------------------------------------------------
# 1. Returns a zeroinfl model object for a normally expressed gene
# ---------------------------------------------------------------------------

test_that("ModelZINB returns a zeroinfl model for an expressed gene", {
  skip_if_no_data()

  model <- ModelZINB(seurat_obj, feature = expressed_gene)
  expect_s3_class(model, "zeroinfl")
})

# ---------------------------------------------------------------------------
# 2. Model has theta and coefficients slots (used by SampleZINB)
# ---------------------------------------------------------------------------

test_that("ModelZINB result has theta and coefficients slots", {
  skip_if_no_data()

  model <- ModelZINB(seurat_obj, feature = expressed_gene)
  expect_true(!is.null(model$theta),         label = "theta slot present")
  expect_true(!is.null(model$coefficients),  label = "coefficients slot present")
})

# ---------------------------------------------------------------------------
# 3. cells_use subsets the data (model$n matches cells + added zero)
# ---------------------------------------------------------------------------

test_that("cells_use limits the cells used for fitting", {
  skip_if_no_data()

  cells_sub <- colnames(seurat_obj)[1:50]
  model     <- ModelZINB(seurat_obj, feature = expressed_gene, cells_use = cells_sub)

  # model$n = 50 cells + 1 artificial zero appended by add_zero = TRUE
  expect_equal(model$n, 51L, label = "model$n equals length(cells_use) + 1")
})

# ---------------------------------------------------------------------------
# 4. add_zero = FALSE: model$n equals exactly length(cells_use)
# ---------------------------------------------------------------------------

test_that("add_zero=FALSE gives model$n equal to length(cells_use)", {
  skip_if_no_data()

  cells_sub <- colnames(seurat_obj)[1:50]
  model     <- ModelZINB(seurat_obj, feature = expressed_gene,
                         cells_use = cells_sub, add_zero = FALSE)

  expect_equal(model$n, 50L, label = "model$n equals length(cells_use) exactly")
})

# ---------------------------------------------------------------------------
# 5. All-zero expression → clear error (not an opaque convergence failure)
# ---------------------------------------------------------------------------

test_that("all-zero expression stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ModelZINB(seurat_zero, feature = expressed_gene),
    regexp = "all expression values.*are zero",
    label  = "all-zero gene gives clear error before zeroinfl is called"
  )
})

# ---------------------------------------------------------------------------
# 6. Feature not in Seurat object → informative error
# ---------------------------------------------------------------------------

test_that("feature not in seurat_obj stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ModelZINB(seurat_obj, feature = "FAKE_GENE_XYZ"),
    regexp = "not found in rownames",
    label  = "bad feature name gives informative error"
  )
})

# ---------------------------------------------------------------------------
# 7. Invalid cells_use barcodes → informative error
# ---------------------------------------------------------------------------

test_that("invalid cells_use barcodes stop with an informative error", {
  skip_if_no_data()

  bad_cells <- c(colnames(seurat_obj)[1:10], "FAKE_BARCODE_1", "FAKE_BARCODE_2")

  expect_error(
    ModelZINB(seurat_obj, feature = expressed_gene, cells_use = bad_cells),
    regexp = "barcodes not found in colnames",
    label  = "invalid barcodes in cells_use give informative error"
  )
})

# ---------------------------------------------------------------------------
# 8. Invalid slot → informative error
# ---------------------------------------------------------------------------

test_that("invalid slot stops with an informative error", {
  skip_if_no_data()

  expect_error(
    ModelZINB(seurat_obj, feature = expressed_gene, slot = "bad_slot"),
    regexp = "Invalid slot",
    label  = "bad slot gives informative error"
  )
})

###############################################################################
# SampleZINB tests
###############################################################################

# ---------------------------------------------------------------------------
# 9. SampleZINB returns a numeric vector of the requested length
# ---------------------------------------------------------------------------

test_that("SampleZINB returns a numeric vector of length ncells", {
  skip_if_no_data()
  set.seed(42)

  model  <- ModelZINB(seurat_obj, feature = expressed_gene)
  yobs   <- as.numeric(exp[expressed_gene, ])
  n      <- 100L
  result <- SampleZINB(model, yobs = yobs, ncells = n)

  expect_length(result, n)
  expect_true(is.numeric(result), label = "output is numeric")
})

# ---------------------------------------------------------------------------
# 10. SampleZINB output is non-negative integer-valued counts
# ---------------------------------------------------------------------------

test_that("SampleZINB output contains non-negative integer-valued counts", {
  skip_if_no_data()
  set.seed(42)

  model  <- ModelZINB(seurat_obj, feature = expressed_gene)
  yobs   <- as.numeric(exp[expressed_gene, ])
  result <- SampleZINB(model, yobs = yobs, ncells = 200L)

  expect_true(min(result) >= 0, label = "no negative simulated counts")
  expect_true(all(result == floor(result)), label = "all values are integer-valued")
})

# ---------------------------------------------------------------------------
# 11. ncells passed explicitly overrides model$n (off-by-one guard)
# ---------------------------------------------------------------------------

test_that("explicit ncells gives exactly that many simulated values", {
  skip_if_no_data()
  set.seed(1)

  model  <- ModelZINB(seurat_obj, feature = expressed_gene)
  yobs   <- as.numeric(exp[expressed_gene, ])
  n_want <- ncol(seurat_obj)
  result <- SampleZINB(model, yobs = yobs, ncells = n_want)

  expect_length(result, n_want,
    label = "simulated vector length equals the explicitly requested ncells")
})

# ---------------------------------------------------------------------------
# 12. All-zero yobs → returns a zero vector without calling rzinegbin
# ---------------------------------------------------------------------------

test_that("all-zero yobs returns a zero vector of length ncells", {
  skip_if_no_data()

  model  <- ModelZINB(seurat_obj, feature = expressed_gene)
  yobs   <- rep(0, 100)
  result <- SampleZINB(model, yobs = yobs, ncells = 50L)

  expect_length(result, 50L, label = "zero-yobs output has correct length")
  expect_true(all(result == 0), label = "all simulated values are zero when yobs is all-zero")
})

# ---------------------------------------------------------------------------
# 13. Full round-trip: ModelZINB → SampleZINB with ncells = length(cells_use)
#     produces the correct number of values (regression for the off-by-one bug)
# ---------------------------------------------------------------------------

test_that("ModelZINB → SampleZINB round-trip with ncells=length(cells_use) is not off-by-one", {
  skip_if_no_data()
  set.seed(7)

  cells_sub <- colnames(seurat_obj)[1:80]
  model     <- ModelZINB(seurat_obj, feature = expressed_gene, cells_use = cells_sub)
  yobs      <- as.numeric(exp[expressed_gene, cells_sub])
  result    <- SampleZINB(model, yobs = yobs, ncells = length(cells_sub))

  expect_length(result, length(cells_sub),
    label = "round-trip output length matches number of input cells, not model$n")
})
