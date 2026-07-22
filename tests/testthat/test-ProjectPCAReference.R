library(testthat)
library(compact)
library(Seurat)


.make_projection_fixture <- function() {
  set.seed(42)
  counts <- matrix(
    rpois(20L * 24L, lambda = 5),
    nrow = 20L,
    dimnames = list(
      paste0("gene", seq_len(20L)),
      paste0("cell", seq_len(24L))
    )
  )

  object <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA")
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  Seurat::VariableFeatures(object) <- rownames(object)
  object <- Seurat::ScaleData(
    object, features = rownames(object), scale.max = 10, verbose = FALSE
  )
  object <- Seurat::RunPCA(
    object, features = rownames(object), npcs = 5,
    approx = FALSE, verbose = FALSE
  )

  reference_data <- tryCatch(
    Seurat::GetAssayData(object, assay = "RNA", layer = "data"),
    error = function(e) Seurat::GetAssayData(object, assay = "RNA", slot = "data")
  )
  query_data <- reference_data
  query_data[1:2, 1:6] <- query_data[1:2, 1:6] + 0.25

  object[["perturbed"]] <- Seurat::CreateAssayObject(counts = counts)
  object <- tryCatch(
    Seurat::SetAssayData(
      object, assay = "perturbed", layer = "data", new.data = query_data
    ),
    error = function(e) Seurat::SetAssayData(
      object, assay = "perturbed", slot = "data", new.data = query_data
    )
  )

  list(object = object, unchanged = colnames(object)[7:24])
}


test_that("ProjectPCAReference stores a fixed-reference reduction", {
  fixture <- .make_projection_fixture()
  result <- ProjectPCAReference(
    fixture$object,
    query_assay = "perturbed",
    reference_assay = "RNA",
    reference_reduction = "pca",
    dims = 1:3,
    unchanged.cells = fixture$unchanged,
    reduction.name = "perturbed_fixed",
    verbose = FALSE
  )

  expect_true("perturbed_fixed" %in% names(result@reductions))
  expect_equal(ncol(Seurat::Embeddings(result, "perturbed_fixed")), 3L)
  expect_equal(
    rownames(Seurat::Loadings(result[["perturbed_fixed"]])),
    rownames(Seurat::Loadings(result[["pca"]]))
  )
})


test_that("untouched cells retain their original PCA coordinates", {
  fixture <- .make_projection_fixture()
  result <- ProjectPCAReference(
    fixture$object,
    query_assay = "perturbed",
    dims = 1:3,
    unchanged.cells = fixture$unchanged,
    reduction.name = "perturbed_fixed",
    verbose = FALSE
  )

  observed <- Seurat::Embeddings(result, "perturbed_fixed")[
    fixture$unchanged, , drop = FALSE
  ]
  expected <- Seurat::Embeddings(result, "pca")[
    fixture$unchanged, 1:3, drop = FALSE
  ]
  colnames(expected) <- colnames(observed)
  expect_equal(observed, expected, tolerance = 1e-6)

  changed <- colnames(result)[1:6]
  changed_fixed <- Seurat::Embeddings(result, "perturbed_fixed")[changed, ]
  changed_reference <- Seurat::Embeddings(result, "pca")[changed, 1:3]
  expect_gt(max(abs(changed_fixed - changed_reference)), 0)
})


test_that("projection provenance records the validated reference model", {
  fixture <- .make_projection_fixture()
  result <- ProjectPCAReference(
    fixture$object,
    query_assay = "perturbed",
    dims = 1:3,
    reduction.name = "perturbed_fixed",
    verbose = FALSE
  )

  info <- result[["perturbed_fixed"]]@misc$fixed_reference
  expect_equal(info$reference_assay, "RNA")
  expect_equal(info$query_assay, "perturbed")
  expect_equal(info$dims, 1:3)
  expect_lte(info$validation$scaling_error, 1e-6)
  expect_lte(info$validation$reference_reconstruction_error, 1e-6)
})


test_that("missing scale.data features recommend a clean side-reference PCA", {
  fixture <- .make_projection_fixture()
  object <- fixture$object
  scaled <- tryCatch(
    Seurat::GetAssayData(object, assay = "RNA", slot = "scale.data"),
    error = function(e) {
      Seurat::GetAssayData(object, assay = "RNA", layer = "scale.data")
    }
  )
  scaled <- scaled[-1, , drop = FALSE]
  object <- tryCatch(
    suppressWarnings(Seurat::SetAssayData(
      object, assay = "RNA", slot = "scale.data", new.data = scaled
    )),
    error = function(e) Seurat::SetAssayData(
      object, assay = "RNA", layer = "scale.data", new.data = scaled
    )
  )

  expect_error(
    ProjectPCAReference(
      object,
      query_assay = "perturbed",
      dims = 1:3,
      reduction.name = "perturbed_fixed",
      verbose = FALSE
    ),
    "Recommended action.*pca_fixed_reference"
  )
})


test_that("ProjectPCAReference validates assays, dimensions, and overwrites", {
  fixture <- .make_projection_fixture()

  expect_error(
    ProjectPCAReference(fixture$object, query_assay = "missing"),
    "Query assay"
  )
  expect_error(
    ProjectPCAReference(
      fixture$object, query_assay = "perturbed", dims = 999
    ),
    "dims"
  )
  expect_error(
    ProjectPCAReference(
      fixture$object, query_assay = "perturbed", dims = c(1, 3)
    ),
    "leading consecutive"
  )

  result <- ProjectPCAReference(
    fixture$object,
    query_assay = "perturbed",
    dims = 1:3,
    reduction.name = "perturbed_fixed",
    verbose = FALSE
  )
  expect_error(
    ProjectPCAReference(
      result,
      query_assay = "perturbed",
      dims = 1:3,
      reduction.name = "perturbed_fixed",
      verbose = FALSE
    ),
    "already exists"
  )
})
