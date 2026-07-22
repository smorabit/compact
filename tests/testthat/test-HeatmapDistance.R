distance_matrix_original <- data.frame(
  DMSO = c(0, 1),
  Bortezomib = c(1, 0),
  row.names = c("DMSO", "Bortezomib")
)

distance_matrix_perturbed <- data.frame(
  DMSO = c(0, 3),
  Bortezomib = c(3, 0),
  row.names = c("DMSO", "Bortezomib")
)

test_that("HeatmapDistance applies explicit shared color limits", {
  p <- HeatmapDistance(
    distance_matrix_original,
    distance_matrix_perturbed,
    min_val = 0,
    max_val = 2
  )

  expect_s3_class(p, "patchwork")
  expect_equal(p[[1]]$scales$get_scales("fill")$limits, c(0, 2))
  expect_equal(p[[2]]$scales$get_scales("fill")$limits, c(0, 2))
})

test_that("HeatmapDistance accepts one explicit limit and derives the other", {
  p_min <- HeatmapDistance(
    distance_matrix_original,
    distance_matrix_perturbed,
    min_val = -1
  )
  p_max <- HeatmapDistance(
    distance_matrix_original,
    distance_matrix_perturbed,
    max_val = 4
  )

  expect_equal(p_min[[1]]$scales$get_scales("fill")$limits, c(-1, 3))
  expect_equal(p_max[[1]]$scales$get_scales("fill")$limits, c(0, 4))
})

test_that("HeatmapDistance validates explicit color limits", {
  expect_error(
    HeatmapDistance(
      distance_matrix_original,
      distance_matrix_perturbed,
      min_val = c(0, 1)
    ),
    "min_val must be NULL or a single finite numeric value",
    fixed = TRUE
  )
  expect_error(
    HeatmapDistance(
      distance_matrix_original,
      distance_matrix_perturbed,
      max_val = Inf
    ),
    "max_val must be NULL or a single finite numeric value",
    fixed = TRUE
  )
  expect_error(
    HeatmapDistance(
      distance_matrix_original,
      distance_matrix_perturbed,
      min_val = 2,
      max_val = 2
    ),
    "min_val must be less than max_val",
    fixed = TRUE
  )
})

test_that("HeatmapDistance squishes values outside explicit limits", {
  p <- HeatmapDistance(
    distance_matrix_original,
    distance_matrix_perturbed,
    min_val = 0,
    max_val = 2
  )
  fill_scale <- p[[2]]$scales$get_scales("fill")

  expect_equal(fill_scale$oob(c(-1, 1, 3), c(0, 2)), c(0, 1, 2))
})
