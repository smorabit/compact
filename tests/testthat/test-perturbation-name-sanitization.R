test_that("perturbation names are converted to valid Seurat assay names", {
  examples <- data.frame(
    input = c(
      "HuMicA-M1_down",
      "module 1_up",
      "M-1/down",
      "1module_down",
      "red_up"
    ),
    expected = c(
      "HuMicA.M1_down",
      "module.1_up",
      "M.1.down",
      "X1module_down",
      "red_up"
    ),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(examples))){
    input <- examples$input[i]
    expected <- examples$expected[i]

    if(identical(input, expected)){
      expect_silent(
        safe_name <- .sanitize_perturbation_name(input)
      )
    } else {
      expect_warning(
        safe_name <- .sanitize_perturbation_name(input),
        paste0("Using '", expected, "' instead"),
        fixed = TRUE
      )
    }

    expect_identical(safe_name, expected)
  }
})

test_that("perturbation names must be a single non-empty string", {
  expect_error(.sanitize_perturbation_name(character()), "non-empty character string")
  expect_error(.sanitize_perturbation_name(c("red_up", "blue_up")), "non-empty character string")
  expect_error(.sanitize_perturbation_name(""), "non-empty character string")
  expect_error(.sanitize_perturbation_name(NA_character_), "non-empty character string")
})
