test_that("TargetModulePerturbation is exported with targeting arguments", {
  expect_true("TargetModulePerturbation" %in% getNamespaceExports("compact"))
  expect_true(is.function(TargetModulePerturbation))
  expect_true(all(
    c("target.by", "target_name") %in% names(formals(TargetModulePerturbation))
  ))
})

test_that("the old standalone function name is not exported", {
  expect_false("ModulePerturbationTargeted" %in% getNamespaceExports("compact"))
})

test_that("TargetModulePerturbation uses message-based three-stage progress", {
  function_text <- paste(deparse(body(TargetModulePerturbation)), collapse = "\n")

  expect_match(
    function_text,
    "Applying primary in-silico perturbation to hub genes",
    fixed = TRUE
  )
  expect_match(
    function_text,
    "Applying log-space signal propagation throughout co-expression network",
    fixed = TRUE
  )
  expect_match(
    function_text,
    "Computing cell-cell transition probabilities based on the perturbation",
    fixed = TRUE
  )
  expect_false(grepl("\\bprint\\s*\\(", function_text))
})
