# The prior helpers and the analysis functions over collections of models: a
# 'modellist' and the windows of an 'expandingwindow'. What is checked is that
# each model gets its own result, the same as a call on it alone; the methods
# for single models are tested in their own files.

test_that("the prior helpers give every model of a collection its own prior", {
  models <- fx_var_modellist()
  windows <- fx_expanding_window()
  for (helper in list(minnesota_prior, ssvs_prior, inclusion_prior)) {
    for (collection in list(models, windows)) {
      priors <- helper(collection)
      expect_length(priors, length(collection))
      for (i in seq_along(collection)) {
        expect_equal(priors[[i]], helper(collection[[i]]))
      }
    }
  }
})

test_that("prior options reach the windows of an expanding window exercise", {
  windows <- fx_expanding_window()
  shrunk <- add_prior_options(windows, shrinkage = "horseshoe")
  expect_s3_class(shrunk, "expandingwindow")
  for (i in seq_along(windows)) {
    expect_identical(shrunk[[i]][["model"]][["shrinkage"]], "horseshoe")
  }
})

test_that("historical decompositions and multipliers are taken model by model", {
  models <- fx_var_modellist()
  hd <- historical_decomposition(models, response = "y")
  expect_length(hd, length(models))
  expect_equal(hd[[2]], historical_decomposition(models[[2]], response = "y"))

  windows <- fx_expanding_window()
  hd <- historical_decomposition(windows, response = "y")
  expect_length(hd, length(windows))

  # A failure names the model it happened in.
  expect_error(historical_decomposition(models, response = "nope"), "of model '")
  expect_error(multipliers(models, impulse = "x", response = "y"), "of model '")
})
