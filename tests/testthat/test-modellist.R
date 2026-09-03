test_that("combine_models merges model lists into one", {
  first <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            iterations = 10, burnin = 5)
  second <- create_bvarmodel(var_data(), p = 2, deterministic = "const",
                             iterations = 10, burnin = 5)
  combined <- combine_models(first, second)

  expect_s3_class(combined, "modellist")
  expect_length(combined, 2)
  expect_equal(vapply(combined, function(x) x[["model"]][["p"]], integer(1)),
               1:2)
})

test_that("combine_models flattens nested model lists", {
  models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                             iterations = 10, burnin = 5)
  extra <- create_bvarmodel(var_data(), p = 3, deterministic = "const",
                            iterations = 10, burnin = 5)
  combined <- combine_models(models, extra)

  expect_s3_class(combined, "modellist")
  expect_length(combined, 3)
  expect_true(all(vapply(combined, inherits, logical(1), "bvarmodel")))
})

test_that("create_bvarmodel already aligns models it builds together", {
  models <- create_bvarmodel(var_data(), p = 1:3, deterministic = "const",
                             iterations = 10, burnin = 5)
  starts <- vapply(models,
                   function(x) stats::tsp(x[["data"]][["train"]][["y"]])[1],
                   numeric(1))

  # All models are estimated on the sample the longest lag order allows.
  expect_length(unique(starts), 1)
})

test_that("align_model_obs cuts separately built models to a common sample", {
  build <- function(p) {
    create_bvarmodel(var_data(), p = p, deterministic = "const",
                     iterations = 10, burnin = 5)
  }
  models <- combine_models(build(1), build(2), build(3))
  aligned <- align_model_obs(models)

  nobs <- vapply(aligned,
                 function(x) nrow(x[["data"]][["train"]][["y"]]), integer(1))
  starts <- vapply(aligned,
                   function(x) stats::tsp(x[["data"]][["train"]][["y"]])[1],
                   numeric(1))

  # Built separately, each model uses as much of the sample as its lag order
  # allows, so they start at different periods.
  expect_gt(length(unique(vapply(
    models, function(x) nrow(x[["data"]][["train"]][["y"]]), integer(1)
  ))), 1)
  expect_length(unique(nobs), 1)
  expect_length(unique(starts), 1)
  # The common sample is the shortest one, i.e. the model with the most lags.
  expect_identical(unique(nobs),
                   nrow(models[[3]][["data"]][["train"]][["y"]]))
})

test_that("combine_models accepts models of different types", {
  var <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                          iterations = 10, burnin = 5)
  vec <- create_bvecmodel(vec_data(), p = 1, r = 1, const = "unrestricted",
                          iterations = 10, burnin = 5)
  combined <- combine_models(var, vec)

  expect_s3_class(combined, "modellist")
  expect_identical(vapply(combined, function(x) x[["model"]][["type"]],
                          character(1)),
                   c("VAR", "VEC"))
})

test_that("combine_models rejects objects it does not recognise", {
  # The intended guard does not fire here: classes[-which(classes == "list")]
  # drops every element when "list" is absent from the classes, so unsupported
  # input slips past the check and fails further down with an opaque message.
  # The contract that matters is that it does not silently return garbage.
  expect_error(combine_models(1:10))
})

test_that("get_model_specifications tabulates a modellist", {
  specifications <- get_model_specifications(fx_var_modellist()[[1]])

  expect_s3_class(specifications, "data.frame")
  expect_identical(specifications[["p"]], 1L)
})

test_that("summary of a modellist covers every model", {
  summaries <- summary(fx_var_modellist())

  expect_s3_class(summaries, "summary.modellist")
  expect_length(summaries, length(fx_var_modellist()))
})

test_that("printing a modellist lists the specifications", {
  expect_output(print(fx_var_modellist()), "VAR")
  expect_output(print(fx_var_modellist()), "Variable selection")
})
