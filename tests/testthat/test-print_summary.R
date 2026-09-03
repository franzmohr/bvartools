test_that("print tabulates the specification of a VAR model", {
  output <- utils::capture.output(print(fx_var_model()))

  expect_match(output[1], "Type")
  expect_match(output[1], "Variable selection")
  expect_match(output[2], "VAR")
  # The number of observations of the estimation sample is reported.
  expect_match(output[2], as.character(nrow(
    fx_var_model()[["data"]][["train"]][["y"]])))
})

test_that("print tabulates the specification of a VEC model", {
  output <- utils::capture.output(print(fx_vec_model()))

  expect_match(output[1], "rank")
  expect_match(output[1], "n_restricted")
  expect_match(output[2], "VEC")
})

test_that("summary of a VAR model reports one block per equation", {
  summarised <- summary(fx_var_fitted())

  expect_s3_class(summarised, "summary.bvarmodel")
  expect_output(print(summarised), "Bayesian VAR model")
  for (variable in fx_var_fitted()[["model"]][["endogen"]]) {
    expect_output(print(summarised), variable)
  }
  expect_named(summarised, c("model", "a", "sigma"))
  expect_output(print(summarised), "Variance-covariance matrix")
})

test_that("the coefficient summary matches the posterior draws", {
  model <- fx_var_fitted()
  summarised <- summary(model)
  draws <- model[["posterior"]][["a"]][["coeffs"]]
  k <- model[["model"]][["k"]]
  n_regressors <- k * model[["model"]][["p"]] + model[["model"]][["n"]]

  # Equations are rows and regressors columns, while the draws hold the
  # coefficient matrix stacked column by column.
  expect_identical(dim(summarised[["a"]][["means"]]), c(k, n_regressors))
  expect_equal(as.numeric(summarised[["a"]][["means"]]),
               unname(colMeans(draws)))
  expect_equal(as.numeric(summarised[["a"]][["sd"]]),
               unname(apply(draws, 2, stats::sd)))
  expect_equal(as.numeric(summarised[["a"]][["median"]]),
               unname(apply(draws, 2, stats::median)))
  expect_equal(as.numeric(summarised[["a"]][["q_lower"]]),
               unname(apply(draws, 2, stats::quantile, probs = 0.025)))
  expect_equal(as.numeric(summarised[["a"]][["q_upper"]]),
               unname(apply(draws, 2, stats::quantile, probs = 0.975)))
})

test_that("the summary reports the covariance, not the precision", {
  model <- fx_var_fitted()
  summarised <- summary(model)
  k <- model[["model"]][["k"]]
  precision <- model[["posterior"]][["u_sigma_inv"]][["coeffs"]]

  covariances <- vapply(seq_len(nrow(precision)),
                        function(i) c(solve(matrix(precision[i, ], k))),
                        numeric(k^2))

  expect_equal(as.numeric(summarised[["sigma"]][["means"]]),
               rowMeans(covariances))
  # A covariance matrix is symmetric.
  expect_equal(summarised[["sigma"]][["means"]],
               t(summarised[["sigma"]][["means"]]))
})

test_that("summary of a VEC model reports the cointegration term", {
  summarised <- summary(fx_vec_fitted())

  expect_s3_class(summarised, "summary.bvecmodel")
  expect_output(print(summarised), "Bayesian VEC model")
  expect_output(print(summarised), "l.R")
})

test_that("the credible interval of the summary can be set", {
  narrow <- summary(fx_var_fitted(), ci = 0.5)
  wide <- summary(fx_var_fitted(), ci = 0.95)

  expect_output(print(narrow), "25%")
  expect_output(print(wide), "2.5%")
})

test_that("printing selection criteria shows all four", {
  output <- utils::capture.output(print(selection_criteria(fx_var_fitted())))

  expect_true(any(grepl("LL", output)))
  expect_true(any(grepl("AIC", output)))
  expect_true(any(grepl("BIC", output)))
  expect_true(any(grepl("HQ", output)))
})

test_that("printing a list of selection criteria shows every model", {
  output <- utils::capture.output(print(selection_criteria(fx_var_modellist())))

  expect_true(any(grepl("Model 1", output)))
  expect_true(any(grepl("Model 2", output)))
})
