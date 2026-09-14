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

test_that("regressor names fall back to the columns of the regressor matrix", {
  object <- fx_var_fitted()
  object[["model"]][["deterministic"]] <- NULL

  # The specification no longer names the deterministic terms, but the data do.
  expect_identical(.get_regressor_names_bvarmodel(object),
                   c("y.l1", "Dp.l1", "r.l1", "const"))
  expect_no_error(summary(object))
})

test_that("regressor names fall back when nothing names the terms", {
  object <- fx_var_fitted()
  object[["model"]][["deterministic"]] <- NULL
  dimnames(object[["data"]][["train"]][["x"]]) <- NULL

  # A short vector of names would be silently misaligned with the coefficients
  # it labels, so the missing ones are generated instead.
  expect_identical(.get_regressor_names_bvarmodel(object),
                   c("y.l1", "Dp.l1", "r.l1", "det.1"))
  expect_no_error(summary(object))
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

test_that("the summary of a period reports the covariance of that period", {
  for (model in list(fx_var_tvp_fitted("sv"), fx_vec_tvp_fitted("sv"))) {
    k <- model[["model"]][["k"]]
    kk <- k^2
    tt <- nrow(model[["data"]][["train"]][["y"]])
    precision <- as.matrix(model[["posterior"]][["u_sigma_inv"]][["coeffs"]])
    covariance_means <- function(period) {
      rowMeans(vapply(seq_len(nrow(precision)),
                      function(i) c(solve(matrix(precision[i, (period - 1) * kk + 1:kk], k))),
                      numeric(kk)))
    }

    expect_equal(as.numeric(summary(model, period = 1)[["sigma"]][["means"]]),
                 covariance_means(1))
    expect_equal(as.numeric(summary(model)[["sigma"]][["means"]]),
                 covariance_means(tt))
  }
})

test_that("a time varying model with one covariance matrix can be summarised and plotted", {
  for (model in list(fx_var_tvp_fitted("gamma"), fx_vec_tvp_fitted("gamma"))) {
    k <- model[["model"]][["k"]]
    precision <- as.matrix(model[["posterior"]][["u_sigma_inv"]][["coeffs"]])
    expect_equal(ncol(precision), k^2)

    covariances <- vapply(seq_len(nrow(precision)),
                          function(i) c(solve(matrix(precision[i, ], k))),
                          numeric(k^2))
    expect_equal(as.numeric(summary(model, period = 1)[["sigma"]][["means"]]),
                 rowMeans(covariances))
    expect_plots(plot(model))
  }
})

test_that("summary of a VEC model reports the cointegration term", {
  summarised <- summary(fx_vec_fitted())

  expect_s3_class(summarised, "summary.bvecmodel")
  expect_output(print(summarised), "Bayesian VEC model")
  expect_output(print(summarised), "l.lr")
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

test_that("covariances that were never estimated carry no significance mark", {
  # A gamma prior on the error variances estimates no covariances, so the
  # off-diagonal entries and both of their credible bounds are exactly zero.
  # Comparing the signs of the two bounds alone marks such a zero as
  # significant, since zero shares its sign with itself.
  output <- utils::capture.output(print(summary(fx_svar_fitted())))
  block <- output[seq(grep("Variance-covariance matrix", output)[1],
                      length(output))]
  marked <- function(row) {
    any(grepl(paste0("^", row, "[[:space:]].*\\*[[:space:]]*$"), block))
  }

  expect_false(marked("y_Dp"))
  expect_false(marked("y_r"))
  expect_false(marked("Dp_r"))
  # A variance is positive by construction, so the ones that are estimated
  # keep their mark.
  expect_true(marked("y_y"))
})

test_that("drawn values of rho are summarised and printed", {
  model <- fx_at_vec_tvp()
  summarised <- summary(model)
  rho <- as.numeric(model[["posterior"]][["beta"]][["rho"]])

  expect_equal(summarised[["rho"]][["means"]], mean(rho))
  expect_equal(summarised[["rho"]][["sd"]], stats::sd(rho))
  expect_equal(summarised[["rho"]][["median"]], stats::median(rho))
  expect_output(print(summarised), "Autocorrelation of the cointegration state equation")

  # The specification says how rho enters before any draw is looked at.
  expect_output(print(model), "drawn from a uniform prior on \\(0.9, 0.999\\), starting at 0.99")
})

test_that("a period is a single period of the sample", {
  model <- fx_at_vec_tvp()
  tt <- nrow(model[["data"]][["train"]][["y"]])

  expect_error(summary(model, period = 1.5), paste0("single integer between 1 and ", tt))
  expect_error(summary(model, period = c(1, 2)), "single integer")
  expect_error(summary(model, period = NA), "single integer")
  expect_error(summary(model, period = tt + 1), "single integer")
  expect_identical(summary(model, period = 2)[["model"]][["period"]], 2L)
})

test_that("the summary of a model with one endogenous variable can be printed", {
  # The variance of a single equation is the only entry of the lower triangle
  # of its covariance table, and selecting it dropped the table to a vector,
  # whose credible bounds print() then looked for as columns that no longer
  # existed. Every AR model and every quantile regression failed this way.
  data("at_macrodata", package = "bvartools", envir = environment())
  y <- at_macrodata[["domestic"]][, "Dp", drop = FALSE] * 100
  set.seed(1234567)
  for (error in c("gamma", "ald")) {
    model <- create_bvarmodel(y, p = 1, deterministic = "const", error = error,
                              iterations = 20, burnin = 5)
    model <- add_priors(model, coef = list(v_i = 0.01),
                        sigma = list(shape = 3, rate = 0.01))
    model <- add_initial_values(model)
    model <- add_posterior_coefficients(model)
    expect_output(print(summary(model)), "Variance:")
  }
})
