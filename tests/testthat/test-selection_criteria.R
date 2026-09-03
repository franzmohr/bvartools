test_that("selection_criteria summarises the four in-sample criteria", {
  criteria <- selection_criteria(fx_var_fitted())

  expect_s3_class(criteria, "selcrit")
  expect_named(criteria, c("model", "LL", "AIC", "BIC", "HQ"))
  for (name in c("LL", "AIC", "BIC", "HQ")) {
    expect_named(criteria[[name]], c("mean", "median", "qlower", "qupper"))
    expect_identical(nrow(criteria[[name]]), 1L)
  }
  expect_identical(attr(criteria, "ci"), c("2.5%", "97.5%"))
})

test_that("the log-likelihood summary matches the stored draws", {
  model <- fx_var_fitted()
  criteria <- selection_criteria(model)
  loglik <- rowSums(model[["posterior"]][["loglik"]])

  expect_equal(criteria[["LL"]][["mean"]], mean(loglik))
  expect_equal(criteria[["LL"]][["median"]], stats::median(loglik))
  expect_equal(criteria[["LL"]][["qlower"]],
               unname(stats::quantile(loglik, 0.025)))
  expect_equal(criteria[["LL"]][["qupper"]],
               unname(stats::quantile(loglik, 0.975)))
})

test_that("the information criteria follow their definitions", {
  model <- fx_var_fitted()
  spec <- model[["model"]]
  criteria <- selection_criteria(model)

  nobs <- nrow(model[["data"]][["train"]][["y"]])
  # The penalty counts the regressors of one equation, not the k * (kp + n)
  # coefficients of the whole system.
  n_coeffs <- spec[["k"]] * spec[["p"]] +
    spec[["m"]] * (spec[["s"]] + 1) + spec[["n"]]
  loglik <- criteria[["LL"]][["mean"]]

  expect_equal(criteria[["AIC"]][["mean"]], -2 * loglik + 2 * n_coeffs)
  expect_equal(criteria[["BIC"]][["mean"]], -2 * loglik + log(nobs) * n_coeffs)
  expect_equal(criteria[["HQ"]][["mean"]],
               -2 * loglik + 2 * log(log(nobs)) * n_coeffs)
})

test_that("the criteria are ordered by how hard they penalise size", {
  criteria <- selection_criteria(fx_var_fitted())

  # With 74 observations log(T) > 2 log(log(T)) > 2, so BIC penalises hardest.
  expect_lt(criteria[["AIC"]][["mean"]], criteria[["HQ"]][["mean"]])
  expect_lt(criteria[["HQ"]][["mean"]], criteria[["BIC"]][["mean"]])
})

test_that("a longer lag order is penalised more for the same fit", {
  criteria <- selection_criteria(fx_var_modellist())
  penalty <- function(x) x[["BIC"]][["mean"]] + 2 * x[["LL"]][["mean"]]

  # The penalty term alone has to grow with the number of coefficients.
  expect_lt(penalty(criteria[[1]]), penalty(criteria[[2]]))
})

test_that("selection_criteria works for a modellist", {
  criteria <- selection_criteria(fx_var_modellist())

  expect_s3_class(criteria, "selcritlist")
  expect_length(criteria, length(fx_var_modellist()))
  expect_true(all(vapply(criteria, inherits, logical(1), "selcrit")))
})

test_that("choose_best_model picks the model with the lowest criterion", {
  criteria <- selection_criteria(fx_var_modellist())

  for (criterion in c("AIC", "BIC", "HQ")) {
    means <- vapply(criteria, function(x) x[[criterion]][["mean"]], numeric(1))
    expect_equal(as.integer(choose_best_model(criteria, criterion = criterion)),
                 which.min(means))
  }
  # The log-likelihood is maximised rather than minimised.
  loglik <- vapply(criteria, function(x) x[["LL"]][["mean"]], numeric(1))
  expect_equal(as.integer(choose_best_model(criteria, criterion = "LL")),
               which.max(loglik))
})

test_that("get_model_specifications describes the underlying model", {
  spec <- get_model_specifications(selection_criteria(fx_var_fitted()))

  expect_s3_class(spec, "data.frame")
  expect_identical(spec[["type"]], "VAR")
  expect_identical(spec[["k"]], 3L)
  expect_identical(spec[["p"]], 1L)
})
