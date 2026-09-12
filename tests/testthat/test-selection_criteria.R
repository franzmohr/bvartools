test_that("selection_criteria summarises the in-sample criteria", {
  criteria <- selection_criteria(fx_var_fitted())

  expect_s3_class(criteria, "selcrit")
  expect_named(criteria, c("model", "LL", "AIC", "BIC", "HQ", "WAIC"))
  for (name in c("LL", "AIC", "BIC", "HQ", "WAIC")) {
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
  # The log-likelihood is the one of the whole system, so the penalty counts
  # the free parameters of the whole system: the regressors of one equation
  # times the k equations, and the k(k + 1)/2 of the error covariance.
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] +
    spec[["m"]] * (spec[["s"]] + 1) + spec[["n"]]) +
    spec[["k"]] * (spec[["k"]] + 1) / 2
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

test_that("the penalty of an error correction model grows with the rank", {
  # Pi = alpha beta' is a k x k_ect matrix of rank r and has r(k + k_ect - r)
  # free elements, so the penalty has to grow by k + k_ect - 2r + 1 from one
  # rank to the next. Counting the rank itself, as the method used to, grows it
  # by one, which the gain in fit from an extra cointegration vector exceeds
  # almost always -- the criterion then prefers full rank whatever the data say.
  data("us_macrodata", envir = environment())

  penalties <- vapply(0:2, function(r) {
    object <- create_bvecmodel(data = us_macrodata, p = 2,
                               const = "unrestricted", r = r,
                               iterations = 20, burnin = 10)
    object <- add_priors(object, coef = list(v_i = 1, v_i_det = 1 / 10),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = 3, scale = 1))
    object <- add_initial_values(object)
    object <- add_posterior_loglik(add_posterior_coefficients(object))

    criteria <- selection_criteria(object)
    k <- object[["model"]][["k"]]
    k_ect <- ncol(object[["data"]][["train"]][["w"]])
    n_x <- ncol(object[["data"]][["train"]][["x"]])

    # Recovered from the criterion rather than read off the model, so that this
    # checks what the penalty actually was.
    nobs <- nrow(object[["data"]][["train"]][["y"]])
    nparams <- (criteria[["BIC"]][["mean"]] +
                  2 * criteria[["LL"]][["mean"]]) / log(nobs)

    expect_equal(nparams, r * (k + k_ect - r) + k * n_x + k * (k + 1) / 2)
    nparams
  }, numeric(1))

  expect_true(all(diff(penalties) > 1))
})

test_that("WAIC penalises by the flexibility a fit used", {
  model <- fx_var_fitted()
  criteria <- selection_criteria(model)
  loglik <- model[["posterior"]][["loglik"]]

  # Watanabe's definition, on the deviance scale.
  lppd <- sum(log(colMeans(exp(loglik))))
  p_waic <- sum(apply(loglik, 2, stats::var))

  expect_equal(criteria[["WAIC"]][["mean"]], -2 * (lppd - p_waic))
  expect_equal(criteria[["WAIC"]][["median"]], criteria[["WAIC"]][["mean"]])

  # The quantile columns are a normal interval around the estimate rather than
  # posterior quantiles, so they are symmetric.
  expect_equal(criteria[["WAIC"]][["mean"]] - criteria[["WAIC"]][["qlower"]],
               criteria[["WAIC"]][["qupper"]] - criteria[["WAIC"]][["mean"]])
})

test_that("a single draw leaves WAIC out rather than reporting nonsense", {
  model <- fx_var_fitted()
  model[["posterior"]][["loglik"]] <-
    model[["posterior"]][["loglik"]][1, , drop = FALSE]

  criteria <- selection_criteria(model)

  # The variance across draws is what WAIC penalises with, and one draw has
  # none, so the criterion is absent instead of zero.
  expect_null(criteria[["WAIC"]])
  expect_false("WAIC" %in% names(criteria))
})
