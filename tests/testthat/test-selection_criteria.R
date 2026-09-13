test_that("selection_criteria summarises the in-sample criteria", {
  criteria <- selection_criteria(fx_var_fitted())

  expect_s3_class(criteria, "selcrit")
  expect_named(criteria, c("model", "LL", "AIC", "BIC", "HQ", "WAIC", "LOOIC"))
  for (name in c("LL", "AIC", "BIC", "HQ", "WAIC", "LOOIC")) {
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

  # The penalties are added to the deviance at the point estimate. The mean of
  # the deviance over the posterior is the larger quantity, by about the
  # effective number of parameters, so using it would charge the complexity of
  # the model once through the penalty and again through the averaging.
  loglik <- model[["posterior"]][["loglik"]]
  deviance <- -2 * sum(colMeans(loglik)) - sum(apply(loglik, 2, stats::var))

  expect_equal(criteria[["AIC"]][["mean"]], deviance + 2 * n_coeffs)
  expect_equal(criteria[["BIC"]][["mean"]], deviance + log(nobs) * n_coeffs)
  expect_equal(criteria[["HQ"]][["mean"]],
               deviance + 2 * log(log(nobs)) * n_coeffs)

  expect_lt(deviance, -2 * mean(rowSums(loglik)))
})

test_that("the parameter counting criteria are reported without a band", {
  criteria <- selection_criteria(fx_var_fitted())

  # A criterion is a function of the data and the estimator, not a parameter
  # with a posterior, so it has a value and no credible interval.
  for (name in c("AIC", "BIC", "HQ")) {
    expect_equal(criteria[[name]][["median"]], criteria[[name]][["mean"]])
    expect_true(is.na(criteria[[name]][["qlower"]]))
    expect_true(is.na(criteria[[name]][["qupper"]]))
  }
})

test_that("the criteria are ordered by how hard they penalise size", {
  criteria <- selection_criteria(fx_var_fitted())

  # With 74 observations log(T) > 2 log(log(T)) > 2, so BIC penalises hardest.
  expect_lt(criteria[["AIC"]][["mean"]], criteria[["HQ"]][["mean"]])
  expect_lt(criteria[["HQ"]][["mean"]], criteria[["BIC"]][["mean"]])
})

test_that("a longer lag order is penalised more for the same fit", {
  criteria <- selection_criteria(fx_var_modellist())
  # Both criteria are the same deviance plus a multiple of the number of
  # parameters, so their difference isolates the penalty without the fit.
  penalty <- function(x) x[["BIC"]][["mean"]] - x[["AIC"]][["mean"]]

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

test_that("the default criterion is the one that stays defined for every model", {
  criteria <- selection_criteria(fx_var_modellist())

  # WAIC penalises by the flexibility the fit used, which is the only penalty
  # that means the same thing for a constant, a shrunk and a time varying
  # specification, so it is what a comparison defaults to.
  expect_equal(choose_best_model(criteria),
               choose_best_model(criteria, criterion = "WAIC"))
  expect_identical(formals(choose_best_model.selcritlist)[["criterion"]], "WAIC")

  # Choosing and plotting should not disagree about what they show by default
  expect_identical(formals(plot.selcritlist)[["criterion"]],
                   formals(choose_best_model.selcritlist)[["criterion"]])
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

    # Recovered from the criteria rather than read off the model, so that this
    # checks what the penalty actually was. BIC and AIC differ by the penalty
    # alone, which leaves the number of parameters behind.
    nobs <- nrow(object[["data"]][["train"]][["y"]])
    nparams <- (criteria[["BIC"]][["mean"]] -
                  criteria[["AIC"]][["mean"]]) / (log(nobs) - 2)

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

test_that("LOOIC reweights the posterior towards the data it has not seen", {
  model <- fx_var_fitted()
  criteria <- selection_criteria(model)
  loglik <- as.matrix(model[["posterior"]][["loglik"]])

  loo <- .psis_loo(loglik)
  expect_equal(criteria[["LOOIC"]][["mean"]], loo[["looic"]])
  expect_equal(criteria[["LOOIC"]][["median"]], criteria[["LOOIC"]][["mean"]])

  # As for WAIC the band is a normal interval around the estimate rather than
  # posterior quantiles, so it is symmetric.
  expect_equal(criteria[["LOOIC"]][["mean"]] - criteria[["LOOIC"]][["qlower"]],
               criteria[["LOOIC"]][["qupper"]] - criteria[["LOOIC"]][["mean"]])

  # Cross validation cannot flatter the model more than the fit it is derived
  # from, and it estimates the same effective number of parameters that WAIC
  # penalises with.
  expect_gt(loo[["looic"]], -2 * sum(apply(loglik, 2, function(z) {
    .log_sum_exp(z) - log(nrow(loglik))
  })))
  expect_equal(loo[["p_loo"]], sum(apply(loglik, 2, stats::var)), tolerance = 0.3)
})

test_that("the importance weights are smoothed and normalised", {
  set.seed(20240101)

  # A heavy tailed set of ratios is what the smoothing exists for: the largest
  # weights are pulled in, so the sum stays one and no single draw dominates.
  smoothed <- .psis_weights(rt(2000, df = 2), tail_len = 134)
  weights <- exp(smoothed[["log_weights"]])

  expect_equal(sum(weights), 1)
  expect_true(is.finite(smoothed[["pareto_k"]]))
  expect_lt(max(weights), 1)

  # A tail with nothing to fit is left alone rather than smoothed into noise,
  # and is flagged by an infinite shape parameter.
  expect_identical(.psis_weights(rep(0, 500), tail_len = 100)[["pareto_k"]], Inf)
  expect_identical(.psis_weights(stats::rnorm(20), tail_len = 4)[["pareto_k"]], Inf)
})

test_that("the Pareto fit recovers the shape of a known tail", {
  set.seed(20240102)
  rgpd <- function(n, k, sigma) sigma * expm1(-k * log1p(-stats::runif(n))) / k

  for (shape in c(0.2, 0.7)) {
    estimates <- replicate(20, .gpd_fit(rgpd(500, shape, 2))[["k"]])
    expect_lt(abs(mean(estimates) - shape), 0.05)
  }
})

test_that("a single draw leaves the cross validated criteria out", {
  model <- fx_var_fitted()
  model[["posterior"]][["loglik"]] <-
    model[["posterior"]][["loglik"]][1, , drop = FALSE]

  criteria <- selection_criteria(model)

  # There is nothing to reweight with a single draw, and that draw is its own
  # point estimate, so the criteria that count parameters remain.
  expect_null(criteria[["LOOIC"]])
  expect_false(is.na(criteria[["AIC"]][["mean"]]))
})

test_that("a pointwise log-likelihood too variable for its correction is reported", {
  model <- fx_var_fitted()
  shape <- dim(as.matrix(model[["posterior"]][["loglik"]]))
  set.seed(20260913)

  # A variance far below the threshold in every period ...
  calm <- matrix(stats::rnorm(prod(shape), mean = -2, sd = 0.1), shape[1])
  model[["posterior"]][["loglik"]] <- calm
  criteria <- selection_criteria(model)
  expect_identical(attr(criteria[["WAIC"]], "n_high_var"), 0L)
  expect_false(any(grepl("log-likelihood variance", utils::capture.output(print(criteria)))))

  # ... and far above it in two periods.
  wild <- calm
  wild[, 1:2] <- stats::rnorm(2 * shape[1], mean = -2, sd = 3)
  model[["posterior"]][["loglik"]] <- wild
  criteria <- selection_criteria(model)
  expect_identical(attr(criteria[["WAIC"]], "n_high_var"), 2L)
  expect_output(print(criteria),
                "2 periods have a pointwise log-likelihood variance above 0.4")

  # A list of criteria names the model concerned.
  models <- fx_var_modellist()
  for (i in seq_along(models)) {
    shape_i <- dim(as.matrix(models[[i]][["posterior"]][["loglik"]]))
    models[[i]][["posterior"]][["loglik"]] <-
      matrix(stats::rnorm(prod(shape_i), mean = -2, sd = if (i == 2) 3 else 0.1), shape_i[1])
  }
  expect_output(print(selection_criteria(models)), "in model 2 \\(")
})
