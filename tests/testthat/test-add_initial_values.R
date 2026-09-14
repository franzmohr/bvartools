test_that("VAR initial values match the prior dimensions", {
  model <- fx_var_initial()
  initial <- model[["initial"]]

  expect_named(initial, c("a", "u_sigma_inv"))
  expect_identical(dim(initial[["a"]]), dim(model[["priors"]][["a"]][["mu"]]))
  expect_identical(dim(initial[["u_sigma_inv"]]),
                   rep(model[["model"]][["k"]], 2L))
  expect_true(all(is.finite(initial[["a"]])))
})

test_that("the initial error precision is symmetric positive definite", {
  sigma_inv <- fx_var_initial()[["initial"]][["u_sigma_inv"]]

  expect_equal(sigma_inv, t(sigma_inv))
  expect_true(all(eigen(sigma_inv, only.values = TRUE)[["values"]] > 0))
})

test_that("the default initial coefficients are the LS estimate", {
  model <- fx_var_initial()
  y <- t(model[["data"]][["train"]][["y"]])
  x <- t(model[["data"]][["train"]][["x"]])

  ls_estimate <- matrix(tcrossprod(y, x) %*% solve(tcrossprod(x)))
  expect_equal(as.numeric(model[["initial"]][["a"]]), as.numeric(ls_estimate))
})

test_that("VEC initial values include the cointegration vectors", {
  model <- fx_vec_initial()
  spec <- model[["model"]]
  initial <- model[["initial"]]

  expect_named(initial, c("beta", "a", "u_sigma_inv"))
  expect_equal(dim(initial[["beta"]]), c(spec[["k_beta"]], spec[["rank"]]))
  expect_identical(dim(initial[["a"]]), dim(model[["priors"]][["a"]][["mu"]]))
  expect_true(all(is.finite(initial[["beta"]])))
})

test_that("initial values can be drawn from the prior instead", {
  # Drawing needs a proper prior, so the uninformative fixture is not usable.
  informative <- add_priors(fx_var_model(),
                            coef = list(v_i = 1, v_i_det = 1),
                            sigma = list(df = 3, scale = 1))

  set.seed(42)
  first <- add_initial_values(informative, method = "prior")
  set.seed(42)
  again <- add_initial_values(informative, method = "prior")
  set.seed(43)
  other <- add_initial_values(informative, method = "prior")

  expect_identical(dim(first[["initial"]][["a"]]),
                   dim(fx_var_initial()[["initial"]][["a"]]))
  expect_true(all(is.finite(first[["initial"]][["a"]])))
  # The draw uses R's RNG, so the seed controls it.
  expect_equal(first[["initial"]][["a"]], again[["initial"]][["a"]])
  expect_false(isTRUE(all.equal(first[["initial"]][["a"]],
                                other[["initial"]][["a"]])))
})

test_that("an unknown method and missing priors are rejected", {
  expect_error(add_initial_values(fx_var_priors(), method = "unknown"),
               "can be 'ols' or 'prior'")
  expect_error(add_initial_values(fx_var_model()), "No information on priors")
  expect_error(add_initial_values(fx_vec_priors(), method = "ols"),
               "can be 'maxlik' or 'prior'")
})

test_that("a time varying model with BVS repeats its draws from the same seed", {
  # The seed is set after add_initial_values(), and not before, on purpose: the
  # state precisions of a TVP model used to be drawn there with no seed, so the
  # chain started somewhere else in every R session.
  for (error in c("ald", "gamma")) {
    model <- create_bvarmodel(var_data(), p = 1, deterministic = "const", tvp = TRUE,
                              error = error, varsel = "bvs",
                              iterations = fx_iterations, burnin = fx_burnin)
    model <- suppressWarnings(
      add_priors(model, coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4),
                 sigma = list(shape = 3, rate = 0.01),
                 varsel = list(inprior = 0.5, covar = FALSE)))

    set.seed(1)
    first <- add_initial_values(model)
    set.seed(2)
    again <- add_initial_values(model)
    expect_identical(first[["initial"]], again[["initial"]])
    expect_equal(diag(first[["initial"]][["a_sigma_inv"]]),
                 as.numeric(model[["priors"]][["a"]][["shape"]] /
                              model[["priors"]][["a"]][["rate"]]))

    set.seed(1000)
    first <- add_posterior_coefficients(first)
    set.seed(1000)
    again <- add_posterior_coefficients(again)
    expect_identical(first[["posterior"]], again[["posterior"]], label = error)
  }
})

test_that("initial values are added to every model of a modellist", {
  models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                             iterations = 10, burnin = 5)
  models <- add_priors(models, coef = list(v_i = 0, v_i_det = 0),
                       sigma = list(df = 1, scale = 0.0001))
  models <- add_initial_values(models)

  expect_s3_class(models, "modellist")
  expect_true(all(vapply(models, function(x) !is.null(x[["initial"]]),
                         logical(1))))
})

test_that("the LS covariance coefficients are stored row by row", {
  # Psi's strict lower triangle is stored row by row, the order the samplers
  # unpack it in. Its LS regressors come out column by column, and the two
  # orders only differ from k = 4 on, so three variables would not catch it.
  levels <- at_domestic()
  data <- stats::ts.intersect(y = diff(levels[, "y"]), Dp = levels[, "Dp"],
                              r = diff(levels[, "r"]), lr = diff(levels[, "lr"])) * 100
  model <- create_bvarmodel(stats::window(data, end = c(1998, 1)), p = 1,
                            deterministic = "const", error = "gamma+covar",
                            iterations = fx_iterations, burnin = fx_burnin)
  model <- add_initial_values(add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                                         sigma = list(shape = 3, rate = 0.0001)))

  y <- t(model[["data"]][["train"]][["y"]])
  x <- t(model[["data"]][["train"]][["x"]])
  u <- y - tcrossprod(y, x) %*% solve(tcrossprod(x)) %*% x
  k <- nrow(u)
  expected <- NULL
  for (i in 2:k) {
    regressors <- -t(u[1:(i - 1), , drop = FALSE])
    expected <- c(expected, solve(crossprod(regressors), crossprod(regressors, u[i, ])))
  }

  expect_equal(k, 4L)
  expect_equal(as.numeric(model[["initial"]][["psi"]]), expected)
})
