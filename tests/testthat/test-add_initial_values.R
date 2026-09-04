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
