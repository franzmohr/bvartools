test_that("VAR coefficient priors match the number of coefficients", {
  model <- fx_var_priors()
  spec <- model[["model"]]
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])

  expect_named(model[["priors"]], c("a", "u_sigma"))
  expect_identical(model[["priors"]][["a"]][["type"]], "normal")
  expect_identical(dim(model[["priors"]][["a"]][["mu"]]), c(n_coeffs, 1L))
  expect_identical(dim(model[["priors"]][["a"]][["v_inv"]]),
                   c(n_coeffs, n_coeffs))
})

test_that("the uninformative coefficient prior has zero precision", {
  model <- fx_var_priors()

  expect_true(all(model[["priors"]][["a"]][["mu"]] == 0))
  expect_true(all(model[["priors"]][["a"]][["v_inv"]] == 0))
})

test_that("coefficient and deterministic precisions are set separately", {
  model <- add_priors(fx_var_model(),
                      coef = list(v_i = 4, v_i_det = 0.25),
                      sigma = list(df = 1, scale = 0.0001))
  spec <- model[["model"]]
  v_inv <- model[["priors"]][["a"]][["v_inv"]]

  # The coefficient vector holds all lag coefficients first, then all
  # deterministic ones.
  n_lagged <- spec[["k"]] * spec[["k"]] * spec[["p"]]
  n_det <- spec[["k"]] * spec[["n"]]

  expect_length(diag(v_inv), n_lagged + n_det)
  expect_true(all(diag(v_inv)[seq_len(n_lagged)] == 4))
  expect_true(all(diag(v_inv)[n_lagged + seq_len(n_det)] == 0.25))
  # No prior correlation between coefficients.
  expect_true(all(v_inv[upper.tri(v_inv)] == 0))
})

test_that("the error variance prior is stored as given", {
  model <- add_priors(fx_var_model(),
                      coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = 3, scale = 2))

  prior <- model[["priors"]][["u_sigma"]]
  expect_identical(prior[["type"]], "wishart")
  expect_equal(prior[["df"]], 3)
  expect_equal(prior[["scale"]], diag(2, model[["model"]][["k"]]))
})

test_that("a degrees of freedom of 'k' resolves to the number of variables", {
  model <- add_priors(fx_var_model(),
                      coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = "k", scale = 1))

  expect_equal(model[["priors"]][["u_sigma"]][["df"]],
               model[["model"]][["k"]])
})

test_that("VEC priors cover beta, the short-run coefficients and sigma", {
  model <- fx_vec_priors()
  spec <- model[["model"]]

  expect_named(model[["priors"]], c("beta", "a", "u_sigma"))
  expect_identical(model[["priors"]][["beta"]][["type"]], "cointspace")
  expect_identical(dim(model[["priors"]][["beta"]][["p_tau_inv"]]),
                   c(spec[["k_beta"]], spec[["k_beta"]]))

  # The a block holds the loadings alpha followed by the short-run
  # coefficients, which for p = 1 are only the unrestricted deterministic terms.
  n_coeffs <- spec[["k"]] * spec[["rank"]] +
    spec[["k"]] * (spec[["k"]] * (spec[["p"]] - 1) + spec[["n"]])
  expect_equal(dim(model[["priors"]][["a"]][["mu"]]), c(n_coeffs, 1))
})

test_that("priors are added to every model of a modellist", {
  models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                             iterations = 10, burnin = 5)
  models <- add_priors(models, coef = list(v_i = 1, v_i_det = 1),
                       sigma = list(df = 1, scale = 0.0001))

  expect_s3_class(models, "modellist")
  expect_true(all(vapply(models, function(x) !is.null(x[["priors"]]),
                         logical(1))))
  # Longer lag order means more coefficients to shrink.
  expect_lt(nrow(models[[1]][["priors"]][["a"]][["mu"]]),
            nrow(models[[2]][["priors"]][["a"]][["mu"]]))
})

test_that("a time varying cointegration prior needs an autocorrelation", {
  model <- create_bvecmodel(vec_data(), p = 2, r = 1, tvp = TRUE,
                            const = "unrestricted", iterations = 10, burnin = 5)
  args <- list(coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 0.0001),
               sigma = list(df = "k", scale = 1))

  expect_error(
    do.call(add_priors, c(list(model), args, list(coint = list(v_i = 0)))),
    "coint[$]rho")
  expect_error(
    do.call(add_priors, c(list(model), args, list(coint = list(rho = 1)))),
    "smaller than 1")
  expect_no_error(
    do.call(add_priors, c(list(model), args, list(coint = list(rho = 0.999)))))
})
