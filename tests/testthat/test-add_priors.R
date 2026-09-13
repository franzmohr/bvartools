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

test_that("the prior support of rho is taken in pairs and has to hold rho", {
  model <- create_bvecmodel(vec_data(), p = 2, r = 1, tvp = TRUE,
                            const = "unrestricted", iterations = 10, burnin = 5)
  args <- list(coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 0.0001),
               sigma = list(df = "k", scale = 1))
  priors <- function(coint) do.call(add_priors, c(list(model), args, list(coint = coint)))

  # One end alone would leave the sampler to invent the other.
  expect_error(priors(list(rho = 0.99, rho_min = 0.9)), "specified together")
  expect_error(priors(list(rho = 0.99, rho_max = 0.999)), "specified together")

  expect_error(priors(list(rho = 0.99, rho_min = 0.999, rho_max = 0.9)),
               "rho_min < rho_max")
  expect_error(priors(list(rho = 0.99, rho_min = 0, rho_max = 0.999)),
               "rho_min < rho_max")

  # Drawn, rho is where the chain starts, so a starting value its own prior
  # gives no weight to is a specification that means two things at once.
  expect_error(priors(list(rho = 0.99, rho_min = 0.995, rho_max = 0.999)),
               "must lie between")

  # Without the pair rho stays a hyperparameter, and the prior carries only it.
  fixed <- priors(list(rho = 0.999))
  expect_equal(fixed[["priors"]][["beta"]][["rho"]], 0.999)
  expect_null(fixed[["priors"]][["beta"]][["rho_min"]])

  drawn <- priors(list(rho = 0.99, rho_min = 0.9, rho_max = 0.999))
  expect_equal(drawn[["priors"]][["beta"]][["rho"]], 0.99)
  expect_equal(drawn[["priors"]][["beta"]][["rho_min"]], 0.9)
  expect_equal(drawn[["priors"]][["beta"]][["rho_max"]], 0.999)

  # The prior on the state before the sample is built from the starting value
  # and does not follow the draw; that is the documented departure from Koop et
  # al. (2011), and it is what makes the draw an exact Gibbs block.
  expect_equal(unique(diag(drawn[["priors"]][["beta"]][["v_inv"]])), 1 - 0.99^2)
})

# --- cointegration space prior centred on the ML estimate --------------------

ml_coint_priors <- function(coint, model = fx_vec_model()) {
  add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10), coint = coint,
             sigma = list(df = "k", scale = 1))
}

# The orthonormal basis of the space spanned by Johansen's estimate and of its
# orthogonal complement.
ml_space <- function(model = fx_vec_model()) {
  beta <- bvartools:::.coint_ml(model)[["beta"]]
  basis <- qr.Q(qr(beta), complete = TRUE)
  r <- ncol(beta)
  list(h = basis[, seq_len(r), drop = FALSE],
       h_perp = basis[, -seq_len(r), drop = FALSE])
}

test_that("the ML cointegration estimate is the one used for initial values", {
  expect_equal(matrix(fx_vec_initial()[["initial"]][["beta"]]),
               matrix(bvartools:::.coint_ml(fx_vec_model())[["beta"]]))
})

test_that("a cointegration space prior can be centred on the ML estimate", {
  model <- ml_coint_priors(list(v_i = 0.01, p_tau_i = "ml"))
  p_tau_inv <- model[["priors"]][["beta"]][["p_tau_inv"]]
  space <- ml_space()
  k_beta <- model[["model"]][["k_beta"]]

  expect_equal(model[["priors"]][["beta"]][["v_inv"]], 0.01)
  expect_identical(dim(p_tau_inv), c(k_beta, k_beta))
  expect_equal(p_tau_inv, t(p_tau_inv))
  # Along the estimated space the prior is that of the uniform one ...
  expect_equal(p_tau_inv %*% space[["h"]], space[["h"]])
  # ... and away from it at least as tight.
  expect_true(all(eigen(p_tau_inv, symmetric = TRUE)$values >= 1 - 1e-8))
})

test_that("the weight of the ML prior scales its precision off the space", {
  space <- ml_space()
  off_space <- function(weight) {
    p_tau_inv <- ml_coint_priors(list(v_i = 0.01, p_tau_i = "ml", weight = weight))[["priors"]][["beta"]][["p_tau_inv"]]
    drop(crossprod(space[["h_perp"]], p_tau_inv %*% space[["h_perp"]]))
  }
  expect_equal(off_space(10) / off_space(1), 10)

  # A prior worth next to nothing is capped at the uniform prior rather than
  # turned into one that favours the complement of the estimated space.
  weak <- ml_coint_priors(list(v_i = 0.01, p_tau_i = "ml", weight = 1e-12))
  expect_equal(weak[["priors"]][["beta"]][["p_tau_inv"]], diag(1, 2))
})

test_that("the loading shrinkage can be set from the ML loadings", {
  model <- ml_coint_priors(list(v_i = "ml", p_tau_i = 1))
  v_inv <- model[["priors"]][["beta"]][["v_inv"]]

  expect_length(v_inv, 1)
  expect_true(is.finite(v_inv) && v_inv > 0)
  expect_equal(model[["priors"]][["beta"]][["p_tau_inv"]], diag(1, 2))
})

test_that("a full p_tau_i matrix is taken as it is", {
  p_tau_i <- matrix(c(2, 0.5, 0.5, 1), 2)
  model <- ml_coint_priors(list(v_i = 0.1, p_tau_i = p_tau_i))
  expect_equal(model[["priors"]][["beta"]][["p_tau_inv"]], p_tau_i)
  # Still the old behaviour for a vector of diagonal elements.
  model <- ml_coint_priors(list(v_i = 0.1, p_tau_i = c(1, 2)))
  expect_equal(model[["priors"]][["beta"]][["p_tau_inv"]], diag(c(1, 2)))
})

test_that("an ML cointegration prior is checked", {
  # With zero shrinkage the sampler never sees p_tau_i.
  expect_error(ml_coint_priors(list(v_i = 0, p_tau_i = "ml")), "coint[$]v_i")
  expect_error(ml_coint_priors(list(v_i = 0.1, p_tau_i = "mle")), "coint[$]p_tau_i")
  expect_error(ml_coint_priors(list(v_i = "mle", p_tau_i = 1)), "coint[$]v_i")
  expect_error(ml_coint_priors(list(v_i = 0.1, p_tau_i = "ml", weight = 0)), "coint[$]weight")
  expect_warning(ml_coint_priors(list(v_i = 0.1, p_tau_i = 1, weight = 2)), "coint[$]weight")
  expect_error(ml_coint_priors(list(v_i = 0.1, p_tau_i = diag(3))), "coint[$]p_tau_i")
  expect_error(ml_coint_priors(list(v_i = 0.1, p_tau_i = matrix(c(1, 0, 1, 1), 2))), "symmetric")
})

test_that("scaling is refused once the space prior has a direction", {
  model <- ml_coint_priors(list(v_i = 0.01, p_tau_i = "ml"))
  expect_error(scale_error_correction(model), "before 'add_priors'")
  # A prior without a direction does not care.
  expect_no_error(scale_error_correction(fx_vec_priors()))
})

test_that("a tight ML prior keeps the draws of beta in the estimated space", {
  space <- ml_space()
  model <- ml_coint_priors(list(v_i = "ml", p_tau_i = "ml", weight = 1e6))
  model <- add_initial_values(model)
  set.seed(314159)
  model <- add_posterior_coefficients(model)

  tilt <- function(object) {
    draws <- object[["posterior"]][["beta"]][["coeffs"]]
    abs(draws %*% space[["h_perp"]]) / abs(draws %*% space[["h"]])
  }
  expect_lt(max(tilt(model)), 1e-3)
  expect_lt(stats::median(tilt(model)), stats::median(tilt(fx_vec_fitted())))
})
