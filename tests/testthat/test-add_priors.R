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

test_that("the error prior of a VEC model is stored as given, without the rank", {
  # The constant VECs with a Wishart prior add the rank to the posterior degrees
  # of freedom themselves, for the prior of the loadings given the error
  # covariance (Koop et al., 2010, eq. 8). Adding it to the prior as well
  # counted it twice.
  data <- stats::window(at_data(), end = c(2019, 4)) * 100
  for (r in 1:2) {
    wishart <- create_bvecmodel(data, p = 2, r = r, const = "unrestricted",
                                iterations = 10, burnin = 5)
    wishart <- add_priors(wishart, coef = list(v_i = 0, v_i_det = 0),
                          coint = list(v_i = 0, p_tau_i = 1),
                          sigma = list(df = "k + 1", scale = 1))
    expect_equal(wishart[["priors"]][["u_sigma"]][["df"]], ncol(data) + 1)

    gamma <- create_bvecmodel(data, p = 2, r = r, const = "unrestricted",
                              error = "gamma", iterations = 10, burnin = 5)
    gamma <- add_priors(gamma, coef = list(v_i = 0, v_i_det = 0),
                        coint = list(v_i = 0, p_tau_i = 1),
                        sigma = list(shape = 3, rate = 1))
    expect_equal(as.numeric(gamma[["priors"]][["u_sigma"]][["shape"]]),
                 rep(3, ncol(data)))
  }
})

test_that("add_priors rejects unknown elements of sigma", {
  expect_error(
    add_priors(fx_var_model(), coef = list(v_i = 0, v_i_det = 0),
               sigma = list(df = 3, scale = 2, sclae = 2)),
    "Element 'sclae' in argument 'sigma' is not recognised"
  )
  expect_error(
    add_priors(fx_vec_model(), coef = list(v_i = 1, v_i_det = 1 / 10),
               coint = list(v_i = 0, p_tau_i = 1),
               sigma = list(df = "k", scale = 1, covar = TRUE)),
    "Element 'covar' in argument 'sigma' is not recognised"
  )
})

test_that("add_priors rejects non-positive Wishart degrees of freedom", {
  # The samplers reject df <= 0, so add_priors() must stop before them.
  msg <- "'sigma$df' must be positive"
  expect_error(
    add_priors(fx_var_model(), coef = list(v_i = 0, v_i_det = 0),
               sigma = list(df = 0, scale = 0.0001)),
    msg, fixed = TRUE
  )
  expect_error(
    add_priors(fx_var_model(), coef = list(v_i = 0, v_i_det = 0),
               sigma = list(df = "k - k", scale = 0.0001)),
    msg, fixed = TRUE
  )
  # The samplers take whole degrees of freedom, so a fraction is refused
  # rather than truncated -- to zero, for one below one.
  expect_error(
    add_priors(fx_var_model(), coef = list(v_i = 0, v_i_det = 0),
               sigma = list(df = 0.5, scale = 0.0001)),
    "whole number"
  )
  expect_error(
    add_priors(fx_vec_model(), coef = list(v_i = 1, v_i_det = 1 / 10),
               coint = list(v_i = 0, p_tau_i = 1),
               sigma = list(df = 3.7, scale = 1)),
    "whole number"
  )
  expect_error(
    add_priors(fx_vec_model(), coef = list(v_i = 1, v_i_det = 1 / 10),
               coint = list(v_i = 0, p_tau_i = 1),
               sigma = list(df = 0, scale = 1)),
    msg, fixed = TRUE
  )
  expect_error(
    add_priors(fx_vec_model(), coef = list(v_i = 1, v_i_det = 1 / 10),
               coint = list(v_i = 0, p_tau_i = 1),
               sigma = list(df = -1, scale = 1)),
    msg, fixed = TRUE
  )

  # The smallest accepted value reaches the sampler for both model classes.
  var <- add_priors(fx_var_model(), coef = list(v_i = 0, v_i_det = 0),
                    sigma = list(df = 1, scale = 0.0001))
  expect_identical(var[["priors"]][["u_sigma"]][["df"]], 1L)
  vec <- add_priors(fx_vec_model(), coef = list(v_i = 1, v_i_det = 1 / 10),
                    coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = 1, scale = 1))
  expect_equal(vec[["priors"]][["u_sigma"]][["df"]], 1)
  vec <- add_posterior_coefficients(add_initial_values(vec))
  expect_false(is.null(vec[["posterior"]]))
})

test_that("add_priors rejects unknown elements of coint", {
  expect_error(
    add_priors(fx_vec_model(), coef = list(v_i = 1, v_i_det = 1 / 10),
               coint = list(v_i = 0, p_tau = 1),
               sigma = list(df = "k", scale = 1)),
    "Element 'p_tau' in argument 'coint' is not recognised"
  )
  expect_error(
    add_priors(fx_vec_model(), coef = list(v_i = 1, v_i_det = 1 / 10),
               coint = list(v_i = 0, p_tau_i = 1, rho_mn = 0.9),
               sigma = list(df = "k", scale = 1)),
    "Element 'rho_mn' in argument 'coint' is not recognised"
  )
})

test_that("add_priors rejects elements of coef the model does not use", {
  vec_coint <- list(v_i = 0, p_tau_i = 1)
  vec_sigma <- list(df = "k", scale = 1)

  # coint_var belongs to the Minnesota prior of a VAR model only.
  expect_error(
    add_priors(fx_vec_model(), coef = list(v_i = 1, coint_var = TRUE),
               coint = vec_coint, sigma = vec_sigma),
    "Element 'coint_var' in argument 'coef' is not recognised"
  )
  expect_no_error(
    add_priors(fx_var_model(), coef = list(v_i = 0, coint_var = TRUE),
               sigma = list(df = 1, scale = 0.0001))
  )

  expect_error(
    add_priors(fx_var_model(),
               coef = list(minnesota = list(kappa1 = 0.5, kappa2 = 0.1,
                                            kappa4 = 5, kapa3 = 1)),
               sigma = list(df = 1, scale = 0.0001)),
    "Element 'kapa3' in argument 'coef\\$minnesota' is not recognised"
  )

  expect_error(
    add_priors(fx_var_model(), coef = list(v_i_det = 1),
               sigma = list(df = 1, scale = 0.0001)),
    "If 'coef$v_i' is not specified, 'coef$minnesota' must be specified.",
    fixed = TRUE
  )
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

test_that("the prior of the loadings of a time varying VEC follows coef$v_i", {
  model <- create_bvecmodel(vec_data(), p = 2, r = 1, tvp = TRUE,
                            const = "unrestricted", iterations = 10, burnin = 5)
  model <- add_priors(model,
                      coef = list(v_i = 1 / 2, v_i_det = 0.1, shape = 3, rate = 0.0001),
                      coint = list(rho = 0.999),
                      sigma = list(df = "k", scale = 1))
  n_alpha <- ncol(model[["data"]][["train"]][["y"]]) * model[["model"]][["rank"]]
  precision <- diag(model[["priors"]][["a"]][["v_inv"]])

  # Koop et al. (2011): the loadings' variance is that of the other
  # coefficients times 1 - rho^2, so alpha beta' keeps the scale coef$v_i asks
  # for once beta's stationary variance 1 / (1 - rho^2) is multiplied in.
  expect_equal(precision[1:n_alpha], rep((1 / 2) / (1 - 0.999^2), n_alpha))
  expect_equal(precision[n_alpha + 1], 1 / 2)
})

test_that("the loadings of a time varying VEC can have a state variance rate of their own", {
  model <- create_bvecmodel(vec_data(), p = 2, r = 1, tvp = TRUE,
                            const = "unrestricted", iterations = 10, burnin = 5)
  model <- add_priors(model,
                      coef = list(v_i = 1, v_i_det = 0.1, shape = 3,
                                  rate = 1e-5, rate_alpha = 1e-10, rate_det = 1e-8),
                      coint = list(rho = 0.999),
                      sigma = list(df = "k", scale = 1))
  k <- ncol(model[["data"]][["train"]][["y"]])
  n_alpha <- k * model[["model"]][["rank"]]
  n_det <- k * model[["model"]][["n"]]
  rate <- as.numeric(model[["priors"]][["a"]][["rate"]])
  n_a <- length(rate)

  # The loadings multiply levels and the deterministic terms shift every
  # period, so both need far less drift than the coefficients of the
  # differenced regressors to leave the residuals alone.
  expect_equal(rate[1:n_alpha], rep(1e-10, n_alpha))
  expect_equal(rate[n_a - n_det + 1:n_det], rep(1e-8, n_det))
  expect_equal(rate[(n_alpha + 1):(n_a - n_det)], rep(1e-5, n_a - n_det - n_alpha))

  # Without it the loadings take coef$rate, as before.
  default <- add_priors(create_bvecmodel(vec_data(), p = 2, r = 1, tvp = TRUE,
                                         const = "unrestricted", iterations = 10, burnin = 5),
                        coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-5, rate_det = 1e-8),
                        coint = list(rho = 0.999), sigma = list(df = "k", scale = 1))
  expect_equal(as.numeric(default[["priors"]][["a"]][["rate"]])[1:n_alpha], rep(1e-5, n_alpha))
})

test_that("a VAR model has no loadings to give a rate", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const", tvp = TRUE,
                            iterations = 10, burnin = 5)
  expect_error(add_priors(model,
                          coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-5, rate_alpha = 1e-10),
                          sigma = list(df = "k", scale = 1)),
               "rate_alpha")
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
  set.seed(314159)
  model <- add_initial_values(model)
  model <- add_posterior_coefficients(model)

  tilt <- function(object) {
    draws <- object[["posterior"]][["beta"]][["coeffs"]]
    abs(draws %*% space[["h_perp"]]) / abs(draws %*% space[["h"]])
  }
  expect_lt(max(tilt(model)), 1e-3)
  expect_lt(stats::median(tilt(model)), stats::median(tilt(fx_vec_fitted())))
})

# --- time varying cointegration space centred on the ML estimate -------------

tvp_coint_model <- function() {
  cached_fixture("tvp_coint_model", create_bvecmodel(
    vec_data(), p = 2, r = 1, tvp = TRUE, const = "unrestricted",
    iterations = 10, burnin = 5))
}

tvp_ml_priors <- function(coint) {
  add_priors(tvp_coint_model(),
             coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4),
             coint = coint, sigma = list(df = "k", scale = 1))
}

test_that("a time varying cointegration space can be centred on the ML estimate", {
  rho <- 0.999
  model <- tvp_ml_priors(list(rho = rho, p_tau_i = "ml", weight = 0.1))
  p_tau <- model[["priors"]][["beta"]][["p_tau"]]
  v_inv <- model[["priors"]][["beta"]][["v_inv"]]
  space <- ml_space(tvp_coint_model())
  h <- space[["h"]]
  h_perp <- space[["h_perp"]]

  expect_equal(p_tau, t(p_tau))
  eigenvalues <- eigen(p_tau, symmetric = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10 & eigenvalues <= 1 + 1e-10))
  # Along the estimated space the transition is rho alone ...
  expect_equal(p_tau %*% h, h)
  # ... and off it strictly less.
  expect_lt(max(eigen(crossprod(h_perp, p_tau %*% h_perp), symmetric = TRUE)$values), 1)

  # The state before the sample has the stationary distribution the transition
  # implies: precision (1 - rho^2) H H' + H_perp (I - rho^2 T^2) H_perp'.
  tt <- crossprod(h_perp, p_tau %*% h_perp)
  expect_equal(v_inv, (1 - rho^2) * tcrossprod(h) +
                 h_perp %*% (diag(nrow(tt)) - rho^2 * tt %*% tt) %*% t(h_perp))
  expect_equal(c(model[["priors"]][["beta"]][["mu"]]), rep(0, nrow(v_inv)))
})

test_that("the weight of the time varying ML prior tightens it, down to rho's floor", {
  rho <- 0.999
  space <- ml_space(tvp_coint_model())
  off_space <- function(object) {
    max(eigen(crossprod(space[["h_perp"]],
                        object[["priors"]][["beta"]][["p_tau"]] %*% space[["h_perp"]]),
              symmetric = TRUE)$values)
  }
  weak <- tvp_ml_priors(list(rho = rho, p_tau_i = "ml", weight = 0.01))
  strong <- tvp_ml_priors(list(rho = rho, p_tau_i = "ml", weight = 0.1))
  expect_lt(off_space(strong), off_space(weak))

  # Asking for more than rho allows stops at T = 0 and says so.
  expect_warning(floor <- tvp_ml_priors(list(rho = rho, p_tau_i = "ml", weight = 1e6)),
                 "coint[$]rho")
  expect_equal(floor[["priors"]][["beta"]][["p_tau"]], tcrossprod(space[["h"]]))

  # A weight so small that T is the identity everywhere is no prior on the
  # direction at all, and leaves the one without p_tau_i in place.
  none <- tvp_ml_priors(list(rho = rho, p_tau_i = "ml", weight = 1e-12))
  plain <- tvp_ml_priors(list(rho = rho))
  expect_null(none[["priors"]][["beta"]][["p_tau"]])
  expect_identical(none[["priors"]][["beta"]], plain[["priors"]][["beta"]])
})

test_that("a time varying ML prior is checked", {
  expect_warning(tvp_ml_priors(list(rho = 0.999, p_tau_i = 1)), "coint[$]p_tau_i")
  expect_warning(tvp_ml_priors(list(rho = 0.999, weight = 2)), "coint[$]weight")
  expect_error(tvp_ml_priors(list(rho = 0.999, p_tau_i = "ml", weight = -1)), "coint[$]weight")

  # The transition has a direction, so scaling afterwards is refused too.
  model <- tvp_ml_priors(list(rho = 0.999, p_tau_i = "ml", weight = 0.1))
  expect_error(scale_error_correction(model), "before 'add_priors'")
})

test_that("a gamma error prior can be given one value per equation", {
  model <- create_bvarmodel(diff(at_data()) * 100, p = 1, deterministic = "const",
                            error = "gamma", iterations = 10, burnin = 5)
  priors <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                       sigma = list(shape = c(3, 4, 5), rate = c(1, 2, 3)))[["priors"]][["u_sigma"]]
  expect_equal(as.numeric(priors[["shape"]]), c(3, 4, 5))
  expect_equal(as.numeric(priors[["rate"]]), c(1, 2, 3))
  expect_error(add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                          sigma = list(shape = 3, rate = c(1, 2))),
               "one per endogenous variable")
  expect_error(add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                          sigma = list(shape = 3, rate = c(1, -2, 3))),
               "larger than 0")
})

test_that("a time varying VEC with a Minnesota prior keeps its state equation", {
  # Its shape and rate, or omega_v, were set only in the branch without a
  # Minnesota prior, so add_initial_values() then failed.
  model <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted", tvp = TRUE,
                            error = "gamma", iterations = 10, burnin = 5)
  minnesota <- list(kappa1 = 2, kappa2 = 0.5, kappa4 = 5)
  prior <- function(coef) {
    add_priors(model, coef = coef, coint = list(rho = 0.999),
               sigma = list(shape = 3, rate = 0.01))
  }
  centred <- prior(list(minnesota = minnesota, shape = 3, rate = 1e-4))
  n <- nrow(centred[["priors"]][["a"]][["v_inv"]])
  expect_identical(dim(centred[["priors"]][["a"]][["shape"]]), c(n, 1L))
  expect_identical(dim(centred[["priors"]][["a"]][["rate"]]), c(n, 1L))
  expect_no_error(add_initial_values(centred))

  noncentred <- prior(list(minnesota = minnesota, omega_v = 1e-4))
  expect_identical(dim(noncentred[["priors"]][["a"]][["omega_v"]]), c(n, 1L))

  # The loadings get the same compensating scale as under a plain prior.
  n_alpha <- model[["model"]][["k"]] * model[["model"]][["rank"]]
  constant <- add_priors(create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted",
                                          iterations = 10, burnin = 5),
                         coef = list(minnesota = minnesota),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = "k", scale = 1))
  expect_equal(diag(centred[["priors"]][["a"]][["v_inv"]])[1:n_alpha],
               diag(constant[["priors"]][["a"]][["v_inv"]])[1:n_alpha] / (1 - 0.999^2))
})

test_that("priors that are not priors are refused", {
  # A negative precision on the deterministic terms ran to the end.
  expect_error(add_priors(fx_var_model(), coef = list(v_i = 1, v_i_det = -5),
                          sigma = list(df = "k", scale = 1)),
               "v_i_det")

  # So did a rho below -1, whose stationary variance is negative, until a
  # misleading error about the regressors.
  model <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted", tvp = TRUE,
                            error = "gamma", iterations = 10, burnin = 5)
  expect_error(add_priors(model, coef = list(v_i = 1, shape = 3, rate = 1e-4),
                          coint = list(rho = -1.5), sigma = list(shape = 3, rate = 0.01)),
               "larger than -1")
})

test_that("a Minnesota prior of a model with exogenous variables needs kappa3", {
  endogen <- stats::window(at_macrodata[["domestic"]][, c("y", "Dp")], start = c(1998, 1)) * 100
  exogen <- stats::window(at_macrodata[["foreign"]][, "Dp.s", drop = FALSE], start = c(1998, 1)) * 100
  model <- create_bvarmodel(endogen, p = 1, exogen = exogen, s = 0,
                            deterministic = "const", iterations = 10, burnin = 5)
  minnesota <- list(kappa1 = 0.2, kappa2 = 0.5, kappa4 = 100)
  expect_error(add_priors(model, coef = list(minnesota = minnesota),
                          sigma = list(df = "k", scale = 1)),
               "kappa3")
  expect_no_error(add_priors(model, coef = list(minnesota = c(minnesota, kappa3 = 1)),
                             sigma = list(df = "k", scale = 1)))
})

test_that("scaling after an informative cointegration prior is refused", {
  # v_i sets the prior of alpha given beta, and scaling rescales beta, so the
  # same v_i means a different prior afterwards -- numeric or "ml".
  for (v_i in list(0.5, "ml")) {
    model <- add_priors(fx_vec_model(), coef = list(v_i = 1, v_i_det = 1 / 10),
                        coint = list(v_i = v_i, p_tau_i = 1),
                        sigma = list(df = "k", scale = 1))
    expect_error(scale_error_correction(model), "before 'add_priors'", info = format(v_i))
  }
})

test_that("an error correction term of one series can be scaled", {
  data <- at_macrodata[["domestic"]][, "lr", drop = FALSE] * 100
  model <- create_bvecmodel(data, p = 2, r = 1, const = "restricted",
                            iterations = 10, burnin = 5)
  scaled <- scale_error_correction(model)
  expect_length(attr(scaled[["data"]][["train"]][["w"]], "scale"), 2)
})
