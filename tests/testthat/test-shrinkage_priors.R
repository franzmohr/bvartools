test_that("minnesota_prior returns prior moments for all coefficients", {
  model <- fx_var_model()
  spec <- model[["model"]]
  prior <- minnesota_prior(model, kappa1 = 0.5, kappa2 = 0.1)
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])

  expect_type(prior, "list")
  expect_named(prior, c("mu", "v_inv", "sigma_inv"))
  expect_identical(dim(prior[["mu"]]), c(n_coeffs, 1L))
  expect_identical(dim(prior[["v_inv"]]), c(n_coeffs, n_coeffs))
  expect_identical(dim(prior[["sigma_inv"]]), c(spec[["k"]], spec[["k"]]))
  # Random walk prior: all prior means are zero for differenced data.
  expect_true(all(prior[["mu"]] == 0))
  expect_true(all(diag(prior[["v_inv"]]) > 0))
  expect_true(all(prior[["v_inv"]][upper.tri(prior[["v_inv"]])] == 0))
})

test_that("minnesota_prior follows the Minnesota variance formula", {
  kappa1 <- 0.5
  kappa2 <- 0.1
  kappa4 <- 5
  spec <- fx_var_model()[["model"]]
  k <- spec[["k"]]
  prior <- minnesota_prior(fx_var_model(), kappa1 = kappa1, kappa2 = kappa2,
                           kappa4 = kappa4)
  precision <- diag(prior[["v_inv"]])

  # sigma_inv holds the inverse residual variances of the univariate AR models
  # the prior variances are scaled with.
  sigma_sq <- 1 / diag(prior[["sigma_inv"]])

  # The k x kp coefficient matrix is stored column by column, so the
  # coefficient of variable j in equation l enters at ((j - 1) * k + l).
  expected <- matrix(NA_real_, k, k)
  for (l in seq_len(k)) {
    for (j in seq_len(k)) {
      expected[l, j] <- if (l == j) {
        1 / kappa1
      } else {
        1 / (kappa1 * kappa2) * sigma_sq[j] / sigma_sq[l]
      }
    }
  }
  expect_equal(precision[seq_len(k^2)], as.numeric(expected))

  # Deterministic terms are scaled with kappa4 and the equation's own variance.
  expect_equal(precision[k^2 + seq_len(k)],
               1 / (kappa1 * kappa4 * sigma_sq))
})

test_that("a smaller kappa2 shrinks cross-variable lags harder", {
  k <- fx_var_model()[["model"]][["k"]]
  loose <- diag(minnesota_prior(fx_var_model(), kappa1 = 0.5,
                                kappa2 = 0.5)[["v_inv"]])
  tight <- diag(minnesota_prior(fx_var_model(), kappa1 = 0.5,
                                kappa2 = 0.05)[["v_inv"]])

  own <- (seq_len(k) - 1) * k + seq_len(k)
  cross <- setdiff(seq_len(k^2), own)

  expect_true(all(tight[cross] > loose[cross]))
  # kappa2 leaves own lags untouched.
  expect_equal(tight[own], loose[own])
})

test_that("minnesota_prior shrinks more distant lags harder", {
  model <- create_bvarmodel(var_data(), p = 2, deterministic = "const",
                            iterations = 10, burnin = 5)
  k <- model[["model"]][["k"]]
  precision <- diag(minnesota_prior(model, kappa1 = 0.5,
                                    kappa2 = 0.5)[["v_inv"]])

  first_lag <- precision[seq_len(k^2)]
  second_lag <- precision[k^2 + seq_len(k^2)]
  expect_true(all(second_lag > first_lag))
})

test_that("minnesota_prior rejects non-positive shrinkage parameters", {
  expect_error(minnesota_prior(fx_var_model(), kappa1 = 0), "must be positive")
  expect_error(minnesota_prior(fx_var_model(), kappa2 = -1), "must be positive")
  expect_error(minnesota_prior(fx_var_model(), kappa3 = 0), "must be positive")
})

# A short training sample can have fewer observations than the VAR form has
# regressors per equation. The default sigma = "AR" never needs that estimate --
# it regresses each variable on its own lags alone -- so it has to keep working
# there. Only sigma = "VAR" is entitled to complain, and it should say why.
short_var_model <- function(end = c(1982, 1), p = 4) {
  create_bvarmodel(stats::window(var_data(), end = end), p = p,
                   deterministic = "const", iterations = 10, burnin = 5)
}

test_that("minnesota_prior works when the VAR form has too few observations", {
  model <- short_var_model()
  spec <- model[["model"]]
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])
  prior <- minnesota_prior(model, kappa1 = 0.5, kappa2 = 0.1)

  expect_identical(dim(prior[["mu"]]), c(n_coeffs, 1L))
  expect_length(diag(prior[["v_inv"]]), n_coeffs)
  expect_true(all(is.finite(diag(prior[["v_inv"]]))))
  expect_true(all(diag(prior[["v_inv"]]) > 0))
})

test_that("minnesota_prior says why sigma = 'VAR' needs a longer sample", {
  expect_error(minnesota_prior(short_var_model(), sigma = "VAR"),
               "regressors per equation")
})

test_that("minnesota_prior says why its own AR regressions need a longer sample", {
  expect_error(minnesota_prior(short_var_model(end = c(1981, 2))),
               "such regressors")
})

test_that("minnesota_prior still returns a VAR error covariance when it can", {
  k <- fx_var_model()[["model"]][["k"]]
  prior <- minnesota_prior(fx_var_model(), kappa1 = 0.5, kappa2 = 0.1,
                           sigma = "VAR")

  expect_identical(dim(prior[["sigma_inv"]]), c(k, k))
  expect_true(all(is.finite(prior[["sigma_inv"]])))
})

test_that("the same sample rules hold for VEC models", {
  short_vec <- create_bvecmodel(stats::window(vec_data(), end = c(1983, 2)),
                                p = 6, r = 1, const = "unrestricted",
                                iterations = 10, burnin = 5)

  expect_type(minnesota_prior(short_vec, kappa1 = 0.5, kappa2 = 0.1), "list")
  expect_error(minnesota_prior(short_vec, kappa1 = 0.5, kappa2 = 0.1,
                               sigma = "VAR"),
               "regressors per equation")
})

test_that("minnesota_prior works for VEC models", {
  prior <- minnesota_prior(fx_vec_model(), kappa1 = 0.5, kappa2 = 0.1)

  expect_type(prior, "list")
  expect_true(all(prior[["mu"]] == 0))
  expect_true(all(diag(prior[["v_i"]]) > 0))
})

test_that("add_priors can build the coefficient prior from Minnesota", {
  model <- add_priors(fx_var_model(),
                      coef = list(minnesota = list(kappa1 = 0.5, kappa2 = 0.1,
                                                   kappa4 = 5)),
                      sigma = list(df = 1, scale = 0.0001))
  direct <- minnesota_prior(fx_var_model(), kappa1 = 0.5, kappa2 = 0.1,
                            kappa4 = 5)

  expect_equal(model[["priors"]][["a"]][["v_inv"]], direct[["v_inv"]])
})

test_that("add_priors caps the Minnesota prior variances of a VAR at max_var", {
  minnesota <- list(kappa1 = 0.5, kappa2 = 0.1, kappa4 = 5)
  sigma <- list(df = 1, scale = 0.0001)
  uncapped <- add_priors(fx_var_model(), coef = list(minnesota = minnesota),
                         sigma = sigma)
  # The smallest prior variance, so that every larger one is cut back.
  max_var <- min(1 / diag(uncapped[["priors"]][["a"]][["v_inv"]]))

  capped <- add_priors(fx_var_model(),
                       coef = list(minnesota = minnesota, max_var = max_var),
                       sigma = sigma)
  direct <- minnesota_prior(fx_var_model(), kappa1 = 0.5, kappa2 = 0.1,
                            kappa4 = 5, max_var = max_var)

  expect_equal(capped[["priors"]][["a"]][["v_inv"]], direct[["v_inv"]])
  expect_false(isTRUE(all.equal(capped[["priors"]][["a"]][["v_inv"]],
                                uncapped[["priors"]][["a"]][["v_inv"]])))
})

test_that("add_priors caps the Minnesota prior variances of a VEC model at max_var", {
  # With p = 1 a VEC model has no short-run lags for the cap to apply to.
  model <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted",
                            iterations = fx_iterations, burnin = fx_burnin)
  minnesota <- list(kappa1 = 0.5, kappa2 = 0.1, kappa4 = 5)
  coint <- list(v_i = 0, p_tau_i = 1)
  sigma <- list(df = "k", scale = 1)
  uncapped <- add_priors(model, coef = list(minnesota = minnesota),
                         coint = coint, sigma = sigma)
  max_var <- min(1 / diag(uncapped[["priors"]][["a"]][["v_inv"]]))

  capped <- add_priors(model,
                       coef = list(minnesota = minnesota, max_var = max_var),
                       coint = coint, sigma = sigma)
  direct <- minnesota_prior(model, kappa1 = 0.5, kappa2 = 0.1, kappa4 = 5,
                            max_var = max_var)

  expect_equal(capped[["priors"]][["a"]][["v_inv"]], direct[["v_i"]])
  expect_false(isTRUE(all.equal(capped[["priors"]][["a"]][["v_inv"]],
                                uncapped[["priors"]][["a"]][["v_inv"]])))
})

test_that("add_priors requires the mandatory Minnesota parameters", {
  expect_error(
    add_priors(fx_var_model(),
               coef = list(minnesota = list(kappa1 = 0.5)),
               sigma = list(df = 1, scale = 0.0001)),
    "kappa1"
  )
})

test_that("ssvs_prior returns spike and slab standard deviations", {
  spec <- fx_var_model()[["model"]]
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])
  prior <- ssvs_prior(fx_var_model(), tau = c(0.05, 10))

  expect_named(prior, c("tau0", "tau1"))
  expect_identical(dim(prior[["tau0"]]), c(n_coeffs, 1L))
  expect_true(all(prior[["tau0"]] == 0.05))
  expect_true(all(prior[["tau1"]] == 10))
  # The spike has to be tighter than the slab for the mixture to select.
  expect_true(all(prior[["tau0"]] < prior[["tau1"]]))
})

test_that("the semiautomatic SSVS prior scales with coefficient uncertainty", {
  prior <- ssvs_prior(fx_var_model(), semiautomatic = c(0.1, 10))

  expect_true(all(prior[["tau0"]] < prior[["tau1"]]))
  # Scaling by the LS standard errors makes the bounds coefficient specific.
  expect_gt(length(unique(as.numeric(prior[["tau0"]]))), 1)
  expect_equal(as.numeric(prior[["tau1"]]) / as.numeric(prior[["tau0"]]),
               rep(100, nrow(prior[["tau0"]])))
})

test_that("the semiautomatic SSVS prior says when the sample is too short", {
  # Unlike the Minnesota prior, this least squares estimate is genuinely
  # needed, so a short sample is refused rather than worked around. The fixed
  # 'tau' values remain available.
  expect_error(ssvs_prior(short_var_model(), semiautomatic = c(0.1, 10)),
               "regressors per equation")
  expect_type(ssvs_prior(short_var_model(), tau = c(0.05, 10)), "list")

  short_vec <- create_bvecmodel(stats::window(vec_data(), end = c(1983, 2)),
                                p = 6, r = 1, const = "unrestricted",
                                iterations = 10, burnin = 5)
  expect_error(ssvs_prior(short_vec, semiautomatic = c(0.1, 10)),
               "regressors per equation")
})

test_that("ssvs_prior is refused for model types it cannot handle", {
  sv_model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                               error = "sv", iterations = 10, burnin = 5)
  tvp_model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                                tvp = TRUE, iterations = 10, burnin = 5)

  expect_error(ssvs_prior(sv_model), "stochastic volatility")
  expect_error(ssvs_prior(tvp_model), "time varying parameter")
})

test_that("inclusion_prior returns probabilities and selectable positions", {
  spec <- fx_var_model()[["model"]]
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])
  prior <- inclusion_prior(fx_var_model(), prob = 0.3)

  expect_named(prior, c("prior", "include"))
  expect_identical(dim(prior[["prior"]]), c(n_coeffs, 1L))
  expect_true(all(prior[["prior"]] == 0.3))
  # Deterministic terms are kept out of the selection by default.
  expect_identical(as.integer(prior[["include"]]),
                   seq_len(spec[["k"]] * spec[["k"]] * spec[["p"]]))
})

test_that("inclusion_prior can keep deterministic terms selectable", {
  spec <- fx_var_model()[["model"]]
  n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])
  prior <- inclusion_prior(fx_var_model(), prob = 0.5,
                           exclude_deterministics = FALSE)

  expect_length(prior[["include"]], n_coeffs)
})

test_that("a Minnesota-like inclusion prior varies across coefficients", {
  prior <- inclusion_prior(fx_var_model(), minnesota_like = TRUE)

  expect_true(all(prior[["prior"]] >= 0 & prior[["prior"]] <= 1))
  expect_gt(length(unique(as.numeric(prior[["prior"]]))), 1)
})

test_that("inclusion_prior rejects probabilities outside the unit interval", {
  expect_error(inclusion_prior(fx_var_model(), prob = 1.5), "between 0 and 1")
  expect_error(inclusion_prior(fx_var_model(), prob = -0.1), "between 0 and 1")
})

test_that("variable selection priors are attached by add_priors", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            varsel = "ssvs", iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1),
                      sigma = list(df = 1, scale = 0.0001),
                      varsel = list(tau = c(0.05, 10), inprior = 0.5))

  prior <- model[["priors"]][["a"]]
  expect_true(all(c("inprior", "include", "tau0", "tau1") %in% names(prior)))
  expect_true(all(prior[["inprior"]] == 0.5))
})

test_that("Bayesian variable selection only needs inclusion probabilities", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            varsel = "bvs", iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1),
                      sigma = list(df = 1, scale = 0.0001),
                      varsel = list(inprior = 0.5))

  prior <- model[["priors"]][["a"]]
  expect_true(all(c("inprior", "include") %in% names(prior)))
  expect_false("tau0" %in% names(prior))
})

test_that("a BVS prior of a VEC model takes v_i_det from v_i when it is not given", {
  model <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted",
                            varsel = "bvs", iterations = 10, burnin = 5)
  coint <- list(v_i = 0, p_tau_i = 1)
  sigma <- list(df = "k", scale = 1)
  varsel <- list(inprior = 0.5)

  implicit <- add_priors(model, coef = list(v_i = 1), coint = coint,
                         sigma = sigma, varsel = varsel)
  explicit <- add_priors(model, coef = list(v_i = 1, v_i_det = 1), coint = coint,
                         sigma = sigma, varsel = varsel)
  expect_equal(implicit[["priors"]], explicit[["priors"]])

  expect_warning(
    add_priors(model, coef = list(v_i = 0), coint = coint,
               sigma = sigma, varsel = varsel),
    "uninformative prior"
  )
})

test_that("a VEC model with variable selection has finite inclusion priors and samples", {
  # The loadings are never selected. Their prior inclusion probability used to
  # be NA, which the sampler refuses, so no VEC model with BVS or SSVS could be
  # estimated.
  varsel <- list(bvs = list(inprior = 0.5),
                 ssvs = list(inprior = 0.5, semiautomatic = c(0.1, 10)))
  for (method in names(varsel)) {
    model <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted",
                              varsel = method, iterations = 10, burnin = 5)
    model <- suppressWarnings(
      add_priors(model, coef = list(v_i = 1, v_i_det = 0.1),
                 coint = list(v_i = 0, p_tau_i = 1),
                 sigma = list(df = "k", scale = 1), varsel = varsel[[method]]))
    inprior <- model[["priors"]][["a"]][["inprior"]]
    expect_true(all(is.finite(inprior)), info = method)
    expect_true(all(inprior >= 0 & inprior <= 1), info = method)

    model <- add_initial_values(model)
    set.seed(20260913)
    expect_no_error(model <- add_posterior_coefficients(model))
    expect_identical(n_draws(model), 10L)
  }
})

test_that("BVS can be combined with a Minnesota prior without v_i", {
  minnesota <- list(kappa1 = 0.5, kappa2 = 0.1, kappa4 = 5)
  varsel <- list(inprior = 0.5)

  var_model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                                varsel = "bvs", iterations = 10, burnin = 5)
  expect_no_warning(
    var_model <- add_priors(var_model, coef = list(minnesota = minnesota),
                            sigma = list(df = 1, scale = 0.0001),
                            varsel = varsel)
  )
  expect_true(all(c("v_inv", "inprior", "include") %in%
                    names(var_model[["priors"]][["a"]])))

  vec_model <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted",
                                varsel = "bvs", iterations = 10, burnin = 5)
  expect_no_warning(
    vec_model <- add_priors(vec_model, coef = list(minnesota = minnesota),
                            coint = list(v_i = 0, p_tau_i = 1),
                            sigma = list(df = "k", scale = 1),
                            varsel = varsel)
  )
  expect_true(all(c("v_inv", "inprior", "include") %in%
                    names(vec_model[["priors"]][["a"]])))
})

test_that("varsel$exclude_det keeps deterministic terms out of the selection", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            varsel = "bvs", iterations = 10, burnin = 5)
  spec <- model[["model"]]
  priors <- function(exclude_det) {
    add_priors(model, coef = list(v_i = 1, v_i_det = 1),
               sigma = list(df = 1, scale = 0.0001),
               varsel = list(inprior = 0.5, exclude_det = exclude_det))
  }

  n_lagged <- spec[["k"]] * spec[["k"]] * spec[["p"]]
  n_det <- spec[["k"]] * spec[["n"]]
  expect_length(priors(TRUE)[["priors"]][["a"]][["include"]], n_lagged)
  expect_length(priors(FALSE)[["priors"]][["a"]][["include"]], n_lagged + n_det)
})

test_that("add_priors rejects unknown elements of varsel", {
  var_ssvs <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                               varsel = "ssvs", iterations = 10, burnin = 5)
  expect_error(
    add_priors(var_ssvs, coef = list(v_i = 1, v_i_det = 1),
               sigma = list(df = 1, scale = 0.0001),
               varsel = list(inprior = 0.5, tau = c(0.05, 10),
                             exclude_deterministic = TRUE)),
    "Element 'exclude_deterministic' in argument 'varsel' is not recognised"
  )

  var_bvs <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                              varsel = "bvs", iterations = 10, burnin = 5)
  expect_error(
    add_priors(var_bvs, coef = list(v_i = 1, v_i_det = 1),
               sigma = list(df = 1, scale = 0.0001),
               varsel = list(inprior = 0.5, exclude_deterministic = TRUE)),
    "Element 'exclude_deterministic' in argument 'varsel' is not recognised"
  )

  vec_bvs <- create_bvecmodel(vec_data(), p = 1, r = 1, const = "unrestricted",
                              varsel = "bvs", iterations = 10, burnin = 5)
  expect_error(
    add_priors(vec_bvs, coef = list(v_i = 1, v_i_det = 1 / 10),
               coint = list(v_i = 0, p_tau_i = 1),
               sigma = list(df = "k", scale = 1),
               varsel = list(inprior = 0.5, inprob = 0.5)),
    "Element 'inprob' in argument 'varsel' is not recognised"
  )
})

test_that("a Minnesota-like inclusion prior of a VEC model has one block per exogenous difference", {
  # The exogenous variables of a VEC model enter with their current difference
  # and s - 1 lags of it. The prior was filled with s lags after the current
  # difference, which ran past the matrix whenever the model had fewer
  # unrestricted deterministic terms than exogenous variables.
  endo <- stats::window(at_domestic()[, c("y", "r")], end = c(2019, 4)) * 100
  exo <- stats::window(bvartools::at_macrodata[["foreign"]][, c("y.s", "r.s", "poil")],
                       end = c(2019, 4)) * 100
  kappa <- c(0.8, 0.5, 0.4, 0.9)

  # The expected probabilities, read off the names of the regressors rather than
  # from the positions of their blocks.
  expected_prior <- function(model) {
    k <- model[["model"]][["k"]]
    x_names <- colnames(model[["data"]][["train"]][["x"]])
    endo_names <- model[["model"]][["endogen"]]
    is_diff <- grepl("^d[.].+[.]l[0-9]+$", x_names)
    variable <- sub("^d[.](.+)[.]l[0-9]+$", "\\1", x_names)
    lag <- rep(NA_integer_, length(x_names))
    lag[is_diff] <- as.integer(sub("^.*[.]l", "", x_names[is_diff]))
    prior <- matrix(NA_real_, k, length(x_names))
    for (j in seq_along(x_names)) {
      if (is_diff[j] && variable[j] %in% endo_names) {
        prior[, j] <- ifelse(endo_names == variable[j], kappa[1], kappa[2]) / lag[j]
      } else if (is_diff[j]) {
        prior[, j] <- kappa[3] / (1 + lag[j])
      } else {
        prior[, j] <- kappa[4]
      }
    }
    c(rep(1, k * model[["model"]][["rank"]]), prior)
  }

  for (s in 1:3) {
    for (const in list(NULL, "unrestricted")) {
      model <- create_bvecmodel(endo, p = 2, exogen = exo, s = s, r = 1, const = const,
                                varsel = "bvs", iterations = 10, burnin = 5)
      prior <- inclusion_prior(model, minnesota_like = TRUE, kappa1 = kappa[1],
                               kappa2 = kappa[2], kappa3 = kappa[3], kappa4 = kappa[4])
      expect_equal(as.numeric(prior[["prior"]]), expected_prior(model),
                   info = paste("s =", s, "and const", if (is.null(const)) "absent" else const))
    }
  }

  # The same prior reached through add_priors()
  model <- create_bvecmodel(endo, p = 2, exogen = exo, s = 2, r = 1,
                            varsel = "bvs", iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 1),
                      coint = list(v_i = 0, p_tau_i = 1),
                      sigma = list(df = "k", scale = 1),
                      varsel = list(inprior = 0.5, minnesota = kappa))
  expect_equal(as.numeric(model[["priors"]][["a"]][["inprior"]]), expected_prior(model))
})

test_that("the semiautomatic SSVS prior of a structural model uses least squares of its recursive equations", {
  # A structural VAR and a structural VEC model of at_macrodata. Equation i
  # regresses variable i on the regressors and on minus the current values of
  # the variables before it. The VAR method refused structural models, and the
  # VEC method fell back to 'tau' for their contemporaneous coefficients.
  per_equation_se <- function(y, x) {
    y <- matrix(as.numeric(y), nrow(y))
    x <- matrix(as.numeric(x), nrow(x))
    lapply(seq_len(ncol(y)), function(i) {
      z <- cbind(x, -y[, seq_len(i - 1), drop = FALSE])
      residuals <- y[, i] - z %*% solve(crossprod(z), crossprod(z, y[, i]))
      sqrt(diag(solve(crossprod(z))) * sum(residuals^2) / (nrow(z) - ncol(z)))
    })
  }
  stack_se <- function(se, n_x, k) {
    c(c(t(sapply(se, function(s) s[seq_len(n_x)]))),
      unlist(lapply(1:(k - 1), function(j) sapply((j + 1):k, function(i) se[[i]][n_x + j]))))
  }

  svar <- create_bvarmodel(diff(at_data()) * 100, p = 2, deterministic = "const",
                           structural = TRUE, error = "gamma", varsel = "ssvs",
                           iterations = 10, burnin = 5)
  y <- as.matrix(svar[["data"]][["train"]][["y"]])
  x <- as.matrix(svar[["data"]][["train"]][["x"]])
  prior <- ssvs_prior(svar, semiautomatic = c(0.1, 10))
  expect_equal(as.numeric(prior[["tau1"]]) / 10, stack_se(per_equation_se(y, x), ncol(x), 3))
  expect_no_error(add_priors(svar, coef = list(v_i = 0), sigma = list(shape = 3, rate = 1),
                             varsel = list(inprior = 0.5, semiautomatic = c(0.1, 10))))

  svec <- create_bvecmodel(at_data() * 100, p = 2, r = 1, const = "unrestricted",
                           structural = TRUE, error = "gamma", varsel = "ssvs",
                           iterations = 10, burnin = 5)
  y <- as.matrix(svec[["data"]][["train"]][["y"]])
  x <- as.matrix(svec[["data"]][["train"]][["x"]])
  prior <- ssvs_prior(svec, semiautomatic = c(0.1, 10))
  expect_equal(as.numeric(prior[["tau1"]])[-(1:3)] / 10, stack_se(per_equation_se(y, x), ncol(x), 3))
})

test_that("a Minnesota-like inclusion prior gives the contemporaneous coefficients kappa2", {
  svar <- create_bvarmodel(diff(at_data()) * 100, p = 1, deterministic = "const",
                           structural = TRUE, error = "gamma", iterations = 10, burnin = 5)
  prior <- inclusion_prior(svar, minnesota_like = TRUE, kappa1 = 0.8, kappa2 = 0.4, kappa4 = 0.9)
  expect_equal(as.numeric(utils::tail(prior[["prior"]], 3)), rep(0.4, 3))
  expect_equal(as.numeric(utils::tail(inclusion_prior(svar, prob = 0.3)[["prior"]], 3)), rep(0.3, 3))

  svec <- create_bvecmodel(at_data() * 100, p = 2, r = 1, const = "unrestricted",
                           structural = TRUE, error = "gamma", iterations = 10, burnin = 5)
  prior <- inclusion_prior(svec, minnesota_like = TRUE, kappa1 = 0.8, kappa2 = 0.4, kappa4 = 0.9)
  expect_equal(as.numeric(utils::tail(prior[["prior"]], 3)), rep(0.4, 3))
})
