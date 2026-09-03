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
