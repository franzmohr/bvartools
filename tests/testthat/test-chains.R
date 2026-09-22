# Several chains: the same simulation once per chain, pooled, and compared by
# chain_diagnostics().

chains_var <- function(thin = 1, iterations = 60) {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            iterations = iterations, burnin = 20, thin = thin)
  model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                      sigma = list(df = "k", scale = 1))
  add_seed(add_initial_values(model), 4711)
}

test_that("one chain draws what the model always drew", {
  model <- chains_var()
  default <- add_posterior_coefficients(model)
  one <- add_posterior_coefficients(model, chains = 1)
  expect_identical(one[["posterior"]], default[["posterior"]])
  expect_null(one[["model"]][["chains"]])
})

test_that("several chains are pooled one after the other, the first with the model's seed", {
  model <- chains_var()
  single <- add_posterior_coefficients(model)
  three <- add_posterior_coefficients(model, chains = 3)
  a <- three[["posterior"]][["a"]][["coeffs"]]
  n <- nrow(single[["posterior"]][["a"]][["coeffs"]])

  expect_identical(three[["model"]][["chains"]], 3L)
  expect_identical(three[["model"]][["seed"]], model[["model"]][["seed"]])
  expect_identical(nrow(a), 3L * n)
  expect_equal(unname(as.matrix(a)[seq_len(n), ]),
               unname(as.matrix(single[["posterior"]][["a"]][["coeffs"]])))
  # The other chains draw other numbers.
  expect_false(isTRUE(all.equal(unname(as.matrix(a)[n + seq_len(n), ]),
                                unname(as.matrix(a)[seq_len(n), ]))))
  expect_identical(nrow(three[["posterior"]][["u_sigma_inv"]][["coeffs"]]), 3L * n)

  # A seeded run reproduces, and the chains' draws label on through each other.
  expect_identical(add_posterior_coefficients(model, chains = 3)[["posterior"]], three[["posterior"]])
  expect_equal(coda::mcpar(a), c(1, 3 * n, 1))
})

test_that("the thinning interval of a chain carries over to the pooled draws", {
  pooled <- add_posterior_coefficients(chains_var(thin = 2, iterations = 30), chains = 2)
  expect_equal(coda::mcpar(pooled[["posterior"]][["a"]][["coeffs"]]), c(2, 120, 2))
})

test_that("forecasts, log-likelihoods and summaries use every chain", {
  model <- add_posterior_coefficients(chains_var(), chains = 2)
  model <- add_posterior_forecasts(add_forecast_input(model, n_ahead = 2))
  model <- add_posterior_loglik(model)
  expect_identical(nrow(model[["posterior"]][["forecast"]][["forecasts"]]), 120L)
  expect_identical(nrow(model[["posterior"]][["loglik"]]), 120L)
  expect_output(print(summary(model)), "Chains: 2\\. Largest split R-hat")
  expect_false(any(grepl("Chains", capture.output(print(summary(add_posterior_coefficients(chains_var())))))))
})

test_that("chain_diagnostics reports split R-hat near one for chains that agree", {
  model <- add_posterior_coefficients(chains_var(iterations = 400), chains = 4)
  diag <- chain_diagnostics(model)
  expect_named(diag, c("block", "parameter", "rhat", "ess"))
  expect_true(all(c("a$coeffs", "u_sigma_inv$coeffs") %in% diag[["block"]]))
  expect_equal(sum(diag[["block"]] == "a$coeffs"), ncol(model[["posterior"]][["a"]][["coeffs"]]))
  expect_true(all(diag[["rhat"]] < 1.05, na.rm = TRUE))
  expect_true(all(diag[["ess"]] > 0))
})

test_that("chain_diagnostics flags chains that describe different distributions", {
  model <- add_posterior_coefficients(chains_var(iterations = 100), chains = 2)
  a <- model[["posterior"]][["a"]][["coeffs"]]
  # Move the second chain of the first coefficient by ten of its standard deviations.
  shifted <- as.matrix(a)
  shifted[101:200, 1] <- shifted[101:200, 1] + 10 * stats::sd(shifted[1:100, 1])
  model[["posterior"]][["a"]][["coeffs"]] <- coda::mcmc(shifted)
  diag <- chain_diagnostics(model)
  rhat <- diag[diag[["block"]] == "a$coeffs", "rhat"]
  expect_gt(rhat[1], 2)
  expect_true(all(rhat[-1] < 1.2))
  expect_output(print(summary(model)), "The chains disagree")
})

test_that("chain_diagnostics splits a chain that is still drifting", {
  model <- add_posterior_coefficients(chains_var(iterations = 100), chains = 2)
  drift <- as.matrix(model[["posterior"]][["a"]][["coeffs"]])
  trend <- rep(seq(0, 10, length.out = 100), 2) * stats::sd(drift[, 1])
  drift[, 1] <- drift[, 1] + trend
  model[["posterior"]][["a"]][["coeffs"]] <- coda::mcmc(drift)
  diag <- chain_diagnostics(model)
  expect_gt(diag[diag[["block"]] == "a$coeffs", "rhat"][1], 1.5)
})

test_that("chains are refused where they cannot be compared", {
  expect_error(add_posterior_coefficients(chains_var(), chains = 0), "whole number")
  expect_error(add_posterior_coefficients(chains_var(), chains = 1.5), "whole number")
  expect_error(chain_diagnostics(add_posterior_coefficients(chains_var())), "single chain")
  expect_error(chain_diagnostics(list()), "bvarmodel")

  discount <- create_bvarmodel(var_data(), p = 1, deterministic = "const", algorithm = "discount",
                               iterations = 10, burnin = 0)
  discount <- add_initial_values(add_priors(discount, coef = list(v_i = 1, v_i_det = 0.1),
                                            sigma = list(df = "k", scale = 1)))
  expect_error(add_posterior_coefficients(discount, chains = 2), "closed-form")

  # Thinning by a factor the length of a chain is not a multiple of used to mix
  # the chains. Each is thinned on its own now, so they stay comparable.
  model <- add_posterior_coefficients(chains_var(), chains = 2)
  expect_no_error(chain_diagnostics(thin(model, thin = 7)))
})

test_that("a VEC model and a list of models take chains as well", {
  vec <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted",
                          iterations = fx_iterations, burnin = fx_burnin)
  vec <- add_priors(vec, coef = list(v_i = 0, v_i_det = 0), coint = list(v_i = 0, p_tau_i = 1),
                    sigma = list(df = "k", scale = 1))
  vec <- add_posterior_coefficients(add_seed(add_initial_values(vec), 99), chains = 2)
  expect_identical(nrow(vec[["posterior"]][["beta"]][["coeffs"]]), 2L * fx_iterations)
  expect_true("beta$coeffs" %in% chain_diagnostics(vec)[["block"]])

  models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                             iterations = 40, burnin = 10)
  models <- add_priors(models, coef = list(v_i = 0, v_i_det = 0), sigma = list(df = "k", scale = 1))
  models <- add_posterior_coefficients(add_seed(add_initial_values(models), 5), chains = 2)
  expect_true(all(sapply(models, function(m) m[["model"]][["chains"]]) == 2L))
  expect_true(all(sapply(models, function(m) nrow(m[["posterior"]][["a"]][["coeffs"]])) == 80L))
})

test_that("the number of chains travels through a file", {
  model <- add_posterior_coefficients(chains_var(), chains = 2)
  path <- temp_h5_file()
  write_to_hdf5(model, path)
  back <- read_model_from_hdf5(path)
  expect_equal(as.integer(back[["model"]][["chains"]]), 2L)
  expect_equal(chain_diagnostics(back)[["rhat"]], chain_diagnostics(model)[["rhat"]])
})
