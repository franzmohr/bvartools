# cointspace_prior() builds the prior on the cointegration space for
# add_priors(), and for packages whose VEC models have a layout of their own.
# The priors it gives are tested with add_priors() in test-add_priors.R; what is
# tested here is what the exported function adds to that.

vec_coint_priors <- function(coint, model = fx_vec_model()) {
  add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10), coint = coint,
             sigma = list(df = "k", scale = 1))
}

test_that("add_priors() stores the prior cointspace_prior() builds", {
  model <- fx_vec_model()
  coints <- list(list(v_i = 0, p_tau_i = 1),
                 list(v_i = 0.1, p_tau_i = c(1, 2)),
                 list(v_i = 0.01, p_tau_i = "ml"),
                 list(v_i = "ml", p_tau_i = "ml", weight = 10))
  for (coint in coints) {
    expect_identical(vec_coint_priors(coint, model)[["priors"]][["beta"]],
                     cointspace_prior(model, coint))
  }
})

test_that("a model without cointegration gets no prior, but its coint is checked", {
  model <- fx_vec_model()
  model[["model"]][["rank"]] <- 0L
  expect_null(cointspace_prior(model, list(v_i = 0, p_tau_i = 1)))
  expect_error(cointspace_prior(model, list(v_i = 0)), "coint[$]p_tau_i")
  expect_error(cointspace_prior(model, list(v_i = 0, p_tau_i = 1, tau = 1)), "not recognised")
})

test_that("the prior is sized from the error correction term, not from model$m", {
  # A model of another layout, such as a sub-model of a global VEC model, whose
  # error correction term holds series that its model$m counts differently from
  # a model of create_bvecmodel().
  model <- fx_vec_model()
  k_w <- NCOL(model[["data"]][["train"]][["w"]])
  model[["model"]][["m"]] <- model[["model"]][["m"]] + 3L

  numeric_prior <- cointspace_prior(model, list(v_i = 0.1, p_tau_i = 1))
  expect_identical(dim(numeric_prior[["p_tau_inv"]]), c(k_w, k_w))

  ml_prior <- cointspace_prior(model, list(v_i = 0.01, p_tau_i = "ml"))
  expect_identical(dim(ml_prior[["p_tau_inv"]]), c(k_w, k_w))

  expect_error(cointspace_prior(model, list(v_i = 0.1, p_tau_i = diag(k_w + 3))),
               paste0(k_w, " x ", k_w))
})

test_that("a cointegration rank without an error correction term is refused", {
  model <- fx_vec_model()
  model[["data"]][["train"]][["w"]] <- NULL
  expect_error(cointspace_prior(model, list(v_i = 0, p_tau_i = 1)), "error correction term")
})
