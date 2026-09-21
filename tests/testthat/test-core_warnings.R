# What the vendored core reports about a model it will run but cannot run well,
# and the priors it refuses outright. The core cannot raise an R condition
# itself; RcppReporter collects its warnings and add_posterior_coefficients()
# raises them once the sampler has returned.

core_model <- function(varsel, coef, vs) {
  m <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                        iterations = fx_iterations, burnin = fx_burnin, varsel = varsel)
  m <- add_priors(m, coef = coef, sigma = list(df = 1, scale = 0.0001), varsel = vs)
  add_initial_values(m)
}

test_that("bvs against a near-flat prior warns, once per block", {
  m <- core_model("bvs", list(v_i = 0.001, v_i_det = 0.001), list(inprior = 0.5))
  w <- NULL
  withCallingHandlers(add_posterior_coefficients(m), warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  expect_length(w, 1)
  expect_match(w, "too flat to select against", fixed = TRUE)
})

test_that("bvs against a prior it can select against says nothing", {
  m <- core_model("bvs", list(v_i = 1, v_i_det = 1), list(inprior = 0.5))
  expect_no_warning(add_posterior_coefficients(m))
})

test_that("the warning does not leave a trace on the model", {
  m <- core_model("bvs", list(v_i = 0.001, v_i_det = 0.001), list(inprior = 0.5))
  fitted <- suppressWarnings(add_posterior_coefficients(m))
  expect_null(fitted[["warnings"]])
})

test_that("ssvs refuses a non-zero prior mean at a selected position", {
  m <- core_model("ssvs", list(v_i = 1, const = "mean"),
                  list(inprior = 0.5, tau = c(0.05, 10)))
  expect_error(add_posterior_coefficients(m), "must be zero at every selected position")
})

test_that("ssvs runs once that position is left out of the selection", {
  m <- core_model("ssvs", list(v_i = 1, const = "mean"),
                  list(inprior = 0.5, tau = c(0.05, 10), exclude_det = TRUE))
  expect_s3_class(add_posterior_coefficients(m), "bvarmodel")
})
