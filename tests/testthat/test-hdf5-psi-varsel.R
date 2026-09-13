skip_if_not_installed("hdf5r")

# BayesTS reads the selection scheme of the covariance block of a time varying
# model from the attribute 'varsel' of /model/priors/psi and from nowhere else.

tvp_covar_model <- function(covar) {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            tvp = TRUE, error = "sv+covar", varsel = "bvs",
                            iterations = 10, burnin = 10)
  model <- add_priors(model,
                      coef = list(v_i = 1, v_i_det = 1 / 10, shape = 3, rate = 1e-4),
                      sigma = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                                   state_variance = 0.05, offset = 1e-8),
                      varsel = list(inprior = 0.5, exclude_det = TRUE, covar = covar))
  add_initial_values(model)
}

psi_varsel_attribute <- function(path) {
  h5 <- hdf5r::h5file(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  if (!h5$exists("model") || !h5[["model"]]$exists("priors") ||
      !h5[["model/priors"]]$exists("psi")) {
    return(NULL)
  }
  hdf5r::h5attr(h5[["model/priors/psi"]], "varsel")
}

test_that("the selection scheme of the covariance block is written where BayesTS reads it", {
  model <- tvp_covar_model(covar = TRUE)
  expect_equal(model[["priors"]][["psi"]][["varsel"]], "bvs")

  path <- temp_h5_file()
  write_to_hdf5(model, filename = path)

  expect_equal(psi_varsel_attribute(path), "bvs")
})

test_that("a covariance block without selection is written as such", {
  model <- tvp_covar_model(covar = FALSE)

  path <- temp_h5_file()
  write_to_hdf5(model, filename = path)

  expect_equal(psi_varsel_attribute(path), model[["priors"]][["psi"]][["varsel"]])
})

test_that("the attribute does not change what is read back", {
  model <- tvp_covar_model(covar = TRUE)

  path <- temp_h5_file()
  write_to_hdf5(model, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_equal(restored[["priors"]][["psi"]][["varsel"]], "bvs")
  expect_null(restored[["model"]][["priors"]])
})
