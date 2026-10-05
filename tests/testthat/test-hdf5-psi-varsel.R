skip_if_not_installed("hdf5r")

# BayesTS reads the selection scheme of the covariance block from the attribute
# 'varsel' of /model/priors/psi and from nowhere else. Every model with such a
# block reads it, with constant coefficients as well as time varying ones, so
# `varsel$covar` decides it everywhere rather than only where the sampler
# happened to have somewhere to put it.

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

constant_covar_model <- function(covar) {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            error = "sv+covar", varsel = "bvs",
                            iterations = 10, burnin = 10)
  model <- add_priors(model,
                      coef = list(v_i = 1, v_i_det = 1 / 10),
                      sigma = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                                   state_variance = 0.05, offset = 1e-8),
                      varsel = list(inprior = 0.5, exclude_det = TRUE, covar = covar))
  add_initial_values(model)
}

test_that("a constant coefficient model may select over the coefficients alone", {

  # This was refused until the samplers with constant coefficients read the
  # covariance block's own scheme: they applied the model's to both blocks and
  # had no way of being told otherwise.
  model <- constant_covar_model(covar = FALSE)

  expect_equal(model[["priors"]][["psi"]][["varsel"]], "none")
  expect_equal(model[["model"]][["varsel"]], "bvs")
  # The coefficients keep the prior selection needs; the covariances get none.
  expect_false(is.null(model[["priors"]][["a"]][["inprior"]]))
  expect_null(model[["priors"]][["psi"]][["inprior"]])

  path <- temp_h5_file()
  write_to_hdf5(model, filename = path)
  expect_equal(psi_varsel_attribute(path), "none")
})

test_that("and selecting over both is still what asking for both gives", {
  model <- constant_covar_model(covar = TRUE)

  expect_equal(model[["priors"]][["psi"]][["varsel"]], "bvs")
  expect_equal(nrow(model[["priors"]][["psi"]][["inprior"]]),
               model[["model"]][["k"]] * (model[["model"]][["k"]] - 1) / 2)

  path <- temp_h5_file()
  write_to_hdf5(model, filename = path)
  expect_equal(psi_varsel_attribute(path), "bvs")
})

test_that("the covariances are estimated rather than switched off when unselected", {

  executable <- Sys.getenv("BAYESTS_EXECUTABLE")
  skip_if(!nzchar(executable) || !file.exists(executable),
          "BAYESTS_EXECUTABLE does not name a BayesTS executable")
  library_path <- Sys.getenv("BAYESTS_LIBRARY_PATH")
  library_path <- if (nzchar(library_path)) {
    strsplit(library_path, .Platform$path.sep, fixed = TRUE)[[1]]
  }
  simulate <- bayests_posterior(executable = executable, library_path = library_path)

  # Drawn, not merely written: a covariance left out of the selection has no
  # draw that is exactly zero, where a selected one has many.
  selected <- add_posterior_coefficients(constant_covar_model(covar = TRUE),
                                         posterior_function = simulate)
  unselected <- add_posterior_coefficients(constant_covar_model(covar = FALSE),
                                           posterior_function = simulate)

  # Only the strict lower triangle is drawn. The rest of the vectorised Psi is
  # the identity and its structural zeros, which are zero whatever is selected.
  k <- selected[["model"]][["k"]]
  drawn <- which(lower.tri(matrix(0, k, k)))

  expect_true(any(selected[["posterior"]][["psi"]][["coeffs"]][, drawn] == 0))
  expect_false(any(unselected[["posterior"]][["psi"]][["coeffs"]][, drawn] == 0))
})
