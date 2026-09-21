# The methods that apply a step to every model of a 'modellist' or every window
# of an 'expandingwindow'. Most of them are an lapply over the single-model
# method, which is tested on its own; what is checked here is that the
# collection comes back a collection, with every member treated, rather than a
# plain list or a list with holes.

test_that("a modellist is forecast and predicted model by model", {

  models <- add_forecast_input(fx_var_modellist(), n_ahead = 3)
  expect_s3_class(models, "modellist")
  expect_true(all(vapply(models, function(m) m[["model"]][["h"]] == 3, logical(1))))

  set.seed(1)
  models <- add_posterior_forecasts(models)
  expect_s3_class(models, "modellist")
  k <- models[[1]][["model"]][["k"]]
  for (m in models) {
    expect_equal(ncol(m[["posterior"]][["forecast"]][["forecasts"]]), 3 * k)
  }

  pred <- predict(models)
  expect_s3_class(pred, "bvarprdlist")
  expect_length(pred, length(models))
  expect_plots(expect_identical(plot(pred), pred))
})


test_that("priors, specifications and predictions reach every window", {

  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            iterations = 10, burnin = 5)
  windows <- use_expanding_window(model, start = c(1997, 2))

  # Priors after the windows rather than before, which is the other order the
  # workflow allows.
  windows <- add_priors(windows, coef = list(v_i = 0, v_i_det = 0),
                        sigma = list(df = 1, scale = 0.0001))
  expect_s3_class(windows, "expandingwindow")
  expect_true(all(vapply(windows, function(w) !is.null(w[["priors"]]), logical(1))))

  specs <- get_model_specifications(windows)
  expect_s3_class(specs, "data.frame")
  expect_equal(nrow(specs), 1L)

  pred <- predict(fx_expanding_forecast())
  expect_s3_class(pred, "expandwindbvarprdlist")
  expect_length(pred, length(fx_expanding_forecast()))
  expect_plots(expect_identical(plot(pred), pred))
})


test_that("the Minnesota prior moments are computed for every model of a list", {

  models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                             iterations = 10, burnin = 5)
  priors <- minnesota_prior(models, kappa1 = 0.5, kappa2 = 0.1)
  expect_length(priors, 2L)
  for (i in seq_along(models)) {
    spec <- models[[i]][["model"]]
    n_coeffs <- spec[["k"]] * (spec[["k"]] * spec[["p"]] + spec[["n"]])
    expect_named(priors[[i]], c("mu", "v_inv", "sigma_inv"))
    expect_identical(dim(priors[[i]][["v_inv"]]), c(n_coeffs, n_coeffs))
  }
})


test_that("sign restrictions and spillovers map over lists and windows", {

  restrictions <- fx_sign_restrictions()

  models <- fx_var_modellist()
  set.seed(2)
  signed <- add_sign_restrictions(models, restrictions)
  expect_s3_class(signed, "modellist")
  for (m in signed) {
    expect_false(is.null(m[["posterior"]][["q"]]))
  }

  spill <- spillover(models, n_ahead = 2)
  expect_length(spill, length(models))
  expect_true(all(vapply(spill, inherits, logical(1), "bvarspillover")))

  windows <- fx_expanding_window()
  set.seed(3)
  signed_windows <- add_sign_restrictions(windows, restrictions)
  expect_s3_class(signed_windows, "expandingwindow")

  path <- spillover(windows, n_ahead = 2)
  expect_s3_class(path, "bvarspilloverts")
  expect_equal(nrow(path), length(windows))
  expect_equal(colnames(path), c("lower", "median", "upper"))
  expect_true(all(path[, "lower"] <= path[, "upper"], na.rm = TRUE))
})


test_that("forecast errors are plotted period by period and by horizon", {

  windows <- fx_expanding_forecast()
  for (criterion in c("FE", "AFE", "RSFE")) {
    expect_plots(expect_identical(plot_forecast_errors_by_period(windows, criterion = criterion),
                                  windows))
  }
  expect_error(plot_forecast_errors_by_period(windows, criterion = "MSE"), "criterion")

  sc <- list(selection_criteria(windows), selection_criteria(windows))
  class(sc) <- c("selcritlist", "list")
  for (criterion in c("FE", "AFE", "RSFE")) {
    expect_plots(expect_identical(plot(sc, criterion = criterion), sc))
  }
})


test_that("RSFE's bands are AFE's up to interpolation, and never below them", {

  sc <- selection_criteria(fx_expanding_forecast())
  cols <- c("median", "qlower", "qupper")
  rsfe <- as.matrix(sc[["RSFE"]][, cols])
  afe <- as.matrix(sc[["AFE"]][, cols])
  # The square root does not reorder the draws, so the two pick the same
  # neighbours; interpolating between squares and rooting cannot come out below
  # interpolating between the roots.
  expect_true(all(rsfe >= afe - 1e-12))
  expect_equal(rsfe, afe, tolerance = 0.1)
  # The root of the mean square is at least the mean absolute value.
  expect_true(all(sc[["RSFE"]][, "mean"] >= sc[["AFE"]][, "mean"] - 1e-12))
})


test_that("the list writers return the paths they wrote", {

  folder <- temp_model_dir()
  paths <- write_to_hdf5(fx_var_modellist(), folder = folder)
  expect_length(paths, length(fx_var_modellist()))
  expect_true(all(file.exists(paths)))

  folder <- temp_model_dir()
  paths <- write_to_hdf5(fx_expanding_window(), folder = folder)
  expect_length(paths, length(fx_expanding_window()))
  expect_true(all(file.exists(paths)))
  expect_true(all(grepl("Window-", basename(paths))))
})
