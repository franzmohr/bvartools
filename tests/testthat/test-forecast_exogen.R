test_that("forecast input needs the exogenous values the lagged regressors reach back to", {
  y <- var_data()
  foreign <- bvartools::at_macrodata[["foreign"]]
  exo <- stats::ts.intersect(dys = diff(foreign[, "y.s"]),
                             dpoil = diff(foreign[, "poil"])) * 100
  model <- create_bvarmodel(y, p = 1, exogen = stats::window(exo, end = stats::end(y)), s = 1,
                            deterministic = "const", iterations = 10, burnin = 5)
  k <- model[["model"]][["k"]]
  freq <- stats::frequency(y)
  last <- stats::tsp(model[["data"]][["train"]][["y"]])[2]
  future <- stats::window(exo, start = last + 1 / freq, end = last + 4 / freq)

  # The forecast periods alone lack the last in-sample period, whose values are
  # the lagged regressors of the first forecast period. This used to fail with
  # R's "number of items to replace is not a multiple of replacement length".
  expect_error(add_forecast_input(model, n_ahead = 4, exogen = future),
               "must cover the periods from")

  # Reaching back one period is enough, and the regressors are then those values:
  # the current ones first, the lagged ones after them.
  covered <- stats::window(exo, start = last, end = last + 4 / freq)
  x <- add_forecast_input(model, n_ahead = 4, exogen = covered)[["data"]][["forecast"]][["x"]]
  expect_equal(x[, k + 1:2], matrix(as.numeric(future), ncol = 2))
  expect_equal(x[1, k + 3:4], as.numeric(stats::window(exo, start = last, end = last)))

  # A series that ends before the last forecast period is refused as well.
  expect_error(add_forecast_input(model, n_ahead = 4,
                                  exogen = stats::window(exo, start = last, end = last + 2 / freq)),
               "must cover the periods from")
})
