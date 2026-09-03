# Plot methods are checked for running cleanly on a null device rather than for
# their appearance.

test_that("posterior draws can be plotted as histograms and traces", {
  expect_plots(plot(fx_var_fitted()))
  expect_plots(plot(fx_var_fitted(), type = "trace"))
})

test_that("VEC posterior draws can be plotted", {
  expect_plots(plot(fx_vec_fitted()))
  expect_plots(plot(fx_vec_fitted(), type = "trace"))
})

test_that("impulse responses and decompositions can be plotted", {
  response <- irf(fx_var_fitted(), impulse = "income", response = "cons",
                  n_ahead = 5)
  decomposition <- fevd(fx_var_fitted(), response = "cons", n_ahead = 5)

  expect_plots(plot(response))
  expect_plots(plot(response, main = "Response", xlab = "Period"))
  expect_plots(plot(decomposition))
})

test_that("forecasts can be plotted", {
  expect_plots(plot(stats::predict(fx_var_forecast(), n_ahead = 5)))
})

test_that("a modellist can be plotted", {
  expect_plots(plot(fx_var_modellist()))
})

test_that("selection criteria can be plotted", {
  expect_plots(plot(selection_criteria(fx_var_modellist())))
})
