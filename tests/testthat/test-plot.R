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

test_that("a decomposition can be plotted with a limited number of groups", {
  decomposition <- fevd(fx_var_fitted(), response = "cons", n_ahead = 5)

  expect_plots(plot(decomposition, max_groups = 2))
  expect_plots(plot(decomposition, max_groups = ncol(decomposition)))
  expect_plots(plot(decomposition, max_groups = 2, main = "Decomposition"))
})

test_that("the plotted groups are the largest contributions plus 'Other'", {
  decomposition <- fevd(fx_var_fitted(), response = "cons", n_ahead = 5)
  largest <- colnames(decomposition)[which.max(colSums(decomposition))]

  # The bars are drawn from the same matrix the legend is labelled with, so the
  # pooling is checked on that matrix rather than on the device.
  pooled <- bvartools:::.limit_fevd_groups(unclass(decomposition), 2)

  expect_identical(colnames(pooled), c(largest, "Other"))
  expect_equal(as.numeric(rowSums(pooled)),
               as.numeric(rowSums(decomposition)))
})

test_that("an implausible number of groups is rejected by the plot method", {
  decomposition <- fevd(fx_var_fitted(), response = "cons", n_ahead = 5)

  expect_error(plot(decomposition, max_groups = 0), "positive integer")
  expect_error(plot(decomposition, max_groups = "two"), "positive integer")
})
