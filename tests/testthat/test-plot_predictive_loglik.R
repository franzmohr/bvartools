# Plotting the terms of the log predictive likelihood. The sum and its band are
# selection_criteria()'s and are tested with the criteria; what is tested here is
# that the terms are read, differenced and accumulated correctly, that the
# refusals are the ones documented, and that "LPL" reaches plot.selcritlist().

# A 'selcrit' carrying an LPL entry built from given per-period densities, which
# is what selection_criteria.expandingwindow() produces. Built by hand so that
# the arithmetic of the plot is checked against numbers the test chose.
lpl_entry <- function(periods, lpd) {
  terms <- data.frame(period = periods, lpd = lpd, nse = rep(0.01, length(lpd)))
  entry <- data.frame(mean = sum(lpd), median = NA_real_,
                      qlower = sum(lpd) - 1, qupper = sum(lpd) + 1)
  attr(entry, "terms") <- terms
  out <- list(LPL = entry)
  class(out) <- append("selcrit", class(out))
  out
}

two_models <- function() {
  x <- list(base = lpl_entry(2020 + (0:3) / 4, c(-1, -1, -1, -1)),
            rival = lpl_entry(2020 + (0:3) / 4, c(-0.5, -0.5, -9, -0.5)))
  class(x) <- append("selcritlist", class(x))
  x
}

test_that("the terms are differenced and accumulated", {

  x <- two_models()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  # Returned invisibly and unchanged, whatever is plotted.
  expect_identical(plot_predictive_loglik(x, baseline = "base"), x)

  # The rival gains 0.5 in each of three periods and loses 8 in the third, so a
  # cumulative difference ends at -6.5 and would read +1.5 after two periods.
  # Recomputed here rather than scraped off the device.
  lpd <- attr(x[["rival"]][["LPL"]], "terms")[["lpd"]]
  base <- attr(x[["base"]][["LPL"]], "terms")[["lpd"]]
  expect_equal(cumsum(lpd - base), c(0.5, 1, -7, -6.5))
  expect_equal(sum(lpd - base), -6.5)
})

test_that("a single model and a missing baseline are both allowed", {

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  one <- lpl_entry(2020 + (0:3) / 4, c(-1, -2, -1, -1))
  expect_s3_class(plot_predictive_loglik(one), "selcrit")
  # No baseline, so the levels are plotted and no differencing is attempted.
  expect_s3_class(plot_predictive_loglik(one, cumulative = FALSE), "selcrit")
  expect_s3_class(plot_predictive_loglik(two_models()), "selcritlist")
})

test_that("the refusals are the documented ones", {

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  x <- two_models()
  expect_error(plot_predictive_loglik(list()), "'selcrit' or 'selcritlist'")
  expect_error(plot_predictive_loglik(x, cumulative = NA), "TRUE or FALSE")
  expect_error(plot_predictive_loglik(x, baseline = "no-such-model"),
               "does not refer to a model")
  expect_error(plot_predictive_loglik(x, baseline = 5), "does not refer to a model")

  # A list whose models were never scored.
  empty <- list(a = structure(list(WAIC = 1), class = c("selcrit", "list")))
  class(empty) <- append("selcritlist", class(empty))
  expect_error(plot_predictive_loglik(empty), "No model in argument")

  # One scored model beside an unscored one: the unscored one is dropped.
  mixed <- list(a = x[["base"]],
                b = structure(list(WAIC = 1), class = c("selcrit", "list")))
  class(mixed) <- append("selcritlist", class(mixed))
  expect_warning(plot_predictive_loglik(mixed), "not plotted")
})

test_that("plot.selcritlist() accepts LPL as a criterion", {

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  # Through the real path, because plot.selcritlist() also reads the model
  # specifications, which a hand-built entry does not carry.
  windows <- add_predictive_loglik(fx_expanding_window())
  one <- selection_criteria(windows)
  expect_false(is.null(one[["LPL"]]))
  expect_false(is.null(attr(one[["LPL"]], "terms")))

  models <- list(first = one, second = one)
  class(models) <- append("selcritlist", class(models))

  # The criterion is accepted where it used to be reported as unavailable.
  expect_invisible(plot(models, criterion = "LPL"))
  expect_error(plot(models, criterion = "no-such-criterion"),
               "not contained in argument")

  # And the terms of that same object plot.
  expect_s3_class(plot_predictive_loglik(models, baseline = 1), "selcritlist")
})
