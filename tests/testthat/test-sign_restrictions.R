# Sign restrictions.
#
# The property that matters is not the shape of the output but that every draw
# the identification kept actually satisfies what was asked of it, and that the
# rotation it was kept for leaves the model observationally equivalent to the
# draw it came from. Most of what follows checks one of those two.

# The responses of every identified draw, one row per draw.
sign_draws <- function(object, impulse, response, n_ahead) {
  irf(object, impulse = impulse, response = response, n_ahead = n_ahead,
      type = "sign", keep_draws = TRUE)
}

test_that("the accepted rotations are orthogonal", {
  object <- fx_var_sign()
  k <- object[["model"]][["k"]]
  rotations <- object[["posterior"]][["q"]][["coeffs"]]

  expect_false(is.null(rotations))
  expect_identical(ncol(rotations), k * k)

  for (i in seq_len(nrow(rotations))) {
    q <- matrix(rotations[i, ], k)
    expect_equal(q %*% t(q), diag(1, k))
  }
})

test_that("a rotated draw describes the data as well as the draw it came from", {
  object <- fx_var_sign()
  k <- object[["model"]][["k"]]
  rotations <- object[["posterior"]][["q"]][["coeffs"]]
  draws <- bvartools:::.collect_draws(object)

  # P Q (P Q)' = Sigma is what makes the rotated model observationally
  # equivalent, and so what makes the identification set valued rather than a
  # second, better estimate.
  for (i in seq_len(nrow(rotations))) {
    impact <- t(chol(draws[[i]][["Sigma"]])) %*% matrix(rotations[i, ], k)
    expect_equal(impact %*% t(impact), draws[[i]][["Sigma"]])
  }
})

test_that("every identified draw satisfies the restrictions", {
  object <- fx_var_sign()

  for (response in c("invest", "cons")) {
    impact <- sign_draws(object, "invest", response, n_ahead = 1)[, 1]
    expect_true(all(impact > 0))
  }
})

test_that("a negative restriction is imposed as asked", {
  restrictions <- data.frame(impulse = "income",
                             response = c("income", "cons"),
                             sign = c(1, -1),
                             horizon = 0)

  set.seed(515)
  object <- add_sign_restrictions(fx_var_fitted(), restrictions)

  expect_true(all(sign_draws(object, "income", "income", 1)[, 1] > 0))
  expect_true(all(sign_draws(object, "income", "cons", 1)[, 1] < 0))
})

test_that("a restriction away from the impact period is imposed at its horizon", {
  restrictions <- data.frame(impulse = "cons",
                             response = "invest",
                             sign = -1,
                             horizon = 2)

  set.seed(616)
  object <- add_sign_restrictions(fx_var_fitted(), restrictions)

  responses <- sign_draws(object, "cons", "invest", n_ahead = 3)

  # The third column is horizon 2, which is the one that was restricted.
  expect_true(all(responses[, 3] < 0))
})

test_that("the horizon column defaults to the impact period", {
  with_horizon <- data.frame(impulse = "invest",
                             response = "cons", sign = 1, horizon = 0)
  without <- data.frame(impulse = "invest", response = "cons", sign = 1)

  set.seed(717)
  a <- add_sign_restrictions(fx_var_fitted(), with_horizon)
  set.seed(717)
  b <- add_sign_restrictions(fx_var_fitted(), without)

  expect_equal(a[["posterior"]][["q"]][["coeffs"]],
               b[["posterior"]][["q"]][["coeffs"]])
})

test_that("the identification is reproducible from a seed", {
  set.seed(818)
  a <- add_sign_restrictions(fx_var_fitted(), fx_sign_restrictions())
  set.seed(818)
  b <- add_sign_restrictions(fx_var_fitted(), fx_sign_restrictions())

  expect_equal(a[["posterior"]][["q"]][["coeffs"]],
               b[["posterior"]][["q"]][["coeffs"]])
})

test_that("a variance decomposition of an identified model adds up", {
  decomp <- fevd(fx_var_sign(), response = "cons", n_ahead = 4, type = "sign")

  # A rotation of the Choleski factor still factorises Sigma, so the shares
  # remain shares.
  expect_equal(as.numeric(rowSums(decomp)), rep(1, nrow(decomp)))
})

test_that("spillover measures of an identified model keep their identity", {
  measures <- spillover(fx_var_sign(), n_ahead = 5, type = "sign",
                        keep_draws = TRUE)

  expect_s3_class(measures, "bvarspillover")

  # The accounting identity holds draw by draw, which is where it is defined;
  # the medians of the summarised measures are medians of different draws and
  # need not add up.
  expect_equal(rowSums(as.matrix(measures[["from"]])),
               as.numeric(measures[["total"]]))
  expect_equal(rowSums(as.matrix(measures[["to"]])),
               as.numeric(measures[["total"]]))
})

test_that("a draw that no rotation was found for is dropped, not counted", {
  # Six restrictions across two shocks, with too few tries to satisfy them
  # every time. The point is the mixture: some draws identified, some not.
  demanding <- data.frame(
    impulse = c(rep("invest", 3), rep("cons", 3)),
    response = c("invest", "income", "cons", "cons", "invest", "income"),
    sign = c(1, -1, -1, 1, -1, -1),
    horizon = 0
  )

  set.seed(20240101)
  object <- add_sign_restrictions(fx_var_fitted(), demanding, max_tries = 10)

  rotations <- object[["posterior"]][["q"]][["coeffs"]]
  identified <- !is.na(rotations[, 1])

  expect_true(any(identified))
  expect_true(any(!identified))

  # A row is either wholly there or wholly absent; a partly filled rotation
  # would be read as a matrix of numbers and quietly used.
  expect_equal(rowSums(is.na(rotations)), ifelse(identified, 0, ncol(rotations)))

  # The responses cover the identified draws and no others.
  responses <- sign_draws(object, "invest", "income", n_ahead = 2)
  expect_identical(nrow(responses), sum(identified))
  expect_true(all(responses[, 1] < 0))

  # And the summary says how much of the posterior that is.
  reported <- summary(object)[["model"]][["sign_restrictions"]]
  expect_identical(reported[["accepted"]], sum(identified))
  expect_identical(reported[["draws"]], nrow(rotations))
})

test_that("a search that finds nothing says so instead of returning an empty model", {
  # Nine restrictions across all three shocks and a single try per draw. This
  # is not a pattern the model cannot produce -- it is one it produces too
  # rarely to stumble on, which is the case the message is written for. The
  # seed is what makes the outcome of the search fixed.
  demanding <- data.frame(
    impulse = rep(c("invest", "income", "cons"), each = 3),
    response = rep(c("invest", "income", "cons"), times = 3),
    sign = c(1, -1, -1, -1, 1, -1, -1, -1, 1),
    horizon = 0
  )

  set.seed(919)
  expect_error(add_sign_restrictions(fx_var_fitted(), demanding, max_tries = 1),
               "No rotation satisfying the restrictions")
})

test_that("thinning carries the rotations along with the draws", {
  object <- thin(fx_var_sign(), thin = 2)

  expect_identical(nrow(object[["posterior"]][["q"]][["coeffs"]]),
                   nrow(object[["posterior"]][["u_sigma_inv"]][["coeffs"]]))
  expect_no_error(irf(object, impulse = "invest", response = "cons",
                      type = "sign"))
})

test_that("the rotations survive a round trip through HDF5", {
  path <- temp_h5_file()
  on.exit(unlink(path), add = TRUE)

  write_to_hdf5(fx_var_sign(), path)
  restored <- read_model_from_hdf5(path)

  expect_equal(restored[["posterior"]][["q"]][["coeffs"]],
               fx_var_sign()[["posterior"]][["q"]][["coeffs"]])
})

test_that("the restriction table is checked before any rotation is drawn", {
  object <- fx_var_fitted()

  expect_error(add_sign_restrictions(object, list(impulse = "invest")),
               "must be a data frame")
  expect_error(add_sign_restrictions(object, fx_sign_restrictions()[0, ]),
               "does not contain any restriction")
  expect_error(add_sign_restrictions(object,
                                     data.frame(impulse = "invest", sign = 1)),
               "must contain the column 'response'")
  expect_error(add_sign_restrictions(object,
                                     data.frame(impulse = "gdp", response = "cons",
                                                sign = 1)),
               "names a variable 'gdp'")
  expect_error(add_sign_restrictions(object,
                                     data.frame(impulse = "invest", response = "cons",
                                                sign = 2)),
               "must be either 1 or -1")
  expect_error(add_sign_restrictions(object,
                                     data.frame(impulse = "invest", response = "cons",
                                                sign = 1, horizon = -1)),
               "non-negative integers")
  expect_error(add_sign_restrictions(object, fx_sign_restrictions(), max_tries = 0),
               "at least 1")

  # Two rows asking opposite things of the same response would spend the whole
  # budget of tries on every draw before failing.
  contradictory <- data.frame(impulse = "invest", response = c("cons", "cons"),
                              sign = c(1, -1), horizon = 0)
  expect_error(add_sign_restrictions(object, contradictory),
               "both signs of the same response")
})

test_that("models that have nothing to rotate are refused", {
  expect_error(add_sign_restrictions(fx_svar_fitted(), fx_sign_restrictions()),
               "not defined for a structural model")

  no_covar <- fx_var_fitted()
  no_covar[["model"]][["error"]] <- "gamma"
  expect_error(add_sign_restrictions(no_covar, fx_sign_restrictions()),
               "error covariances are estimated")
})

test_that("an unidentified model is refused by the analysis functions", {
  object <- fx_var_fitted()

  expect_error(irf(object, impulse = "invest", response = "cons", type = "sign"),
               "add_sign_restrictions")
  expect_error(fevd(object, response = "cons", type = "sign"),
               "add_sign_restrictions")
  expect_error(spillover(object, type = "sign"),
               "add_sign_restrictions")
})
