test_that("irf returns a time series of credible bands", {
  response <- irf(fx_var_fitted(), impulse = "Dp", response = "r",
                  n_ahead = 6)

  expect_s3_class(response, "bvarirf")
  expect_s3_class(response, "ts")
  # The horizon includes the impact period.
  expect_identical(nrow(response), 7L)
  expect_identical(colnames(response), c("2.5%", "50%", "97.5%"))
  expect_equal(stats::start(response)[1], 0)
})

test_that("the credible band is ordered", {
  response <- irf(fx_var_fitted(), impulse = "Dp", response = "r",
                  n_ahead = 6)

  expect_true(all(response[, "2.5%"] <= response[, "50%"]))
  expect_true(all(response[, "50%"] <= response[, "97.5%"]))
})

test_that("a forecast error impulse response has no cross-variable impact", {
  cross <- irf(fx_var_fitted(), impulse = "Dp", response = "r",
               n_ahead = 4)
  own <- irf(fx_var_fitted(), impulse = "r", response = "r", n_ahead = 4)

  # Without orthogonalisation the impact matrix is the identity, so a shock to
  # one variable does not move another on impact, and moves itself by one.
  expect_equal(as.numeric(cross[1, ]), rep(0, 3))
  expect_equal(as.numeric(own[1, ]), rep(1, 3))
})

test_that("the median impulse response matches the recursion by hand", {
  model <- fx_var_fitted()
  k <- model[["model"]][["k"]]
  n_ahead <- 4
  draws <- model[["posterior"]][["a"]][["coeffs"]]

  # Reproduce the forecast error impulse response of r to Dp for every
  # draw and compare the median.
  by_draw <- vapply(seq_len(nrow(draws)), function(i) {
    a <- matrix(draws[i, seq_len(k * k)], k)
    phi <- diag(1, k)
    out <- numeric(n_ahead + 1)
    out[1] <- phi[3, 2]
    for (h in seq_len(n_ahead)) {
      phi <- a %*% phi
      out[h + 1] <- phi[3, 2]
    }
    out
  }, numeric(n_ahead + 1))

  expected <- apply(by_draw, 1, stats::median)
  response <- irf(model, impulse = "Dp", response = "r",
                  n_ahead = n_ahead)
  expect_equal(as.numeric(response[, "50%"]), expected)
})

test_that("an orthogonalised response has a non-zero impact", {
  oir <- irf(fx_var_fitted(), impulse = "Dp", response = "r",
             n_ahead = 4, type = "oir")

  expect_s3_class(oir, "bvarirf")
  # The Cholesky factor is lower triangular, so a shock to the second variable
  # moves the third one on impact.
  expect_false(oir[1, "50%"] == 0)
})

test_that("the reduced form impulse response types are all available", {
  for (type in c("feir", "oir", "gir")) {
    response <- irf(fx_var_fitted(), impulse = "Dp", response = "r",
                    n_ahead = 3, type = type)
    expect_s3_class(response, "bvarirf")
    expect_true(all(is.finite(response)))
  }
})

test_that("the credible interval width is controlled by ci", {
  narrow <- irf(fx_var_fitted(), impulse = "Dp", response = "r",
                n_ahead = 4, ci = 0.5)
  wide <- irf(fx_var_fitted(), impulse = "Dp", response = "r",
              n_ahead = 4, ci = 0.95)

  expect_identical(colnames(narrow), c("25%", "50%", "75%"))
  expect_identical(colnames(wide), c("2.5%", "50%", "97.5%"))
  # A narrower interval is contained in the wider one.
  expect_true(all(narrow[, "25%"] >= wide[, "2.5%"]))
  expect_true(all(narrow[, "75%"] <= wide[, "97.5%"]))
  expect_equal(as.numeric(narrow[, "50%"]), as.numeric(wide[, "50%"]))
})

test_that("structural impulse responses need a structural model", {
  for (type in c("sir", "sgir")) {
    expect_error(irf(fx_var_fitted(), impulse = "Dp", response = "r",
                     n_ahead = 3, type = type),
                 "structural model")
  }
})

test_that("a structural model rejects the reduced-form types", {
  for (type in c("feir", "oir", "gir")) {
    expect_error(irf(fx_svar_fitted(), impulse = "Dp", response = "r",
                     n_ahead = 3, type = type),
                 "not defined for a structural model")
  }
})

test_that("a structural model still produces structural responses", {
  for (type in c("sir", "sgir")) {
    response <- irf(fx_svar_fitted(), impulse = "Dp", response = "r",
                    n_ahead = 3, type = type)

    expect_s3_class(response, "bvarirf")
    expect_identical(nrow(response), 4L)
    expect_true(all(is.finite(response)))
  }
})

test_that("unknown variable names are rejected", {
  expect_error(irf(fx_var_fitted(), impulse = "nonexistent", response = "r"))
  expect_error(irf(fx_var_fitted(), impulse = "Dp", response = "nope"))
})

test_that("a converted VEC model can be used for impulse responses", {
  response <- irf(vec_to_var(fx_vec_fitted()), impulse = "lr", response = "Dp",
                  n_ahead = 4)

  expect_s3_class(response, "bvarirf")
  expect_identical(nrow(response), 5L)
  expect_true(all(is.finite(response)))
})

test_that("an impulse response at horizon zero is the impact period alone", {
  model <- fx_var_fitted()

  cross <- irf(model, impulse = "Dp", response = "r", n_ahead = 0)
  own <- irf(model, impulse = "r", response = "r", n_ahead = 0)

  expect_s3_class(cross, "bvarirf")
  expect_s3_class(cross, "ts")
  expect_identical(nrow(cross), 1L)
  expect_identical(colnames(cross), c("2.5%", "50%", "97.5%"))
  # Phi_0 is the identity and there is no recursion to run, so the forecast
  # error response on impact is one for the own response and zero across
  # variables.
  expect_equal(as.numeric(cross[1, ]), rep(0, 3))
  expect_equal(as.numeric(own[1, ]), rep(1, 3))
})

test_that("the horizon zero response is the impact row of a longer horizon", {
  model <- fx_var_fitted()

  impact <- irf(model, impulse = "Dp", response = "r", n_ahead = 0,
                type = "oir")
  longer <- irf(model, impulse = "Dp", response = "r", n_ahead = 4,
                type = "oir")

  expect_equal(as.numeric(impact[1, ]), as.numeric(longer[1, ]))
})

test_that("an orthogonalised response at horizon zero is the Choleski factor", {
  model <- fx_var_fitted()
  k <- model[["model"]][["k"]]
  sigma_draws <- model[["posterior"]][["u_sigma_inv"]][["coeffs"]]

  # The impact matrix of an orthogonalised response is the lower Choleski
  # factor with its columns scaled to unit shocks, so the response of r to
  # Dp on impact is L[3, 2] / L[2, 2] of that draw.
  by_draw <- vapply(seq_len(nrow(sigma_draws)), function(i) {
    l <- t(chol(solve(matrix(sigma_draws[i, ], k))))
    l[3, 2] / l[2, 2]
  }, numeric(1))

  oir <- irf(model, impulse = "Dp", response = "r", n_ahead = 0,
             type = "oir")

  expect_identical(nrow(oir), 1L)
  expect_equal(as.numeric(oir[1, "50%"]), stats::median(by_draw))
})

test_that("a cumulative response at horizon zero keeps one column of draws", {
  model <- fx_var_fitted()
  store <- nrow(model[["posterior"]][["u_sigma_inv"]][["coeffs"]])

  draws <- irf(model, impulse = "r", response = "r", n_ahead = 0,
               cumulative = TRUE, keep_draws = TRUE)

  # Nothing accumulates over a single period, so this is the impact response
  # itself -- in the draws x horizon shape and not its transpose.
  expect_identical(dim(draws), c(store, 1L))
  expect_equal(as.numeric(draws), rep(1, store))
})

test_that("a negative horizon is rejected", {
  expect_error(irf(fx_var_fitted(), impulse = "Dp", response = "r",
                   n_ahead = -1),
               "at least 0")
})
