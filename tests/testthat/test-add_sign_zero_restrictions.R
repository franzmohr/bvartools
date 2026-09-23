# Identification by sign and zero restrictions.
#
# The algorithm itself is tested in test-arias_rubio_ramirez_waggoner.R. What
# is checked here is the R side of it: that the table of restrictions is read
# the way the documentation says, that every draw the function returns
# satisfies what was asked of it, and that the resample leaves an object the
# rest of the package can go on working with.

szr_model <- function() {
  cached_fixture("szr_model", {
    object <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                               iterations = 60, burnin = 10)
    object <- add_priors(object, coef = list(v_i = 1, v_i_det = 0.1),
                         sigma = list(df = "k", scale = 1))
    object <- add_initial_values(object)
    add_posterior_coefficients(add_seed(object, 216))
  })
}

# The first variable of the fixture is output growth, whose shock is the only
# one with room for a zero restriction on top of it.
szr_restrictions <- function() {
  data.frame(impulse = "y", response = c("Dp", "r"), sign = c(0, 1), horizon = 0)
}

# The impact responses of one identified draw.
szr_impact <- function(object, i) {
  draws <- bvartools:::.collect_draws(object)
  k <- object[["model"]][["k"]]
  q <- matrix(object[["posterior"]][["q"]][["coeffs"]][i, ], k)
  t(solve(solve(chol(draws[[i]][["Sigma"]]), q)))
}

test_that("every draw that comes back satisfies the restrictions", {
  set.seed(1234)
  object <- add_sign_zero_restrictions(szr_model(), szr_restrictions())

  rotations <- object[["posterior"]][["q"]][["coeffs"]]
  expect_false(anyNA(rotations))

  k <- object[["model"]][["k"]]
  for (i in seq_len(nrow(rotations))) {
    q <- matrix(rotations[i, ], k)
    expect_equal(t(q) %*% q, diag(1, k))

    impact <- szr_impact(object, i)
    expect_equal(impact[2, 1], 0)   # inflation does not respond on impact
    expect_gt(impact[3, 1], 0)      # the interest rate rises
  }
})

test_that("the draws are resampled and the object still works", {
  set.seed(1234)
  object <- add_sign_zero_restrictions(szr_model(), szr_restrictions())

  identification <- object[["model"]][["sign_zero_restrictions"]]
  expect_equal(identification[["candidates"]],
               nrow(szr_model()[["posterior"]][["a"]][["coeffs"]]))
  expect_lte(identification[["accepted"]], identification[["candidates"]])
  expect_lte(identification[["effective_sample_size"]], identification[["accepted"]])

  # Every block of draws is resampled together, so the coefficients, the
  # covariance and the rotation of a row still belong to the same draw.
  kept <- identification[["effective_sample_size"]]
  for (block in c("a", "u_sigma_inv", "q")) {
    expect_equal(nrow(object[["posterior"]][[block]][["coeffs"]]), kept)
  }
  expect_identical(object[["model"]][["chains"]], 1L)

  # The rotations are read by everything that takes type = "sign", without
  # either of the two identifications knowing about the other.
  ir <- irf(object, impulse = "y", response = "r", n_ahead = 3, type = "sign")
  expect_equal(nrow(ir), 4)
  expect_gt(as.numeric(ir[1, 2]), 0)
  expect_silent(fevd(object, response = "r", n_ahead = 3, type = "sign"))
})

test_that("the number of resampled draws can be chosen", {
  set.seed(1234)
  object <- add_sign_zero_restrictions(szr_model(), szr_restrictions(), draws = 25)

  expect_equal(nrow(object[["posterior"]][["q"]][["coeffs"]]), 25)
  expect_equal(nrow(object[["posterior"]][["a"]][["coeffs"]]), 25)
})

test_that("a long-run zero restriction is imposed on the long run", {
  set.seed(1234)
  restrictions <- data.frame(impulse = "y", response = "Dp", sign = 0, horizon = Inf)
  object <- add_sign_zero_restrictions(szr_model(), restrictions)

  k <- object[["model"]][["k"]]
  draws <- bvartools:::.collect_draws(object)
  for (i in seq_len(min(5, nrow(object[["posterior"]][["q"]][["coeffs"]])))) {
    q <- matrix(object[["posterior"]][["q"]][["coeffs"]][i, ], k)
    a0 <- solve(chol(draws[[i]][["Sigma"]]), q)
    aplus <- t(cbind(draws[[i]][["A"]], 0)) %*% a0
    long_run <- t(solve(a0 - aplus[1:k, ]))
    expect_equal(long_run[2, 1], 0)
  }
})

test_that("a table with no zero restriction is sent to the cheaper function", {
  restrictions <- data.frame(impulse = "y", response = c("Dp", "r"),
                             sign = c(1, 1), horizon = 0)
  expect_error(add_sign_zero_restrictions(szr_model(), restrictions),
               "add_sign_restrictions")
})

test_that("an ordering that leaves a shock no column to draw is refused", {
  # The third variable's shock has two columns taken before it, so one zero
  # restriction is one too many.
  restrictions <- data.frame(impulse = "r", response = c("y", "Dp"),
                             sign = c(0, 1), horizon = 0)
  expect_error(add_sign_zero_restrictions(szr_model(), restrictions),
               "leaves nothing of its column to draw")
})

test_that("the table is checked the way its sibling's is", {
  object <- szr_model()

  expect_error(add_sign_zero_restrictions(object, list(impulse = "y")),
               "must be a data frame")
  expect_error(add_sign_zero_restrictions(object, szr_restrictions()[0, ]),
               "does not contain any restriction")
  expect_error(add_sign_zero_restrictions(object, data.frame(impulse = "y", sign = 0)),
               "must contain the column 'response'")
  expect_error(add_sign_zero_restrictions(object,
                                          data.frame(impulse = "gdp", response = "r",
                                                     sign = 0)),
               "which the model does not contain")
  expect_error(add_sign_zero_restrictions(object,
                                          data.frame(impulse = "y", response = "r",
                                                     sign = 2)),
               "must be -1, 0 or 1")
  expect_error(add_sign_zero_restrictions(object,
                                          data.frame(impulse = "y", response = "r",
                                                     sign = 0, horizon = -1)),
               "non-negative integers")

  contradictory <- data.frame(impulse = "y", response = c("Dp", "Dp"),
                              sign = c(0, 1), horizon = 0)
  expect_error(add_sign_zero_restrictions(object, contradictory),
               "more than one thing of the same response")

  expect_error(add_sign_zero_restrictions(object, szr_restrictions(), draws = 0),
               "single positive integer")
  expect_error(add_sign_zero_restrictions(object, szr_restrictions(), one_sided = "yes"),
               "TRUE or FALSE")
})

test_that("a model whose covariance cannot be rotated is refused", {
  expect_error(add_sign_zero_restrictions(fx_svar_fitted(), szr_restrictions()),
               "not defined for a structural model")

  no_covar <- szr_model()
  no_covar[["model"]][["error"]] <- "gamma"
  expect_error(add_sign_zero_restrictions(no_covar, szr_restrictions()),
               "error covariances are estimated")
})

test_that("the summary reports the effective sample size", {
  set.seed(1234)
  object <- add_sign_zero_restrictions(szr_model(), szr_restrictions())

  printed <- utils::capture.output(print(summary(object)))
  expect_true(any(grepl("Sign and zero restrictions", printed)))
  expect_true(any(grepl("Draws satisfying the signs", printed)))
  expect_true(any(grepl("Effective sample size", printed)))
})
