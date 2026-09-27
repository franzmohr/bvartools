# The persistence profiles of the cointegrating relations. What the statistic
# means is checked here -- that it starts at one, that it falls for a relation
# that is genuinely cointegrating, and that it is refused where there is no
# relation to follow. The algebra of the level VAR it is built on belongs to
# test-vec_to_var.R.

test_that("a persistence profile starts at one", {
  # PP(0) is beta' Sigma beta over itself, so it is one for every draw and
  # every relation whatever the model says.
  model <- fx_vec_fitted()
  pp <- persistence_profiles(model, n_ahead = 8, keep_draws = TRUE)

  expect_type(pp, "list")
  expect_length(pp, model[["model"]][["rank"]])
  for (relation in pp) {
    expect_equal(unname(relation[, 1]), rep(1, nrow(relation)))
  }
})

test_that("the summary carries the median and the band", {
  model <- fx_vec_fitted()
  pp <- persistence_profiles(model, n_ahead = 6, ci = 0.9)

  expect_s3_class(pp, "bvecpp")
  for (relation in pp) {
    expect_identical(colnames(relation), c("median", "lower", "upper"))
    expect_identical(nrow(relation), 7L)
    expect_equal(unname(relation[1, "median"]), 1)
    expect_true(all(relation[, "lower"] <= relation[, "median"]))
    expect_true(all(relation[, "median"] <= relation[, "upper"]))
    expect_true(all(relation >= 0))
  }
})

test_that("a stationary relation reverts and a random walk does not", {
  # Two systems built by hand, run through the same recursion the method uses,
  # so the estimand is checked rather than the plumbing: a relation that is
  # stationary has a profile that falls to zero, and one that is a random walk
  # has a profile that stays at one.
  profile <- function(A, beta, n_ahead) {
    k <- nrow(A)
    sigma <- diag(k)
    psi <- diag(k)
    out <- numeric(n_ahead + 1)
    scale <- as.numeric(crossprod(beta, sigma %*% beta))
    for (h in 0:n_ahead) {
      m <- crossprod(beta, psi)
      out[h + 1] <- as.numeric(m %*% sigma %*% t(m)) / scale
      psi <- A %*% psi
    }
    out
  }

  # y1 and y2 share a stochastic trend and y1 - y2 is stationary.
  A <- matrix(c(0.5, 0.5, 0.5, 0.5), 2)
  reverting <- profile(A, matrix(c(1, -1)), 12)
  expect_equal(reverting[1], 1)
  expect_lt(reverting[13], 0.01)
  expect_true(all(diff(reverting) <= 1e-12))

  # The sum of two random walks is not stationary and does not revert.
  walk <- profile(diag(2), matrix(c(1, 1)), 12)
  expect_equal(walk, rep(1, 13))
})

test_that("the method refuses where there is no relation to follow", {
  expect_error(persistence_profiles(fx_var_fitted()),
               "no cointegrating relations")

  model <- fx_vec_fitted()
  model[["model"]][["rank"]] <- 0L
  expect_error(persistence_profiles(model), "has no persistence profiles")

  model <- fx_vec_fitted()
  model[["posterior"]][["beta"]] <- NULL
  expect_error(persistence_profiles(model), "no posterior draws of 'beta'")

  expect_error(persistence_profiles(fx_vec_fitted(), ci = 0), "'ci'")
  expect_error(persistence_profiles(fx_vec_fitted(), n_ahead = 0), "'n_ahead'")
})
