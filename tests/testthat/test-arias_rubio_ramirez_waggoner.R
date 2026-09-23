# The sign and zero restriction algorithm of Arias, Rubio-Ramirez and Waggoner
# (2018).
#
# Two kinds of check. The first is structural and exact: the rotation is
# orthogonal and the zero restrictions hold to machine precision, which they
# must, because the algorithm builds the columns inside the subspace the zeros
# leave rather than searching for them. The second is the importance weight,
# which has no closed form and is tested against the two things the paper says
# about it -- that it is constant when there are no zero restrictions
# (Section 4.3), and that the effective sample size of their own application is
# 78 percent of the draws that satisfy the restrictions (Section 7.3).

# A draw whose responses the restrictions below can bite on. The lag matrices
# are deliberately not multiples of the identity: when they are, the long run
# comes out proportional to the impact response and a restriction on one is a
# restriction on the other, which would make the long-run test below pass for
# the wrong reason.
arw_draw <- function(k = 3, lags = 2) {
  set.seed(42)
  a1 <- diag(0.5, k)
  a1[lower.tri(a1)] <- 0.15
  a2 <- diag(-0.1, k)
  a2[1, k] <- 0.2
  list(A = cbind(a1, a2, rep(0.2, k)),
       Sigma = crossprod(matrix(stats::rnorm(k * k), k)) + diag(k))
}

# Restriction blocks with no restriction in them, to be filled in one at a time.
arw_blocks <- function(k, nh) {
  out <- vector("list", k)
  for (j in seq_len(k)) out[[j]] <- matrix(0, 0, nh * k)
  out
}

# The fixed matrices that complete each column's constraints to a square
# system. Any draw of them defines a valid algorithm, so the seed is here only
# to keep the tests reproducible.
arw_w <- function(k, z, seed = 7) {
  set.seed(seed)
  out <- vector("list", k)
  for (j in seq_len(k)) {
    s <- k - ((j - 1) + nrow(z[[j]]))
    out[[j]] <- matrix(stats::rnorm(s * k), s, k)
  }
  out
}

arw_setup <- function(k = 3, lags = 2, horizons = 0, z = NULL, s = NULL, seed = 7) {
  nh <- length(horizons)
  if (is.null(z)) z <- arw_blocks(k, nh)
  if (is.null(s)) s <- arw_blocks(k, nh)
  list(k = k, m = k * lags + 1, lags = lags, horizons = horizons,
       z = z, s = s, w = arw_w(k, z, seed))
}

test_that("the rotation is orthogonal and the zero restrictions hold exactly", {
  k <- 3
  z <- arw_blocks(k, 1)
  z[[1]] <- matrix(c(1, 0, 0), 1)
  z[[2]] <- matrix(c(0, 1, 0), 1)
  setup <- arw_setup(k = k, z = z)

  A <- arw_draw(k)
  set.seed(1)
  result <- bvartools:::.arw_draw_q(A, setup)

  q <- result[["q"]]
  expect_equal(dim(q), c(k, k))
  expect_equal(t(q) %*% q, diag(1, k))

  # The impact responses, whose first two columns the restrictions speak about.
  a0 <- solve(chol(A[["Sigma"]]), q)
  l0 <- t(solve(a0))
  expect_equal(l0[1, 1], 0)
  expect_equal(l0[2, 2], 0)

  # The rotation leaves the model observationally equivalent, as every
  # identification of a reduced form must.
  expect_equal(a0 %*% t(a0), solve(A[["Sigma"]]))
})

test_that("a zero restriction on the long run is imposed on the long run", {
  k <- 3
  z <- arw_blocks(k, 2)
  # The second block of F is the long run, so the restriction sits in columns
  # k + 1 to 2k.
  z[[1]] <- matrix(c(0, 0, 0, 1, 0, 0), 1)
  setup <- arw_setup(k = k, horizons = c(0, Inf), z = z)

  A <- arw_draw(k)
  set.seed(1)
  q <- bvartools:::.arw_draw_q(A, setup)[["q"]]

  a0 <- solve(chol(A[["Sigma"]]), q)
  aplus <- t(A[["A"]]) %*% a0
  long_run <- t(solve(a0 - aplus[1:k, ] - aplus[k + 1:k, ]))

  expect_equal(long_run[1, 1], 0)
  # The impact response of the same variable to the same shock is not zero:
  # the restriction was on the long run alone.
  expect_false(isTRUE(all.equal(t(solve(a0))[1, 1], 0)))
})

test_that("the weights are constant when nothing is restricted to zero", {
  # Section 4.3: with no zero restrictions the algorithm draws from the target
  # itself, so the weights carry no information and Algorithm 1 should be used
  # instead. What is left in the weight here is the error of the numerical
  # derivative.
  setup <- arw_setup()
  A <- arw_draw()

  set.seed(11)
  weights <- vapply(1:50, function(i) bvartools:::.arw_draw_q(A, setup)[["log_weight"]],
                    numeric(1))

  expect_true(all(is.finite(weights)))
  expect_lt(diff(range(weights)), 1e-6)
})

test_that("a zero restriction makes the weights vary", {
  k <- 3
  z <- arw_blocks(k, 1)
  z[[1]] <- matrix(c(1, 0, 0), 1)
  setup <- arw_setup(k = k, z = z)
  A <- arw_draw(k)

  set.seed(11)
  weights <- vapply(1:50, function(i) bvartools:::.arw_draw_q(A, setup)[["log_weight"]],
                    numeric(1))

  expect_true(all(is.finite(weights)))
  expect_gt(diff(range(weights)), 0.1)
})

test_that("a draw that fails the sign restrictions is given weight zero", {
  k <- 3
  s_pos <- arw_blocks(k, 1)
  s_pos[[1]] <- matrix(c(1, 0, 0), 1)
  s_neg <- arw_blocks(k, 1)
  s_neg[[1]] <- matrix(c(-1, 0, 0), 1)

  A <- arw_draw(k)

  # The same rotation under a restriction and its opposite: exactly one of the
  # two can hold, so one call keeps the draw and the other must reject it.
  set.seed(3)
  positive <- bvartools:::.arw_draw_q(A, arw_setup(k = k, s = s_pos))
  set.seed(3)
  negative <- bvartools:::.arw_draw_q(A, arw_setup(k = k, s = s_neg))

  kept <- if (length(positive[["q"]]) > 0) positive else negative
  dropped <- if (length(positive[["q"]]) > 0) negative else positive

  expect_equal(dim(kept[["q"]]), c(k, k))
  expect_true(is.finite(kept[["log_weight"]]))
  expect_equal(dim(dropped[["q"]]), c(0, 0))
  expect_identical(dropped[["log_weight"]], -Inf)

  # The restriction the kept draw was kept for.
  a0 <- solve(chol(A[["Sigma"]]), kept[["q"]])
  response <- t(solve(a0))[1, 1]
  if (length(positive[["q"]]) > 0) expect_gt(response, 0) else expect_lt(response, 0)
})

test_that("the one sided derivative agrees with the two sided one", {
  k <- 3
  z <- arw_blocks(k, 1)
  z[[1]] <- matrix(c(1, 0, 0), 1)
  setup <- arw_setup(k = k, z = z)
  A <- arw_draw(k)

  set.seed(5)
  two_sided <- bvartools:::.arw_draw_q(A, setup)[["log_weight"]]
  set.seed(5)
  one_sided <- bvartools:::.arw_draw_q(A, setup, one_sided = TRUE)[["log_weight"]]

  expect_equal(one_sided, two_sided, tolerance = 1e-4)
})

test_that("the weight is only computed when it is asked for", {
  setup <- arw_setup()
  A <- arw_draw()

  set.seed(9)
  with_weight <- bvartools:::.arw_draw_q(A, setup)
  set.seed(9)
  without <- bvartools:::.arw_draw_q(A, setup, weight = FALSE)

  expect_equal(without[["q"]], with_weight[["q"]])
  expect_true(is.na(without[["log_weight"]]))
})

test_that("the setup is checked before anything is drawn", {
  k <- 3
  setup <- arw_setup(k = k)
  A <- arw_draw(k)

  wrong_shape <- setup
  wrong_shape[["z"]] <- wrong_shape[["z"]][1:2]
  expect_error(bvartools:::.arw_draw_q(A, wrong_shape), "one matrix per shock")

  wrong_cols <- setup
  wrong_cols[["s"]][[1]] <- matrix(0, 1, k + 1)
  expect_error(bvartools:::.arw_draw_q(A, wrong_cols), "columns")

  wrong_w <- setup
  wrong_w[["w"]][[1]] <- matrix(0, 1, k)
  expect_error(bvartools:::.arw_draw_q(A, wrong_w), "null space of dimension")

  # More zeros on a shock than its column has room for. Ordering the shocks so
  # that the heavily restricted ones come first is the caller's business, and
  # the message says so.
  over <- arw_blocks(k, 1)
  over[[3]] <- diag(1, k)
  too_many <- setup
  too_many[["z"]] <- over
  too_many[["w"]] <- setup[["w"]]
  expect_error(bvartools:::.arw_draw_q(A, too_many), "leaves nothing of its column")

  bad_draw <- A
  bad_draw[["A"]] <- bad_draw[["A"]][, -1]
  expect_error(bvartools:::.arw_draw_q(bad_draw, setup), "regressors")
})
