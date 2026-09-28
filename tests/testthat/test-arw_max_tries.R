# More than one rotation per draw in the sign and zero restricted sampler.
#
# The claim to test is that the importance sample stays exact: each draw's
# weight is multiplied by an unbiased estimate of the probability that a try
# succeeds, so the expected weight of a draw does not depend on max_tries. It
# is checked directly on one draw, against the single try of the paper, whose
# weight is the volume element times the indicator of success -- itself an
# unbiased estimate of the same thing. With more than one try a column whose
# signs are all reversed is flipped, which multiplies that probability by 2^s
# for the s shocks carrying signs, the same factor for every draw. Here s = 2,
# so the expected weight must be exactly four times the single try's, and a
# flip that changed the distribution of the rotation kept would break that.

mt_draw <- function(k = 3) {
  set.seed(42)
  a1 <- diag(0.5, k)
  a1[lower.tri(a1)] <- 0.15
  list(A = cbind(a1, rep(0.2, k)),
       Sigma = crossprod(matrix(stats::rnorm(k * k), k)) + diag(k))
}

# One zero restriction on the first shock and signs on the other two, which a
# single rotation satisfies about one time in twenty.
mt_setup <- function(k = 3) {
  z <- lapply(seq_len(k), function(j) matrix(0, 0, k))
  z[[1]] <- matrix(c(0, 1, 0), 1)
  s <- lapply(seq_len(k), function(j) matrix(0, 0, k))
  s[[1]] <- rbind(c(1, 0, 0), c(0, 0, 1))
  s[[2]] <- rbind(c(1, 0, 0), c(0, -1, 0))
  set.seed(7)
  w <- lapply(seq_len(k), function(j) {
    d <- k - ((j - 1) + nrow(z[[j]]))
    matrix(stats::rnorm(d * k), d, k)
  })
  list(k = k, m = k + 1, lags = 1, horizons = 0, z = z, s = s, w = w)
}

mt_weights <- function(n, max_tries, seed) {
  A <- mt_draw()
  setup <- mt_setup()
  set.seed(seed)
  out <- t(vapply(seq_len(n), function(i) {
    r <- bvartools:::.arw_draw_q(A, setup, max_tries = max_tries)
    c(weight = exp(r[["log_weight"]]), tries = r[["tries"]])
  }, numeric(2)))
  out[!is.finite(out[, "weight"]), "weight"] <- 0
  out
}

test_that("a single try is the paper's draw, and reports one try", {
  A <- mt_draw()
  setup <- mt_setup()
  set.seed(11)
  default <- bvartools:::.arw_draw_q(A, setup)
  set.seed(11)
  one <- bvartools:::.arw_draw_q(A, setup, max_tries = 1)
  expect_identical(default, one)
  expect_identical(one[["tries"]], 1L)
})

test_that("more tries leave the expected weight of a draw unchanged", {
  # A single try estimates the weight with a Bernoulli factor, so it needs
  # many more draws than the others for the same precision.
  single <- mt_weights(20000, 1, seed = 1)
  p <- mean(single[, "weight"] > 0)
  # The setup is only a test of the estimate if a single try often fails.
  expect_gt(p, 0.02)
  expect_lt(p, 0.5)

  for (max_tries in c(3, 50)) {
    several <- mt_weights(2000, max_tries, seed = 2)
    difference <- mean(several[, "weight"]) / 4 - mean(single[, "weight"])
    se <- sqrt(stats::var(several[, "weight"]) / 16 / nrow(several) +
                 stats::var(single[, "weight"]) / nrow(single))
    expect_lt(abs(difference), 4 * se)
    # And they find a rotation far more often than one try does.
    expect_gt(mean(several[, "weight"] > 0), p)
    expect_lte(max(several[, "tries"]), max_tries)
  }
})

test_that("the rotation kept satisfies the restrictions", {
  A <- mt_draw()
  setup <- mt_setup()
  set.seed(5)
  r <- bvartools:::.arw_draw_q(A, setup, max_tries = 1000)
  expect_equal(dim(r[["q"]]), c(3, 3))
  impact <- t(solve(solve(chol(A[["Sigma"]]), r[["q"]])))
  expect_equal(impact[2, 1], 0, tolerance = 1e-10)
  expect_true(all(c(impact[1, 1], impact[3, 1], impact[1, 2], -impact[2, 2]) > 0))
})

test_that("the worker refuses fewer than one try", {
  expect_error(bvartools:::.arw_draw_q(mt_draw(), mt_setup(), max_tries = 0), "at least 1")
})

test_that("a model the paper's single try cannot identify is identified with more", {
  object <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                             iterations = 60, burnin = 10)
  object <- add_priors(object, coef = list(v_i = 1, v_i_det = 0.1),
                       sigma = list(df = "k", scale = 1))
  object <- add_posterior_coefficients(add_seed(add_initial_values(object), 216))
  # Output growth does not move inflation on impact, and every shock carries
  # a sign on every variable it is not restricted on.
  restrictions <- data.frame(
    impulse = rep(c("y", "Dp", "r"), c(3, 3, 3)),
    response = rep(c("y", "Dp", "r"), 3),
    sign = c(1, 0, 1, -1, 1, 1, 1, 1, -1))

  set.seed(1)
  one <- tryCatch(suppressWarnings(add_sign_zero_restrictions(object, restrictions)),
                  error = function(e) NULL)
  set.seed(1)
  many <- suppressWarnings(add_sign_zero_restrictions(object, restrictions, max_tries = 500))
  record <- many[["model"]][["sign_zero_restrictions"]]

  expect_identical(record[["max_tries"]], 500L)
  expect_gt(record[["tries"]], record[["candidates"]])
  expect_gt(record[["accepted"]],
            if (is.null(one)) 0 else one[["model"]][["sign_zero_restrictions"]][["accepted"]])
  expect_output(print(summary(many)), "Rotations drawn")

  expect_error(add_sign_zero_restrictions(object, restrictions, max_tries = 0),
               "single positive integer")
  expect_error(add_sign_zero_restrictions(object, restrictions, max_tries = 2.5),
               "single positive integer")
})
