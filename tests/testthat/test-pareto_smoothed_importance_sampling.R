# Pareto smoothed importance sampling, the reweighting behind
# add_sign_zero_restrictions().
#
# What this file checks is the estimator itself, against distributions whose
# tail shape is known by construction: that the fitted shape recovers it, that
# the smoothed weights remain a probability distribution, and that inputs too
# short or too degenerate to say anything about a tail are declined rather than
# guessed at. What add_sign_zero_restrictions() does with the result -- the
# `smooth` argument, the warnings, what summary() prints -- is in
# test-add_sign_zero_restrictions.R.

# Generalised Pareto draws with shape k: the tail this estimator exists to
# describe. For k > 0 the qth moment is infinite once q >= 1 / k, so k = 0.5 is
# the edge of a finite variance and k = 0.8 is past it.
rgpd <- function(n, k, sigma = 1) sigma * (stats::runif(n)^(-k) - 1) / k

test_that("the fitted shape recovers the tail it was drawn from", {
  set.seed(216)

  for (k in c(0.2, 0.5, 0.8)) {
    estimates <- replicate(10, .psis_smooth(log(rgpd(4000, k)))$pareto_k)
    # Calibrated rather than asserted: over twenty replications the estimator
    # had a standard deviation of about 0.1 at each of these shapes, so a
    # tolerance of 0.12 on the mean of ten is loose enough not to flake and
    # tight enough to fail if the fit is wrong.
    expect_equal(mean(estimates), k, tolerance = 0.12)
  }

  # Exponential weights have no Pareto tail at all, which is k = 0.
  expect_lt(abs(mean(replicate(10, .psis_smooth(log(stats::rexp(4000)))$pareto_k))), 0.12)
})

test_that("the smoothed weights are still a probability distribution", {
  set.seed(216)
  out <- .psis_smooth(log(rgpd(1000, 0.7)))

  weights <- exp(out$log_weights)
  expect_equal(sum(weights), 1)
  expect_true(all(weights >= 0))
  expect_length(out$log_weights, 1000)
})

test_that("smoothing rescues a sample one draw was carrying", {
  set.seed(216)

  share <- function(w) max(w) / sum(w)
  ess <- function(w) { w <- w / sum(w); 1 / sum(w^2) }

  # The pathology this exists for: one draw far out in the tail. Raw, it takes
  # essentially the whole sample; smoothed, it takes a fifth of it.
  lw <- c(stats::rnorm(500), 15)
  raw <- exp(lw - max(lw))
  smoothed <- exp(.psis_smooth(lw)$log_weights)

  expect_gt(share(raw), 0.9)
  expect_lt(share(smoothed), 0.3)
  expect_lt(ess(raw), 2)
  expect_gt(ess(smoothed), 10)

  # A heavy tail without a single dominating draw gains less, and gains it in
  # the effective sample size rather than in the largest share -- normalising
  # can leave the biggest weight a larger share of a smaller total, so the
  # share is not the invariant here.
  heavy <- log(rgpd(1000, 0.9))
  expect_gt(ess(exp(.psis_smooth(heavy)$log_weights)),
            ess(exp(heavy - max(heavy))))

  # And where the weights are already even, it is close to a no-op.
  mild <- stats::rnorm(1000) * 0.1
  expect_equal(ess(exp(.psis_smooth(mild)$log_weights)),
               ess(exp(mild - max(mild))), tolerance = 0.01)
})

test_that("a draw that carries no weight stays out of the sample", {
  set.seed(216)
  lw <- c(stats::rnorm(200), rep(-Inf, 100))
  out <- .psis_smooth(lw)

  # A rejected draw is not a draw with a tiny weight; it is not in the sample,
  # and letting it into the tail fit would describe a distribution that has no
  # draws in it.
  expect_true(all(is.infinite(out$log_weights[201:300])))
  expect_true(all(out$log_weights[201:300] < 0))
  expect_equal(sum(exp(out$log_weights[1:200])), 1)

  # The finite draws are smoothed as though the others had never been passed.
  expect_equal(out$log_weights[1:200], .psis_smooth(lw[1:200])$log_weights)
})

test_that("an input with no tail to speak of is declined rather than guessed at", {
  # Too few draws for a tail of at least five.
  expect_true(is.na(.psis_smooth(stats::rnorm(10))$pareto_k))
  # Every weight identical: there is no shape to fit, and fitting one would
  # divide by zero rather than fail loudly.
  expect_true(is.na(.psis_smooth(rep(1, 500))$pareto_k))
  # Nothing at all.
  expect_true(is.na(.psis_smooth(rep(-Inf, 50))$pareto_k))

  # In each case the weights still come back usable, just unsmoothed.
  out <- .psis_smooth(rep(1, 500))
  expect_equal(sum(exp(out$log_weights)), 1)
  expect_equal(length(unique(round(out$log_weights, 12))), 1)
})

test_that("the scale of the log weights does not change the answer", {
  set.seed(216)
  lw <- log(rgpd(1000, 0.6))

  # Importance weights are unnormalised, so adding a constant to every log
  # weight must leave both the smoothed weights and the fitted shape alone --
  # including a constant large enough that exp() of it would overflow.
  shifted <- .psis_smooth(lw + 700)
  expect_equal(shifted$log_weights, .psis_smooth(lw)$log_weights)
  expect_equal(shifted$pareto_k, .psis_smooth(lw)$pareto_k)
  expect_false(anyNA(shifted$log_weights))
})

test_that("the implementation agrees with the reference one in loo", {
  # The numerics here are a translation of loo::psis(), which the authors of
  # the method maintain, rather than an independent derivation -- so the thing
  # worth checking is that the translation is faithful, on weight distributions
  # whose tails run from lighter than exponential to heavier than Cauchy.
  #
  # loo is in Suggests rather than Imports: the package does not use it, and
  # depending on it to check sixty lines of arithmetic would be the tail
  # wagging the dog. Where it is not installed this check is skipped and the
  # properties tested above are what stands behind the implementation.
  skip_if_not_installed("loo")

  set.seed(216)
  cases <- list(
    "normal log ratios" = stats::rnorm(1000),
    "wider" = stats::rnorm(300) * 3,
    "heavy tailed" = log(abs(stats::rt(500, df = 1))),
    "exponential" = stats::rexp(2000),
    "generalised Pareto" = log(rgpd(1500, 0.7))
  )

  for (name in names(cases)) {
    ours <- .psis_smooth(cases[[name]])
    # loo warns where the shape is above its own threshold, which three of
    # these cases are on purpose: a heavy tail is what this is checking.
    theirs <- suppressWarnings(loo::psis(cases[[name]], r_eff = 1))
    reference <- as.vector(stats::weights(theirs, log = TRUE, normalize = TRUE))

    expect_equal(ours$pareto_k, theirs$diagnostics$pareto_k,
                 tolerance = 1e-8, info = name)
    expect_equal(ours$log_weights, reference, tolerance = 1e-10, info = name)
  }
})
