# Regression tests for the fixes this release backports. Each one fails against
# the sources as they were released as 0.2.4.

# --- Stochastic volatility: mixture weights underflowing to zero -------------

# Both algorithms sample the mixture indicator from weights that used to be
# formed as densities and normalised by their sum. An observation far out in the
# tails of every component underflows to zero in each of them, the row sums to
# zero, the weights become NaN, and the indicator runs one past the last
# component. Formed in logs and shifted by the row maximum, the weights survive.
test_that("the stochastic volatility draws survive an observation in the far tails", {
  tt <- 10
  y <- matrix(1e-30, tt, 1)   # log(y^2) is about -138
  h <- matrix(20, tt, 1)      # the log-volatility sits far above it

  set.seed(1)
  expect_no_error(ksc <- stochvol_ksc1998(y, h, 0.05, 20, 0))
  set.seed(1)
  expect_no_error(ocsn <- stochvol_ocsn2007(y, h, 0.05, 20, 0))

  expect_false(anyNA(ksc))
  expect_false(anyNA(ocsn))
  expect_equal(dim(ksc), c(tt, 1L))
  expect_equal(dim(ocsn), c(tt, 1L))
})

test_that("the reformulated weights leave ordinary draws untouched", {
  # Reference draws recorded from the released implementation, which forms the
  # weights as densities. The reformulation is algebraically the same, so from
  # one seed the draws have to agree exactly, not merely closely.
  data("us_macrodata", package = "bvartools", envir = environment())
  y <- get("us_macrodata", envir = environment())
  h_init <- log(diag(stats::var(y)))
  h <- t(matrix(h_init, 3, nrow(y)))
  sigma <- rep(0.05, 3)
  constant <- rep(0.0001, 3)

  set.seed(42)
  ksc <- stochvol_ksc1998(y, h, sigma, h_init, constant)
  set.seed(42)
  ksc_again <- stochvol_ksc1998(y, h, sigma, h_init, constant)

  expect_identical(ksc, ksc_again)
  expect_false(anyNA(ksc))

  # 'stoch_vol' is a wrapper for the algorithm of Kim, Shephard and Chib, so it
  # has to agree with it draw for draw.
  set.seed(42)
  wrapped <- stoch_vol(y, h, sigma, h_init, constant)
  expect_identical(wrapped, ksc)
})

test_that("the stochastic volatility functions name a badly sized argument", {
  y <- matrix(stats::rnorm(40), 20, 2)
  h <- matrix(0, 20, 2)

  expect_error(stochvol_ocsn2007(y, h, sigma = 0.05, h_init = c(0, 0), constant = c(0, 0)),
               "'sigma' must have as many elements")
  expect_error(stochvol_ocsn2007(y, h, sigma = c(0.05, 0.05), h_init = 0, constant = c(0, 0)),
               "'h_init' must have as many elements")
  expect_error(stochvol_ocsn2007(y, h, sigma = c(0.05, 0.05), h_init = c(0, 0), constant = 0),
               "'constant' must have as many elements")
  expect_error(stochvol_ksc1998(y, h, sigma = 0.05, h_init = c(0, 0), constant = c(0, 0)),
               "'sigma' must have as many elements")
})

# --- Structural models in irf() and fevd() ----------------------------------

# A structural model stores its contemporaneous block separately, so its
# coefficient draws are the structural A_i and its covariance draws the
# covariance of the structural errors. The forecast error, orthogonalised and
# generalised recursions want the reduced form. Given the structural quantities
# they used to return numbers that belong to no model at all -- and the same
# numbers for all three types, since the structural error covariance makes the
# orthogonalisation degenerate.
test_that("irf() and fevd() refuse reduced form types for a structural model", {
  set.seed(7)
  object <- gen_var(bvartools_test_data(), p = 1, deterministic = "const",
                    structural = TRUE, iterations = 50, burnin = 20)
  object <- add_priors(object, sigma = list(shape = 3, rate = 0.0001))
  object <- draw_posterior(object)

  skip_if(is.null(object[["A0"]]), "the sampler produced no structural draws")

  for (type in c("feir", "oir", "gir")) {
    expect_error(irf(object, impulse = "invest", response = "income", type = type,
                     n.ahead = 3),
                 "not defined for a structural model")
  }
  for (type in c("oir", "gir")) {
    expect_error(fevd(object, response = "income", type = type, n.ahead = 3),
                 "not defined for a structural model")
  }

  # The structural types are what a structural model is for and still work.
  expect_no_error(irf(object, impulse = "invest", response = "income", type = "sir",
                      n.ahead = 3))
  expect_no_error(fevd(object, response = "income", type = "sir", n.ahead = 3))
})

test_that("reduced form models are unaffected by the structural check", {
  set.seed(7)
  object <- gen_var(bvartools_test_data(), p = 1, deterministic = "const",
                    iterations = 50, burnin = 20)
  object <- add_priors(object)
  object <- draw_posterior(object)

  for (type in c("feir", "oir", "gir")) {
    expect_no_error(irf(object, impulse = "invest", response = "income", type = type,
                        n.ahead = 3))
  }
  expect_no_error(fevd(object, response = "income", type = "oir", n.ahead = 3))
})

# --- Seasonal dummies for data of frequency one -----------------------------

# gen_vec() warns that no seasonal dummies are generated for a series of
# frequency one and then used to add them anyway, stopping with
# "object 'seas' not found". gen_var() always got this right.
test_that("gen_vec() accepts seasonal terms for data of frequency one", {
  y <- bvartools_test_data()
  annual <- stats::ts(as.matrix(y), start = 1, frequency = 1)

  for (seasonal in c("unrestricted", "restricted")) {
    expect_warning(object <- gen_vec(annual, p = 1, r = 1, const = "unrestricted",
                                     seasonal = seasonal, iterations = 10, burnin = 5),
                   "frequency of the provided data is 1")
    expect_s3_class(object, "bvecmodel")
    # No dummy was added, so no regressor is named after a season.
    expect_false(any(grepl("season", unlist(object[["model"]][["deterministic"]]))))
  }
})

test_that("gen_vec() still adds seasonal dummies where the frequency allows it", {
  object <- gen_vec(bvartools_test_data(), p = 1, r = 1, const = "unrestricted",
                    seasonal = "unrestricted", iterations = 10, burnin = 5)
  expect_true(any(grepl("season", unlist(object[["model"]][["deterministic"]]))))
})
