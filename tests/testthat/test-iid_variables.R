# Endogenous variables whose equations carry no coefficients.
#
# The numerics are the vendored BayesTS core, which tests them upstream against
# a hand-reduced system. What is checked here is the R side: that `iid` names
# variables rather than counts them, that a data set in the wrong order is told
# what to do about it, and that a model estimated under the restriction is an
# ordinary 'bvarmodel' everywhere else.

# Three series where the first is white noise correlated with the errors of the
# other two -- the shape a high-frequency surprise gives a model.
iid_data <- function() {
  set.seed(4321)
  tt <- 120
  surprise <- stats::rnorm(tt)
  y <- stats::rnorm(tt) + 0.6 * surprise
  p <- stats::rnorm(tt) - 0.4 * surprise
  stats::ts(cbind(surprise = surprise, y = y, p = p), start = 1)
}

iid_model <- function(...) {
  object <- create_bvarmodel(iid_data(), p = 1, deterministic = "const",
                             iterations = fx_iterations, burnin = fx_burnin, ...)
  object <- add_priors(object, coef = list(v_i = 1, v_i_det = 0.1),
                       sigma = list(df = "k", scale = 1))
  add_initial_values(object)
}

# The positions of `a` that belong to the restricted equations. Coefficients are
# stored as vec of the k by n_x matrix, so equation i owns every k-th one.
restricted_positions <- function(object) {
  k <- object[["model"]][["k"]]
  n_iid <- object[["model"]][["n_iid"]]
  unlist(lapply(seq_len(n_iid), function(i) seq(i, ncol(object[["posterior"]][["a"]][["coeffs"]]), by = k)))
}

test_that("the named equations draw no coefficients at all", {
  object <- add_posterior_coefficients(add_seed(iid_model(iid = "surprise"), 216))

  a <- object[["posterior"]][["a"]][["coeffs"]]
  positions <- restricted_positions(object)

  expect_equal(max(abs(a[, positions])), 0)
  expect_gt(max(abs(a[, -positions])), 0)

  # The restricted variable is still a variable: it is a regressor in the other
  # equations and the error covariance covers it.
  expect_equal(ncol(object[["posterior"]][["u_sigma_inv"]][["coeffs"]]),
               object[["model"]][["k"]]^2)
  sigma <- solve(matrix(colMeans(object[["posterior"]][["u_sigma_inv"]][["coeffs"]]), 3))
  expect_gt(abs(sigma[1, 2]), 0)

  # And everything downstream treats it as the VAR it is.
  ir <- irf(object, impulse = "surprise", response = "y", n_ahead = 3, type = "oir")
  expect_equal(nrow(ir), 4)
  expect_silent(fevd(object, response = "y", n_ahead = 3))
})

test_that("two restricted variables are two restricted equations", {
  object <- add_posterior_coefficients(
    add_seed(iid_model(iid = c("surprise", "y")), 216))

  expect_identical(object[["model"]][["n_iid"]], 2L)
  expect_equal(max(abs(object[["posterior"]][["a"]][["coeffs"]][, restricted_positions(object)])), 0)

  # Naming them in the other order is the same restriction: which of the
  # leading variables comes first makes no difference.
  swapped <- iid_model(iid = c("y", "surprise"))
  expect_identical(swapped[["model"]][["n_iid"]], 2L)
})

test_that("a model without the restriction is untouched", {
  object <- iid_model()
  expect_null(object[["model"]][["n_iid"]])

  specifications <- get_model_specifications(object)
  expect_false("n_iid" %in% names(specifications))

  # And the restriction is what the specification says it is.
  expect_equal(get_model_specifications(iid_model(iid = "surprise"))[["n_iid"]], 1)
})

test_that("the restriction survives a file", {
  object <- add_posterior_coefficients(add_seed(iid_model(iid = "surprise"), 216))

  file <- tempfile(fileext = ".h5")
  on.exit(unlink(file))
  write_to_hdf5(object, file)
  back <- read_model_from_hdf5(file)

  expect_identical(as.integer(back[["model"]][["n_iid"]]), 1L)
  expect_equal(max(abs(back[["posterior"]][["a"]][["coeffs"]][, restricted_positions(back)])), 0)
})

test_that("variables that are not the leading columns are refused", {
  data <- iid_data()

  expect_error(create_bvarmodel(data, p = 1, iid = "y", iterations = 10, burnin = 5),
               "not the first 1 columns")
  # The message names what is there instead, so that the fix is to hand.
  expect_error(create_bvarmodel(data, p = 1, iid = "p", iterations = 10, burnin = 5),
               "'surprise'")

  # Reordering the data is the fix, and it works.
  reordered <- data[, c("p", "surprise", "y")]
  reordered <- stats::ts(reordered, start = stats::start(data),
                         frequency = stats::frequency(data))
  object <- create_bvarmodel(reordered, p = 1, iid = "p", iterations = 10, burnin = 5)
  expect_identical(object[["model"]][["n_iid"]], 1L)
})

test_that("the argument is checked before anything is built", {
  data <- iid_data()
  build <- function(...) create_bvarmodel(data, p = 1, iterations = 10, burnin = 5, ...)

  expect_error(build(iid = 1), "must name endogenous variables")
  expect_error(build(iid = character(0)), "must name endogenous variables")
  expect_error(build(iid = "gdp"), "which 'data' does not contain")
  expect_error(build(iid = c("surprise", "surprise")), "more than once")
  expect_error(build(iid = c("surprise", "y", "p")), "leaves no equation with dynamics")
  # A structural model with a Wishart error term is refused before this for a
  # reason of its own, so the combination is asked for where it can be reached.
  expect_error(build(iid = "surprise", structural = TRUE, error = "gamma"),
               "'iid' and 'structural' cannot be combined")
  expect_error(build(iid = "surprise", varsel = "bvs"), "'iid' and 'varsel' cannot be combined")
  expect_error(build(iid = "surprise", tvp = TRUE), "'iid' and 'tvp' cannot be combined")
})

test_that("the samplers that cannot honour the restriction never see it", {
  # Refusing a time-varying model where it is written rather than where it is
  # estimated is the whole point of the check above; this is the other side of
  # it, that the four algorithms which do read it are the ones a model with
  # `iid` can be given.
  for (error in c("wishart", "gamma", "sv", "ald")) {
    object <- create_bvarmodel(iid_data(), p = 1, deterministic = "const",
                               iid = "surprise", error = error,
                               iterations = 10, burnin = 5)
    expect_identical(object[["model"]][["n_iid"]], 1L)
    expect_true(object[["model"]][["algorithm"]] %in%
                  c("VarNormalWishart", "VarNormalGamma", "VarNormalStochvol", "VarNormalAld"))
  }
})
