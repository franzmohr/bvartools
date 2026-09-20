# The two discounted models, which are not samplers.
#
# What is checked here is the layout -- the two /model attributes, the matrix
# normal prior under its own names, the fixed cointegration matrix, the SUR
# matrix that must not be there -- the refusals, each of which BayesTS would
# otherwise raise against a file that had already been written, and the shape of
# what the filter returns. The filter itself is the vendored BayesTS core's and
# is pinned upstream; what these tests own is the boundary, which is where the
# closed form stops looking like a chain.

discount_vec <- function(r = 1, ...) {
  data("e6", package = "bvartools", envir = environment())
  create_bvecmodel(e6 * 100, p = 2, r = r, const = "unrestricted",
                   algorithm = "discount", iterations = 10, burnin = 0,
                   thin = 1, ...)
}

discount_var <- function(...) {
  data("e1", package = "bvartools", envir = environment())
  create_bvarmodel(diff(log(e1)) * 100, p = 1, deterministic = "const",
                   algorithm = "discount", iterations = 10, burnin = 0,
                   thin = 1, ...)
}

discount_vec_with_priors <- function(r = 1, ...) {
  add_priors(discount_vec(r = r, ...),
             coef = list(v_i = 1, v_i_det = 1 / 10),
             sigma = list(df = "k", scale = 1))
}

fitted_discount_vec <- function(r = 1, ...) {
  add_initial_values(discount_vec_with_priors(r = r, ...))
}


test_that("the discount algorithm is named and its discounts are carried", {

  model <- discount_vec(delta_beta = 0.98, delta_sigma = 0.95)
  expect_equal(model[["model"]][["algorithm"]], "VecTvpDiscount")
  expect_equal(model[["model"]][["delta_beta"]], 0.98)
  expect_equal(model[["model"]][["delta_sigma"]], 0.95)
  # The coefficients of this model drift, whatever delta_beta happens to be.
  expect_true(model[["model"]][["tvp"]])

  var_model <- discount_var(delta_beta = 0.99)
  expect_equal(var_model[["model"]][["algorithm"]], "VarTvpDiscount")
  expect_equal(var_model[["model"]][["delta_sigma"]], 1)
})


test_that("a vector of discounts or ranks gives one model per combination", {

  models <- discount_vec(r = 0:2, delta_beta = c(0.98, 1),
                         delta_sigma = c(0.95, 1))
  expect_s3_class(models, "modellist")
  expect_length(models, 12)

  grid <- t(vapply(models, function(m) {
    c(m[["model"]][["rank"]], m[["model"]][["delta_beta"]],
      m[["model"]][["delta_sigma"]])
  }, numeric(3)))
  expect_equal(nrow(unique(grid)), 12)
  expect_equal(sort(unique(grid[, 1])), c(0, 1, 2))
})


test_that("the specifications the discounted models do not have are refused", {

  expect_error(discount_vec(error = "sv"), "wishart")
  expect_error(discount_vec(varsel = "bvs"), "[Vv]ariable selection")
  expect_error(discount_vec(structural = TRUE), "[Ss]tructural")
  expect_error(discount_vec(burnin = 100), "burnin")
  expect_error(discount_vec(thin = 2), "thin")
  expect_error(discount_vec(delta_beta = 0), "delta_beta")
  expect_error(discount_vec(delta_beta = 1.1), "delta_beta")

  # The degrees of freedom settle at 1 / (1 - delta_sigma), and e6 has two
  # variables, so anything at or below 0.5 leaves the inverse Wishart improper.
  expect_error(discount_vec(delta_sigma = 0.4), "too small")
  expect_s3_class(discount_vec(delta_sigma = 0.9), "bvecmodel")
})


test_that("the coefficient prior is a matrix normal under its own names", {

  model <- add_priors(discount_vec(r = 1),
                      coef = list(v_i = 4, v_i_det = 1 / 10),
                      sigma = list(df = "k", scale = 1))

  k <- model[["model"]][["k"]]
  n_design <- 1 + ncol(model[["data"]][["train"]][["x"]])

  prior <- model[["priors"]][["a"]]
  expect_equal(dim(prior[["mean"]]), c(n_design, k))
  expect_equal(dim(prior[["cov"]]), c(n_design, n_design))
  # The names of a sampler's prior are deliberately absent: a file that carried
  # them would be read as having no coefficient prior at all.
  expect_null(prior[["mu"]])
  expect_null(prior[["v_inv"]])

  # v_i is a precision, and the unrestricted constant is the last column of the
  # design.
  expect_equal(diag(prior[["cov"]])[1], 1 / 4)
  expect_equal(diag(prior[["cov"]])[n_design], 10)

  expect_equal(model[["priors"]][["u_sigma"]][["df"]], k)
  expect_equal(dim(model[["priors"]][["u_sigma"]][["scale"]]), c(k, k))

  # The space is not drawn, so there is no prior over it.
  expect_null(model[["priors"]][["beta"]])
})


test_that("the prior refuses what the discounted models cannot use", {

  model <- discount_vec(r = 1)
  sigma <- list(df = "k", scale = 1)

  expect_error(add_priors(model, coef = list(v_i = 1, shape = 3, rate = 1e-4),
                          sigma = sigma),
               "delta_beta")
  # A flat prior has no covariance to write into the file.
  expect_error(add_priors(model, coef = list(v_i = 0), sigma = sigma),
               "positive")
  expect_error(add_priors(model, coef = list(v_i = 1),
                          coint = list(v_i = 0, p_tau_i = 1), sigma = sigma),
               "fixed cointegration matrix")
  expect_error(add_priors(model, coef = list(v_i = 1), sigma = sigma,
                          varsel = list(inprior = 0.5)),
               "[Vv]ariable selection")
})


test_that("the cointegration matrix is fixed rather than started from", {

  model <- fitted_discount_vec(r = 1)
  k_beta <- ncol(model[["data"]][["train"]][["w"]])

  # One space, not one per period: this is the model, and it does not move.
  expect_equal(length(model[["initial"]][["beta"]]), k_beta * 1L)
  expect_null(model[["initial"]][["beta_init"]])
  # Nothing iterates, so nothing else is started.
  expect_null(model[["initial"]][["a"]])
  expect_null(model[["initial"]][["u_sigma_inv"]])

  # A space of the caller's own, which is what a grid over candidate spaces
  # varies.
  own <- matrix(c(1, -1), k_beta, 1)
  given <- add_initial_values(discount_vec_with_priors(r = 1), beta = own)
  expect_equal(as.vector(given[["initial"]][["beta"]]), as.vector(own))

  expect_error(add_initial_values(discount_vec_with_priors(r = 1),
                                  beta = matrix(1, k_beta + 1, 1)),
               "cointegration matrix")

  # A model of rank zero is a VAR in differences and conditions on nothing.
  expect_null(fitted_discount_vec(r = 0)[["initial"]][["beta"]])
})


test_that("the written file is the one the discounted models read", {

  model <- fitted_discount_vec(r = 1)
  path <- tempfile(fileext = ".h5")
  write_to_hdf5(model, filename = path)
  on.exit(unlink(path), add = TRUE)

  h5 <- hdf5r::h5file(path, mode = "r")
  on.exit(if (h5$is_valid) h5$close_all(), add = TRUE)
  objects <- hdf5r::list.objects(h5)

  expect_true(all(c("priors/a/mean", "priors/a/cov", "initial/beta",
                    "data/train/x", "data/train/w") %in% objects))
  # The filter runs against one coefficient matrix per period rather than a SUR
  # design, and refuses a file that carries nothing but the SUR matrix.
  expect_false("data/train/z" %in% objects)

  attributes <- hdf5r::h5attributes(h5[["model"]])
  expect_equal(attributes[["algorithm"]], "VecTvpDiscount")
  expect_equal(attributes[["delta_beta"]], 1)
  expect_equal(attributes[["delta_sigma"]], 1)
  expect_equal(attributes[["burnin"]], 0L)

  # n_design x k and n_design x n_design as the file sees them, which is as R
  # wrote them.
  expect_equal(h5[["priors/a/mean"]]$dims,
               c(nrow(model[["priors"]][["a"]][["mean"]]),
                 ncol(model[["priors"]][["a"]][["mean"]])))
})


test_that("a posterior of one row is read as one row", {

  # What a discounted model writes: per-period paths without coda's mcpar, and
  # a log-likelihood of a single row. as.matrix() of a dataset whose dimensions
  # hdf5r drops makes a column of it, which read that row as many draws of one
  # period.
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  h5 <- hdf5r::h5file(path, mode = "a")
  group_model <- h5$create_group("model")
  hdf5r::h5attr(group_model, "algorithm") <- "VarTvpDiscount"
  hdf5r::h5attr(group_model, "k") <- 2L
  posterior <- h5$create_group("posterior")
  group_a <- posterior$create_group("a")
  group_a[["mean"]] <- matrix(1:12 * 1.0, 6, 2)
  posterior[["loglik"]] <- matrix(seq_len(6) * -1.0, 1, 6)
  h5$close_all()

  model <- read_model_from_hdf5(path)
  expect_equal(dim(model[["posterior"]][["loglik"]]), c(1L, 6L))
  expect_equal(as.vector(model[["posterior"]][["loglik"]]), -1 * seq_len(6))
  # A path is not a chain, so it is not labelled as one.
  expect_false(inherits(model[["posterior"]][["a"]][["mean"]], "mcmc"))
  expect_equal(dim(model[["posterior"]][["a"]][["mean"]]), c(6L, 2L))
})


test_that("a discounted VAR is estimated in R, as a posterior and not a chain", {

  model <- add_initial_values(add_priors(discount_var(delta_beta = 0.99, delta_sigma = 0.98),
                                         coef = list(v_i = 1, v_i_det = 1 / 10),
                                         sigma = list(df = "k", scale = 1)))

  fitted <- add_posterior_coefficients(model)

  expect_s3_class(fitted, "bvarmodel")

  tt <- nrow(fitted[["data"]][["train"]][["y"]])
  k <- fitted[["model"]][["k"]]
  n_x <- ncol(fitted[["data"]][["train"]][["x"]])

  posterior <- fitted[["posterior"]]
  expect_equal(dim(posterior[["a"]][["mean"]]), c(tt, n_x * k))
  expect_equal(dim(posterior[["a"]][["scale"]]), c(tt, n_x * k))
  expect_equal(dim(posterior[["a"]][["cov"]]), c(tt, n_x^2))
  expect_equal(dim(posterior[["u_sigma"]][["scale"]]), c(tt, k^2))
  expect_equal(dim(posterior[["df"]]), c(tt, 1L))
  expect_true(all(is.finite(posterior[["a"]][["mean"]])))

  # A period is not a draw, so none of it is labelled as one.
  expect_null(posterior[["a"]][["coeffs"]])
  expect_null(posterior[["u_sigma_inv"]])
  expect_false(inherits(posterior[["a"]][["mean"]], "mcmc"))

  # One row, the exact pointwise log marginal likelihood of each period.
  scored <- add_posterior_loglik(fitted)
  expect_equal(dim(scored[["posterior"]][["loglik"]]), c(1L, tt))
  expect_true(all(is.finite(scored[["posterior"]][["loglik"]])))
  expect_false(inherits(scored[["posterior"]][["loglik"]], "mcmc"))

  # LML is that row summed, and it is the only criterion the model carries.
  criteria <- selection_criteria(scored)
  expect_equal(criteria[["LML"]][["mean"]], sum(scored[["posterior"]][["loglik"]]))
  expect_null(criteria[["WAIC"]])
})


test_that("estimating a discounted model consumes no random numbers", {

  model <- add_initial_values(add_priors(discount_var(delta_beta = 0.98),
                                         coef = list(v_i = 1, v_i_det = 1 / 10),
                                         sigma = list(df = "k", scale = 1)))

  # Two runs from different states of the generator, which a sampler would not
  # survive: the posterior is closed form, so it is the same arithmetic twice.
  set.seed(1)
  first <- add_posterior_coefficients(model)
  set.seed(2)
  second <- add_posterior_coefficients(model)

  expect_identical(first[["posterior"]], second[["posterior"]])
})


test_that("a discounted VEC estimates against the space it conditions on", {

  model <- fitted_discount_vec(r = 1, delta_beta = 0.98)
  fitted <- add_posterior_coefficients(model)

  expect_s3_class(fitted, "bvecmodel")

  tt <- nrow(fitted[["data"]][["train"]][["y"]])
  k <- fitted[["model"]][["k"]]
  n_design <- fitted[["model"]][["rank"]] + ncol(fitted[["data"]][["train"]][["x"]])

  posterior <- fitted[["posterior"]]
  expect_equal(dim(posterior[["a"]][["mean"]]), c(tt, n_design * k))
  expect_equal(dim(posterior[["a"]][["cov"]]), c(tt, n_design^2))

  # The loadings mean nothing without the space they load on, so it travels
  # with them -- one row, because it did not move and was not estimated.
  expect_equal(as.vector(posterior[["beta"]][["coeffs"]]),
               as.vector(fitted[["initial"]][["beta"]]))

  expect_equal(dim(add_posterior_loglik(fitted)[["posterior"]][["loglik"]]), c(1L, tt))

  # Rank zero is a VAR in differences, which the model accepts and which has no
  # space to carry.
  rank_zero <- add_posterior_coefficients(fitted_discount_vec(r = 0))
  expect_null(rank_zero[["posterior"]][["beta"]])
})


test_that("a discounted model forecasts and scores what the horizon realised", {

  data("e1", package = "bvartools", envir = environment())
  series <- diff(log(e1)) * 100
  train <- stats::window(series, end = stats::time(series)[nrow(series) - 4])

  model <- add_initial_values(add_priors(
    create_bvarmodel(train, p = 1, deterministic = "const", algorithm = "discount",
                     delta_beta = 0.99, iterations = 15, burnin = 0, thin = 1),
    coef = list(v_i = 1, v_i_det = 1 / 10),
    sigma = list(df = "k", scale = 1)))

  fitted <- add_forecast_input(add_posterior_coefficients(model), n_ahead = 4)
  forecast <- add_posterior_forecasts(fitted)

  k <- fitted[["model"]][["k"]]
  # 'iterations' has lost its chain and kept its name: it is how many i.i.d.
  # draws from the closed form a forecast takes.
  expect_equal(dim(forecast[["posterior"]][["forecast"]][["forecasts"]]), c(15L, 4L * k))
  expect_true(all(is.finite(forecast[["posterior"]][["forecast"]][["forecasts"]])))

  scored <- add_predictive_loglik(forecast, test_sample = series)
  loglik <- scored[["posterior"]][["forecast"]][["loglik"]]
  # One row: the filter carries itself through the realised values, so there is
  # nothing to average over.
  expect_equal(nrow(loglik), 1L)
  expect_true(all(is.finite(loglik)))
})


test_that("the criterion of a discounted model is its marginal likelihood", {

  model <- fitted_discount_vec(r = 1)
  tt <- nrow(model[["data"]][["train"]][["y"]])
  model[["posterior"]] <- list(loglik = matrix(-2, 1, tt))

  criteria <- selection_criteria(model)
  expect_s3_class(criteria, "selcrit")
  expect_equal(criteria[["LML"]][["mean"]], -2 * tt)
  # None of the criteria that need a chain, and no LL: the marginal likelihood
  # integrates the parameters out where LL conditions on them.
  expect_null(criteria[["LL"]])
  expect_null(criteria[["WAIC"]])
  expect_null(criteria[["BIC"]])

  # Maximised, like the other two densities.
  worse <- model
  worse[["posterior"]] <- list(loglik = matrix(-3, 1, tt))
  both <- list(selection_criteria(worse), criteria)
  class(both) <- c("selcritlist", class(both))
  expect_equal(choose_best_model(both, criterion = "LML"), 2L)
})


# Runs only where BAYESTS_EXECUTABLE names a BayesTS build, and
# BAYESTS_LIBRARY_PATH, if needed, the directories of its runtime libraries.
# Everything above checks the file the R side writes; this checks that BayesTS
# agrees, which is the only thing that can.
bayests_files_for_tests <- function() {
  executable <- Sys.getenv("BAYESTS_EXECUTABLE")
  skip_if(!nzchar(executable) || !file.exists(executable),
          "BAYESTS_EXECUTABLE does not name a BayesTS executable")
  library_path <- Sys.getenv("BAYESTS_LIBRARY_PATH")
  library_path <- if (nzchar(library_path)) {
    strsplit(library_path, .Platform$path.sep, fixed = TRUE)[[1]]
  }
  bayests_files(executable = executable, library_path = library_path)
}


test_that("BayesTS estimates a grid of ranks and scores it exactly", {

  run <- bayests_files_for_tests()

  folder <- file.path(tempdir(), "bvartools-discount-grid")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  models <- discount_vec(r = 0:2, delta_beta = c(0.98, 1))
  models <- add_priors(models, coef = list(v_i = 1, v_i_det = 1 / 10),
                       sigma = list(df = "k", scale = 1))
  models <- add_initial_values(models)
  write_to_hdf5(models, folder = folder)

  stored <- open_models(folder)
  expect_equal(nrow(stored[["manifest"]]), 6L)
  expect_equal(sort(unique(stored[["manifest"]][["rank"]])), 0:2)
  expect_setequal(stored[["manifest"]][["delta_beta"]], c(0.98, 1))

  # A failed run raises; what comes back is the paths it was given.
  expect_length(run(model_files(stored)), 6L)

  criteria <- selection_criteria(stored)
  lml <- vapply(criteria, function(x) x[["LML"]][["mean"]], numeric(1))
  expect_true(all(is.finite(lml)))

  # One model read whole: the closed form is one column per period, and the
  # score is one row over the sample rather than one row per draw.
  model <- read_model_from_hdf5(file.path(folder, stored[["manifest"]][6, "file"]))
  tt <- nrow(model[["data"]][["train"]][["y"]])
  k <- model[["model"]][["k"]]
  rank <- model[["model"]][["rank"]]
  n_design <- rank + ncol(model[["data"]][["train"]][["x"]])

  posterior <- model[["posterior"]]
  expect_equal(dim(posterior[["a"]][["mean"]]), c(tt, n_design * k))
  expect_equal(dim(posterior[["a"]][["cov"]]), c(tt, n_design^2))
  expect_equal(dim(posterior[["u_sigma"]][["scale"]]), c(tt, k^2))
  expect_equal(dim(posterior[["df"]]), c(tt, 1L))
  expect_equal(dim(posterior[["loglik"]]), c(1L, tt))
  # No draws of the coefficients: joining one draw per period would look like a
  # sampled path and is not one.
  expect_null(posterior[["a"]][["coeffs"]])
  # The space the run conditioned on travels with the loadings, which mean
  # nothing without it.
  expect_equal(as.vector(posterior[["beta"]][["coeffs"]]),
               as.vector(model[["initial"]][["beta"]]))

  # LML is that row summed, which is the exact log marginal likelihood.
  expect_equal(criteria[[6]][["LML"]][["mean"]], sum(posterior[["loglik"]]))
})


test_that("the two exported helpers say what a discounted model is", {

  expect_true(is_discount_model(discount_vec()))
  expect_true(is_discount_model(discount_vec()[["model"]]))
  expect_true(is_discount_model("VarTvpDiscount"))
  expect_false(is_discount_model("VecNormalWishart"))
  expect_false(is_discount_model(NULL))

  expect_null(check_discount_specification(k = 3, delta_sigma = 0.9))
  expect_error(check_discount_specification(k = 3, error = "sv"), "wishart")
  expect_error(check_discount_specification(k = 3, delta_sigma = 0.5), "too small")
})
