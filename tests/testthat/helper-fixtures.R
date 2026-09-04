# Shared fixtures.
#
# Posterior simulation is the expensive part of every workflow test, so the
# fitted models are built once per test run and reused across files. Everything
# is deliberately tiny: one lag, few draws. The tests check structure and
# invariants, not the statistical quality of a sampler that would need thousands
# of draws to judge.

.fixtures <- new.env(parent = emptyenv())

# Evaluate `expr` the first time `name` is requested, then serve it from cache.
cached_fixture <- function(name, expr) {
  if (!exists(name, envir = .fixtures, inherits = FALSE)) {
    assign(name, expr, envir = .fixtures)
  }
  get(name, envir = .fixtures, inherits = FALSE)
}

# Number of draws kept by the fixture samplers.
fx_iterations <- 30L
fx_burnin <- 5L

# West German investment/income/consumption, log-differences in percent.
var_data <- function() {
  stats::window(diff(log(bvartools::e1)) * 100, end = c(1978, 4))
}

# Danish interest rate and inflation, in percent.
vec_data <- function() {
  bvartools::e6 * 100
}

# --- VAR fixtures -----------------------------------------------------------

fx_var_model <- function() {
  cached_fixture("var_model", create_bvarmodel(
    var_data(), p = 1, deterministic = "const",
    iterations = fx_iterations, burnin = fx_burnin
  ))
}

fx_var_priors <- function() {
  cached_fixture("var_priors", add_priors(
    fx_var_model(),
    coef = list(v_i = 0, v_i_det = 0),
    sigma = list(df = 1, scale = 0.0001)
  ))
}

fx_var_initial <- function() {
  cached_fixture("var_initial", add_initial_values(fx_var_priors()))
}

# Coefficient draws plus the log-likelihood, i.e. everything the evaluation and
# application functions need.
fx_var_fitted <- function() {
  cached_fixture("var_fitted", {
    set.seed(987654)
    object <- add_posterior_coefficients(fx_var_initial())
    add_posterior_loglik(object)
  })
}

# --- VEC fixtures -----------------------------------------------------------

fx_vec_model <- function() {
  cached_fixture("vec_model", create_bvecmodel(
    vec_data(), p = 1, r = 1, const = "unrestricted",
    iterations = fx_iterations, burnin = fx_burnin
  ))
}

fx_vec_priors <- function() {
  cached_fixture("vec_priors", add_priors(
    fx_vec_model(),
    coef = list(v_i = 1, v_i_det = 1 / 10),
    coint = list(v_i = 0, p_tau_i = 1),
    sigma = list(df = "k", scale = 1)
  ))
}

fx_vec_initial <- function() {
  cached_fixture("vec_initial", add_initial_values(fx_vec_priors()))
}

fx_vec_fitted <- function() {
  cached_fixture("vec_fitted", {
    set.seed(456789)
    object <- add_posterior_coefficients(fx_vec_initial())
    add_posterior_loglik(object)
  })
}

# --- utilities --------------------------------------------------------------

# Run plotting code without producing files or opening a device.
expect_plots <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  expect_no_error(force(code))
}

# Number of draws stored in a posterior element.
n_draws <- function(object, par = "a") {
  NROW(object[["posterior"]][[par]][["coeffs"]])
}

# A fitted list of two VAR models, used by the modellist and model selection
# tests.
fx_var_modellist <- function() {
  cached_fixture("var_modellist", {
    models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                               iterations = fx_iterations, burnin = fx_burnin)
    models <- add_priors(models, coef = list(v_i = 0, v_i_det = 0),
                         sigma = list(df = 1, scale = 0.0001))
    models <- add_initial_values(models)
    set.seed(135791)
    models <- add_posterior_coefficients(models)
    add_posterior_loglik(models)
  })
}

# A fitted VAR with simulated forecasts attached.
fx_var_forecast <- function() {
  cached_fixture("var_forecast", {
    object <- add_forecast_input(fx_var_fitted(), n_ahead = 5)
    set.seed(24680)
    add_posterior_forecasts(object)
  })
}

# A short expanding window exercise over the last few periods of the sample.
fx_expanding_window <- function() {
  cached_fixture("expanding_window", {
    model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                              iterations = 10, burnin = 5)
    model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                        sigma = list(df = 1, scale = 0.0001))
    windows <- use_expanding_window(model, start = c(1978, 1))
    windows <- add_initial_values(windows)
    set.seed(864209)
    windows <- add_posterior_coefficients(windows)
    add_posterior_loglik(windows)
  })
}

# The expanding window exercise with forecasts and forecast errors against the
# periods held back from the estimation sample.
fx_expanding_forecast <- function() {
  cached_fixture("expanding_forecast", {
    full <- var_data()
    train <- stats::window(full, end = c(1977, 4))

    model <- create_bvarmodel(train, p = 1, deterministic = "const",
                              iterations = 10, burnin = 5)
    model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                        sigma = list(df = 1, scale = 0.0001))
    windows <- use_expanding_window(model, start = c(1977, 1))
    windows <- add_initial_values(windows)
    set.seed(112358)
    windows <- add_posterior_coefficients(windows)
    windows <- add_forecast_input(windows, n_ahead = 2)
    windows <- add_posterior_forecasts(windows)
    add_forecast_errors(windows, test_sample = full)
  })
}

# Temporary paths inside the session temp directory, which R removes on exit.
temp_h5_file <- function() {
  tempfile(fileext = ".h5")
}

temp_model_dir <- function() {
  path <- tempfile()
  dir.create(path)
  path
}
