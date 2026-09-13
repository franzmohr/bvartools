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

# --- structural VAR fixture -------------------------------------------------

# A recursively identified model: A_0 is lower triangular with a unit diagonal
# and the error covariance is diagonal, which is what `structural = TRUE`
# together with a gamma prior on the error variances produces.
fx_svar_fitted <- function() {
  cached_fixture("svar_fitted", {
    model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                              structural = TRUE, error = "gamma",
                              iterations = fx_iterations, burnin = fx_burnin)
    model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
                        sigma = list(shape = 1e-6, rate = 1e-6))
    model <- add_initial_values(model)
    set.seed(202401)
    add_posterior_coefficients(model)
  })
}

# --- sign restricted fixture ------------------------------------------------

# The restrictions the fixture below imposes: a shock named after investment
# that raises investment and consumption on impact. Loose enough that every
# draw is identified, which keeps the fixture about the identification rather
# than about the acceptance rate.
fx_sign_restrictions <- function() {
  data.frame(impulse = "invest",
             response = c("invest", "cons"),
             sign = c(1, 1),
             horizon = 0)
}

fx_var_sign <- function() {
  cached_fixture("var_sign", {
    set.seed(20240301)
    add_sign_restrictions(fx_var_fitted(), fx_sign_restrictions())
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

# --- time varying fixtures ----------------------------------------------------

# The prior of a stochastic volatility or gamma error term of a time varying
# model.
tvp_sigma_prior <- function(error) {
  switch(error,
         sv = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                   state_variance = 0.05, offset = 1e-8),
         list(shape = 3, rate = 0.01))
}

# A VAR with time varying parameters, by error specification. Under "sv" the
# posterior holds one error precision per period, under "gamma" a single one --
# the pair that summary(), plot() and the draws of the impulse responses have
# to tell apart.
fx_var_tvp_fitted <- function(error) {
  cached_fixture(paste0("var_tvp_", error), {
    model <- create_bvarmodel(var_data(), p = 1, deterministic = "const", tvp = TRUE,
                              error = error, iterations = fx_iterations, burnin = fx_burnin)
    model <- add_priors(model, coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4),
                        sigma = tvp_sigma_prior(error))
    model <- add_initial_values(model)
    set.seed(314271)
    model <- add_posterior_coefficients(model)
    add_posterior_loglik(model)
  })
}

# The same for a VEC model, with rho drawn, so that the posterior also holds a
# block of draws that is neither a path over the periods nor a coefficient.
fx_vec_tvp_fitted <- function(error) {
  cached_fixture(paste0("vec_tvp_", error), {
    model <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "unrestricted", tvp = TRUE,
                              error = error, iterations = fx_iterations, burnin = fx_burnin)
    model <- add_priors(model, coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4),
                        coint = list(rho = 0.99, rho_min = 0.9, rho_max = 0.999),
                        sigma = tvp_sigma_prior(error))
    model <- add_initial_values(model)
    set.seed(271828)
    model <- add_posterior_coefficients(model)
    add_posterior_loglik(model)
  })
}

# Every block of draws of a posterior, flattened to "block$element" names.
draw_blocks <- function(posterior, prefix = NULL) {
  blocks <- list()
  for (name in names(posterior)) {
    element <- posterior[[name]]
    label <- paste(c(prefix, name), collapse = "$")
    if (is.list(element) && !inherits(element, "mcmc")) {
      blocks <- c(blocks, draw_blocks(element, label))
    } else if (!is.null(element)) {
      blocks[[label]] <- element
    }
  }
  blocks
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

# --- at_macrodata fixtures ----------------------------------------------------

# Austrian output, inflation and short-term interest rate in levels, from the
# domestic series of the at_macrodata data set that new tests are written
# against. The data set is either one matrix of all series or a list of its
# domestic ('endogen') and foreign ('exogen') series; both are taken.
at_data <- function() {
  data <- bvartools::at_macrodata
  if (is.list(data) && !stats::is.ts(data)) {
    data <- data[["endogen"]]
  }
  data[, c("y", "Dp", "r")]
}

# A VEC model with time varying parameters and stochastic volatility on
# at_data(), with rho drawn and the cointegration space centred on the ML
# estimate, and its log-likelihood: a posterior with every kind of block and a
# prior with every kind of element, scalars and characters among them.
fx_at_vec_tvp <- function() {
  cached_fixture("at_vec_tvp", {
    model <- create_bvecmodel(at_data(), p = 2, r = 1, const = "unrestricted", tvp = TRUE,
                              error = "sv", iterations = fx_iterations, burnin = fx_burnin)
    # The ML prior may be floored by rho, which is not what these tests are about.
    model <- suppressWarnings(
      add_priors(model, coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4),
                 coint = list(rho = 0.99, rho_min = 0.9, rho_max = 0.999,
                              p_tau_i = "ml", weight = 0.1),
                 sigma = tvp_sigma_prior("sv")))
    model <- add_initial_values(model)
    set.seed(179)
    model <- add_posterior_coefficients(model)
    add_posterior_loglik(model)
  })
}

# A VAR model with constant coefficients and a gamma error term with covariance
# block on at_data(), and its log-likelihood.
fx_at_var <- function() {
  cached_fixture("at_var", {
    model <- create_bvarmodel(at_data(), p = 2, deterministic = "const",
                              error = "gamma+covar", iterations = fx_iterations,
                              burnin = fx_burnin)
    model <- add_priors(model, coef = list(v_i = 1, v_i_det = 0.1),
                        sigma = list(shape = 3, rate = 0.01))
    model <- add_initial_values(model)
    set.seed(180)
    model <- add_posterior_coefficients(model)
    add_posterior_loglik(model)
  })
}
