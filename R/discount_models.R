# The two discounted models, which are not samplers.
#
# VarTvpDiscount and VecTvpDiscount are the matrix normal dynamic linear model
# of West & Harrison (1997, ch. 16) with the discounted Wishart of Uhlig (1997).
# Both posteriors are closed form: one pass over the sample, no chain, and no
# random numbers consumed. That is why so much of the rest of this package has
# to be told about them by name rather than by the `tvp` flag they share with
# the time varying samplers.
#
# What differs, and where each difference is dealt with:
#
#   the file's coefficient prior is /priors/a/mean and /priors/a/cov, a matrix
#     normal, rather than the /priors/a/mu and /priors/a/v_inv of a vectorised
#     normal                                       -- add_priors.bvarmodel(),
#                                                     add_priors.bvecmodel()
#   the regressors are the compact /data/train/x, and a file carrying only the
#     SUR matrix /data/train/z is refused           -- write_to_hdf5.*()
#   `burnin` must be zero and `thin` one            -- the checks below
#   two discounts live at /model/delta_beta and /model/delta_sigma
#                                                   -- create_bvarmodel(),
#                                                      create_bvecmodel()
#   a VEC conditions on a fixed cointegration matrix at /initial/beta rather
#     than drawing one                              -- add_initial_values.bvecmodel()
#   the posterior is per period and has no /posterior/a/coeffs
#                                                   -- read_model_from_hdf5()
#
# Estimation itself is BayesTS's. Nothing here runs the filter: the R side
# writes the file, the command line reads it, and what comes back is read by
# read_model_from_hdf5().

# The algorithms this file is about. One list, so that a function that has to
# ask names it once.
.discount_algorithms <- c("VarTvpDiscount", "VecTvpDiscount")

#' The Discounted Models
#'
#' Whether a model is one of the two discounted models, and what a discounted
#' model cannot be. Both are exported for the packages that build a model of
#' their own on this one and have to make the same distinction.
#'
#' @param x a model, a specification list -- the \code{model} element of one --
#' or the name of an algorithm.
#' @param k the number of endogenous variables.
#' @param error,varsel,structural,burnin,thin the corresponding arguments of
#' \code{\link{create_bvarmodel}} or \code{\link{create_bvecmodel}}.
#' @param delta_beta,delta_sigma the two discount factors.
#'
#' @details
#' \code{VarTvpDiscount} and \code{VecTvpDiscount} are the matrix normal
#' dynamic linear model of West & Harrison (1997, ch. 16) with the discounted
#' Wishart of Uhlig (1997). Their posterior is closed form -- one pass over the
#' sample, no chain, no random numbers consumed -- which is what every refusal
#' below follows from. See \code{\link{create_bvecmodel}} for what they are and
#' how they are set up.
#'
#' \code{check_discount_specification} raises an error for a specification a
#' discounted model cannot carry and returns nothing otherwise. Each refusal is
#' one BayesTS would raise as well, against a file that had already been
#' written; raising it here means one error for a grid of models rather than one
#' per file.
#'
#' @return \code{is_discount_model} returns a single logical.
#' \code{check_discount_specification} returns \code{NULL} invisibly.
#'
#' @references
#'
#' Uhlig, H. (1997). Bayesian vector autoregressions with stochastic volatility.
#' \emph{Econometrica, 65}(1), 59--73. \doi{10.2307/2171813}
#'
#' West, M., & Harrison, J. (1997). \emph{Bayesian forecasting and dynamic models}
#' (2nd ed.). New York: Springer.
#'
#' @examples
#'
#' data("e6")
#' model <- create_bvecmodel(e6 * 100, p = 2, r = 1, const = "unrestricted",
#'                           algorithm = "discount", delta_beta = 0.98,
#'                           iterations = 10, burnin = 0, thin = 1)
#' is_discount_model(model)
#' is_discount_model("VecNormalWishart")
#'
#' @name discount_models
#' @export
is_discount_model <- function(x) {

  if (is.null(x)) {
    return(FALSE)
  }
  algorithm <- if (is.list(x) && !is.null(x[["model"]])) {
    x[["model"]][["algorithm"]]
  } else if (is.list(x)) {
    x[["algorithm"]]
  } else {
    x
  }

  isTRUE(algorithm %in% .discount_algorithms)
}

# The short spelling the rest of the package uses.
.is_discount <- function(x) {
  is_discount_model(x)
}

# The two discounts, checked.
#
# `delta_beta` governs the coefficients and `delta_sigma` the error covariance.
# One is not a neutral default but a model: at `delta_beta = 1` the coefficients
# do not move, and at `delta_sigma = 1` neither does the error covariance. The
# range is therefore half open at that end.
#
# `k` is the number of endogenous variables, and is checked against
# `delta_sigma` because the degrees of freedom of the discounted Wishart settle
# at 1 / (1 - delta_sigma) whatever the prior asks for. Below `k` the inverse
# Wishart is improper and every scale read off it is meaningless, which BayesTS
# refuses rather than warns about -- and refuses on the file, after it has been
# written. Refused here too, where the specification is still in hand.
.check_discount_deltas <- function(delta_beta, delta_sigma, k) {

  for (name in c("delta_beta", "delta_sigma")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) == 0 || anyNA(value)) {
      stop("Argument '", name, "' must be a numeric vector without missing values.")
    }
    if (any(value <= 0 | value > 1)) {
      stop("Argument '", name, "' must lie in (0, 1]. One is the model in which ",
           if (name == "delta_beta") "the coefficients" else "the error covariance",
           " does not move.")
    }
  }

  too_small <- delta_sigma < 1 & 1 / (1 - delta_sigma) < k
  if (any(too_small)) {
    stop("Argument 'delta_sigma' is too small for ", k, " endogenous variables: ",
         "the degrees of freedom of the discounted Wishart settle at ",
         "1 / (1 - delta_sigma), which must be at least the number of variables. ",
         "The smallest value this model admits is ", signif(1 - 1 / k, 3), ".")
  }

  return(invisible(NULL))
}

# What a discounted model cannot be, checked before anything is built.
#
# Each of these is refused by BayesTS as well. They are refused here so that the
# message arrives while the specification is still an argument rather than a
# file, which for a grid of models is the difference between one error and one
# per file.
#' @rdname discount_models
#' @export
check_discount_specification <- function(k, error = "wishart", varsel = "none",
                                         structural = FALSE, burnin = 0,
                                         thin = 1, delta_beta = 1,
                                         delta_sigma = 1) {

  .check_discount_specification(error, varsel, structural, burnin, thin)
  .check_discount_deltas(delta_beta, delta_sigma, k)

  return(invisible(NULL))
}

.check_discount_specification <- function(error, varsel, structural, burnin, thin) {

  if (!identical(error, "wishart")) {
    stop("The discounted models estimate an inverse Wishart error covariance, ",
         "so argument 'error' must be \"wishart\". The discounting of the ",
         "covariance is asked for with 'delta_sigma' instead.")
  }
  if (!identical(varsel, "none")) {
    stop("Variable selection is not available for the discounted models: they ",
         "have no draws for an inclusion indicator to be drawn alongside.")
  }
  if (isTRUE(structural)) {
    stop("A structural model cannot be estimated against an unrestricted error ",
         "covariance, which an inverse Wishart posterior is.")
  }
  if (!identical(as.integer(burnin), 0L)) {
    stop("The discounted models have no chain to burn in, so argument 'burnin' ",
         "must be 0. Argument 'iterations' keeps its meaning: it is how many ",
         "i.i.d. draws a forecast takes from the closed form.")
  }
  if (!identical(as.integer(thin), 1L)) {
    stop("The discounted models draw i.i.d. rather than sweeping, so there is ",
         "nothing to thin: argument 'thin' must be 1.")
  }

  return(invisible(NULL))
}


# The priors of a discounted model, which are a different pair of objects.
#
# A sampler here is given a normal prior over the vectorised coefficients:
# /priors/a/mu, a vector, and /priors/a/v_inv, a precision over the whole of it.
# The discounted models are given a matrix normal prior over the coefficient
# matrix instead: /priors/a/mean, which is n_design x k, and /priors/a/cov,
# which is the n_design x n_design regressor side of its covariance. The full
# prior covariance of the coefficients is Sigma kronecker cov, so `cov` is in
# units of the error covariance, exactly as the natural conjugate prior of a
# constant coefficient Wishart model is -- and, like it, this is what makes the
# posterior closed form.
#
# The names are deliberately not the sampler's. A file that brought `mu` and
# `v_inv` along would be read as having no coefficient prior at all, and the
# model would run against the default one without saying so.
#
# `coef` is read as it is elsewhere: `v_i` is a prior precision, so the diagonal
# of `cov` is its reciprocal, `v_i_det` applies to the unrestricted
# deterministic terms and `v_i_alpha` to the loadings of a VEC. `const` centres
# the intercept. What is not read is `shape` and `rate`: the state variances
# they are the prior of do not exist here, the drift being governed by
# `delta_beta` instead, and a model that took them would be a different model
# from the one the file describes.
.add_priors_discount <- function(object, coef, sigma) {

  allowed <- c("v_i", "v_i_det", "v_i_alpha", "const")
  for (i in names(coef)) {
    if (!i %in% allowed) {
      stop("Element '", i, "' in argument 'coef' is not used by the discounted ",
           "models. They take ", paste0("'", allowed, "'", collapse = ", "),
           ". The drift of the coefficients is governed by 'delta_beta' of ",
           "create_bvarmodel() or create_bvecmodel(), not by a prior on a state ",
           "variance, so 'shape' and 'rate' have nothing to be the prior of.")
    }
  }

  if (is.null(coef[["v_i"]]) || !is.numeric(coef[["v_i"]]) ||
      length(coef[["v_i"]]) != 1 || is.na(coef[["v_i"]]) || coef[["v_i"]] <= 0) {
    stop("Argument 'coef$v_i' must be a single positive number for the ",
         "discounted models. It is a prior precision and the prior it belongs ",
         "to is stated as a covariance, so a flat prior -- 'v_i = 0', which the ",
         "samplers accept -- has no covariance to be written into the file.")
  }
  for (name in c("v_i_det", "v_i_alpha")) {
    value <- coef[[name]]
    if (!is.null(value) &&
        (!is.numeric(value) || length(value) != 1 || is.na(value) || value <= 0)) {
      stop("Argument 'coef$", name, "' must be a single positive number.")
    }
  }

  k <- object[["model"]][["k"]]
  rank <- object[["model"]][["rank"]]
  if (is.null(rank)) {
    rank <- 0L
  }
  x <- object[["data"]][["train"]][["x"]]
  n_x <- if (is.null(x)) 0L else NCOL(x)
  n_design <- rank + n_x
  if (n_design == 0) {
    stop("The model has no regressors at all, so there is nothing for a ",
         "coefficient prior to be the prior of.")
  }

  # The mean: one row per column of the design, one column per equation, and
  # the design is the `rank` error correction columns in front of the compact
  # regressors -- the layout BayesTS builds and validate() checks against.
  mean <- matrix(0, n_design, k)
  if (!is.null(coef[["const"]]) && n_x > 0) {
    position <- which(dimnames(x)[[2]] == "const")
    if (length(position) == 1) {
      position <- position + rank
      value <- coef[["const"]]
      if (is.character(value)) {
        y <- object[["data"]][["train"]][["y"]]
        value <- switch(value,
                        "first" = y[1, ],
                        "mean" = colMeans(y),
                        stop("Invalid specification of coef$const."))
      }
      if (!length(value) %in% c(1, k)) {
        stop("When a numeric is provided in argument 'coef$const', it must be ",
             "either a single number or a vector of the same length as the ",
             "number of endogenous variables in the model.")
      }
      mean[position, ] <- value
    }
  }

  # The covariance. `v_i` is a precision, as everywhere else in add_priors(),
  # so the diagonal is its reciprocal.
  cov <- diag(1 / coef[["v_i"]], n_design)
  if (rank > 0 && !is.null(coef[["v_i_alpha"]])) {
    diag(cov)[seq_len(rank)] <- 1 / coef[["v_i_alpha"]]
  }
  n_det <- object[["model"]][["n"]]
  if (!is.null(n_det) && n_det > 0 && !is.null(coef[["v_i_det"]])) {
    diag(cov)[n_design - n_det + seq_len(n_det)] <- 1 / coef[["v_i_det"]]
  }

  object[["priors"]][["a"]] <- list(type = "matrixnormal",
                                    mean = mean,
                                    cov = cov)

  # The error covariance. An inverse Wishart, whatever `delta_sigma` is: the
  # discount says how fast the sample is forgotten, not what family the
  # posterior belongs to.
  for (i in names(sigma)) {
    if (!i %in% c("df", "scale")) {
      stop("Element '", i, "' in argument 'sigma' is not used by the discounted ",
           "models. They take 'df' and 'scale', the inverse Wishart prior of ",
           "the error covariance.")
    }
  }
  if (is.null(sigma[["df"]]) || is.null(sigma[["scale"]])) {
    stop("Arguments 'sigma$df' and 'sigma$scale' must both be specified.")
  }
  df <- sigma[["df"]]
  if (is.character(df)) {
    if (!grepl("^[k0-9 +*/().-]+$", df)) {
      stop("Use no other letter than 'k' in 'sigma$df' to indicate the number ",
           "of endogenous variables.")
    }
    df <- eval(parse(text = df), list(k = k))
  }
  if (any(df <= 0)) {
    stop("Current specification implies non-positive prior degrees of freedom ",
         "of the error term. 'sigma$df' must be positive.")
  }

  object[["priors"]][["u_sigma"]] <- list(type = "wishart",
                                          df = df,
                                          scale = diag(sigma[["scale"]], k))

  return(object)
}


# The starting values of a discounted model, of which there is at most one.
#
# Nothing iterates, so there is nowhere for a starting value to be the start of,
# and /initial is empty for a discounted VAR and for a discounted VEC of rank
# zero. What a discounted VEC of positive rank needs is not a starting value at
# all: it is the cointegration matrix the run conditions on, which is read from
# /initial/beta because that is where every VEC in this package keeps its space
# and in the same layout -- vec of a k_beta x rank matrix. Reusing the path is
# what lets a file written for a sampling VEC be pointed at this model by
# changing the algorithm and adding the two discounts.
#
# The default space is Johansen's maximum likelihood estimate, which is what
# `method = "maxlik"` means here as it does elsewhere. `beta` overrides it with
# a space of the caller's own -- a vector from theory, a restriction, or the
# posterior mean of a constant coefficient VEC -- and a grid over candidate
# spaces is a grid over this argument.
#
# The seed is still drawn. Estimation consumes no random numbers, so two runs of
# it agree to the bit whatever the seed is, but a forecast takes i.i.d. draws
# from the closed form and those are a chain's worth of random numbers.
.add_initial_values_discount <- function(object, method, beta = NULL, ...) {

  rank <- object[["model"]][["rank"]]
  if (is.null(rank)) {
    rank <- 0L
  }

  if (rank > 0) {

    k_beta <- NCOL(object[["data"]][["train"]][["w"]])

    if (is.null(beta)) {
      if (!method %in% c("maxlik", "ols")) {
        stop("The discounted VEC conditions on a cointegration matrix rather ",
             "than drawing one, so there is no prior to draw a starting value ",
             "from: use method = \"maxlik\" for Johansen's estimate, or pass the ",
             "space in argument 'beta'.")
      }
      beta <- .coint_ml(object)[["beta"]]
    }

    beta <- as.matrix(beta)
    if (nrow(beta) != k_beta || ncol(beta) != rank) {
      stop("Argument 'beta' must be the ", k_beta, " x ", rank,
           " cointegration matrix the model conditions on, and is ",
           nrow(beta), " x ", ncol(beta), ".")
    }
    if (!all(is.finite(beta))) {
      stop("Argument 'beta' holds a value that is not finite.")
    }

    # Not put on the scale of a state equation, which is what
    # .tvp_initial_beta_scale() does for the sampling time varying VECs. There
    # is no state equation for the space here: it does not move, so there is no
    # stationary variance for its norm to match, and rescaling it would only
    # rescale the loadings the filter estimates against it.
    object[["initial"]][["beta"]] <- matrix(beta)
  }

  if (is.null(object[["model"]][["seed"]])) {
    object[["model"]][["seed"]] <- .draw_model_seed()
  }

  return(object)
}


# The four stages of a discounted model, each of which runs the filter in the
# vendored BayesTS core rather than a sampler.
#
# They are separate from the switch() of each step because the steps around them
# do not apply. A sampler's posterior is a chain, so every step wraps what comes
# back in coda's mcpar -- read off `u_sigma_inv`, which these do not have -- and
# a discounted posterior is one column per period of a closed form, with no
# start, no end and nothing thinned. Each step therefore hands the model over
# here and returns what comes back, rather than falling through to machinery
# that would label periods as draws. The forecasts are the exception: they are
# draws, and `.discount_forecasts()` labels them as such itself.
#
# Nothing in `.discount_coefficients()` or the VAR's log-likelihood consumes the
# random number generator: the posterior is closed form, so two runs agree to
# the bit and the seed has nothing to repeat. Forecasting does consume it, being
# i.i.d. draws from that posterior, as does scoring a VEC forecast, and both are
# therefore run under the model's seed like every sampler here.

.discount_coefficients <- function(object) {

  # The entry points replace the posterior of the object they were given rather
  # than rebuilding it from named elements, as the samplers' do: a discounted
  # model need not have an `initial` at all -- nothing iterates -- and naming an
  # element that is absent is an error. What does not survive the round trip is
  # the class, as it does not for any of the C++ entry points here, so it is put
  # back by hand.
  .with_discount_class(object, switch(object[["model"]][["algorithm"]],
                                      VarTvpDiscount = .VarTvpDiscountCoefficients(object),
                                      VecTvpDiscount = .VecTvpDiscountCoefficients(object)))
}

# The class of the model the C++ side was given, put back on what it returned.
.with_discount_class <- function(object, result) {
  class(result) <- class(object)
  result
}

.discount_loglik <- function(object) {

  if (is.null(object[["posterior"]][["a"]][["mean"]])) {
    stop("Object does not contain an estimated posterior in posterior$a$mean. ",
         "Use add_posterior_coefficients() first.")
  }

  # Not wrapped in an 'mcmc': one row, the exact pointwise log marginal
  # likelihood of each period, rather than one row per draw.
  .with_discount_class(object, switch(object[["model"]][["algorithm"]],
                                      VarTvpDiscount = .VarTvpDiscountLogLik(object),
                                      VecTvpDiscount = .VecTvpDiscountLogLik(object)))
}

.discount_forecasts <- function(object) {

  if (is.null(object[["model"]][["h"]])) {
    stop("Model specification does not contain forecast horizon 'h'. Consider using ",
         "function add_forecast_input().")
  }

  object <- .with_discount_class(object,
                                 .with_model_seed(object[["model"]][["seed"]],
                                                  switch(object[["model"]][["algorithm"]],
                                                         VarTvpDiscount = .VarTvpDiscountForecasts(object),
                                                         VecTvpDiscount = .VecTvpDiscountForecasts(object))))

  # The one part of a discounted posterior that is draws, and so the one part
  # labelled as them: i.i.d. from the closed form, one row each. Everything
  # downstream of a forecast -- add_forecast_errors(), the writer's
  # start/end/thin, predict() -- reads coda's mcpar off it, as it does off a
  # sampler's, and a plain matrix has none. There being no chain, the draws are
  # numbered one to S and nothing is thinned, which is what BayesTS writes for
  # them too.
  forecasts <- object[["posterior"]][["forecast"]][["forecasts"]]
  object[["posterior"]][["forecast"]][["forecasts"]] <- coda::mcmc(forecasts, start = 1,
                                                                  end = nrow(forecasts),
                                                                  thin = 1)

  object
}

.discount_score <- function(object) {

  .with_discount_class(object,
                       .with_model_seed(object[["model"]][["seed"]],
                                        switch(object[["model"]][["algorithm"]],
                                               VarTvpDiscount = .VarTvpDiscountScore(object),
                                               VecTvpDiscount = .VecTvpDiscountScore(object))))
}


# The criteria of a discounted model, of which there is one worth having.
#
# The samplers are compared by criteria that estimate what a model would say
# about data it has not seen, out of a chain: WAIC and LOO penalise by the
# flexibility the fit used, and BIC by a count of parameters standing in for it.
# None of them applies here. There is no chain to estimate an effective number
# of parameters from, and a count of parameters would be the wrong penalty for a
# model whose coefficients drift.
#
# What the file carries instead is better than any of them. /posterior/loglik is
# one row -- the one step ahead predictive density of each observation under the
# parameters integrated out exactly -- so its sum is the log marginal likelihood
# of the sample, not an estimate of it. Two discounted models of the same data
# are compared by it directly, whatever they differ in: the rank, the
# cointegration matrix, the lag order or the two discounts.
#
# It is reported as LML rather than as LL, because it is not the log likelihood
# the samplers report and the two are not comparable: LL conditions on the
# parameters and LML integrates them out, so LML is the smaller of the two for
# the same model and the difference is the complexity it paid for.
.selection_criteria_discount <- function(object, ci) {

  if (ci < 0 || ci > 1) {
    stop("Argument 'ci' is not within the permitted range of 0 and 1.")
  }
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low

  loglik <- object[["posterior"]][["loglik"]]
  predictive <- .model_predictive_densities(object)

  if (is.null(loglik) && is.null(predictive)) {
    stop("The model carries neither /posterior/loglik nor a scored forecast. ",
         "Run `bayests posterior <file.h5>` over it, or `bayests loglik` alone.")
  }

  result <- NULL
  result[["model"]] <- object[["model"]]

  if (!is.null(loglik)) {
    loglik <- as.matrix(loglik)
    if (nrow(loglik) != 1) {
      stop("The log-likelihood of a discounted model is one row, the exact ",
           "pointwise log marginal likelihood, and this one has ", nrow(loglik),
           ". The file was not written by one of the discounted models.")
    }
    result[["LML"]] <- .point_criterion(sum(loglik))
  }

  if (!is.null(predictive)) {
    result[["LPL"]] <- .lpl_entry(predictive[["densities"]],
                                  predictive[["periods"]], ci_low, ci_high)
  }

  attr(result, "ci") <- c(paste0(ci_low * 100, "%"), paste0(ci_high * 100, "%"))
  class(result) <- c("selcrit", class(object))

  return(result)
}
