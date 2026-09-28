#' Adaptive Priors, Stationarity and the Steady State
#'
#' Adds the options of the prior that act on top of the one
#' \code{\link{add_priors}} sets: hyperparameters estimated from the data, a
#' restriction of the coefficients to a stationary model, a prior on the
#' unconditional mean in place of one on the intercept, and the prior of the
#' error of soft constraints.
#'
#' @param object an object of class 'bvarmodel' or 'modellist', usually the
#' result of \code{\link{add_priors}}.
#' @param shrinkage a character, \code{"none"} (default), \code{"minnesota"},
#' \code{"horseshoe"} or \code{"normal_gamma"}, or a list with element
#' \code{type} naming one of them and elements \code{group}, \code{shape},
#' \code{rate} and, for \code{"normal_gamma"}, \code{theta}. See 'Details'.
#' @param stationary logical. If \code{TRUE}, every draw of the coefficients
#' is a stationary model. Defaults to \code{FALSE}.
#' @param steady_state an optional list with elements \code{mu}, the prior mean
#' of the unconditional mean of the endogenous variables, and \code{v_i}, its
#' prior precision, each a number or one value per variable. See 'Details'.
#' @param constraints an optional list with elements \code{shape} and
#' \code{rate}, the gamma prior of the error precision of soft constraints, a
#' number or one value per series named in argument \code{soft} of
#' \code{\link{create_bvarmodel}}. Defaults to \code{list(shape = 3, rate = 0.01)}
#' where the model has soft constraints.
#'
#' @details Argument \code{shrinkage} lets the data decide how tightly the prior
#' of \code{\link{add_priors}} shrinks the coefficients of the lagged endogenous
#' variables. Their prior variances are multiplied by scales that are drawn with
#' the coefficients:
#' \itemize{
#'  \item \code{"minnesota"}: one scale per group of coefficients with an inverse
#'  gamma prior of \code{shape} and \code{rate}, the hierarchical Minnesota prior
#'  of Chan (2021). By default own lags form group 1 and the lags of the other
#'  variables group 2, with \code{shape = 3} and \code{rate = 2}, a prior mean of
#'  one: the prior of \code{add_priors} is the centre of the one estimated.
#'  \item \code{"horseshoe"}: the horseshoe prior of Carvalho, Polson and Scott
#'  (2010), a half-Cauchy scale per group and one per coefficient, drawn as in
#'  Makalic and Schmidt (2016). By default all lags form one group.
#'  \item \code{"normal_gamma"}: the normal-gamma prior of Griffin and Brown
#'  (2010) as Huber and Feldkircher (2019) put it on a VAR. Each coefficient's
#'  scale \eqn{\psi_j} has a gamma prior of shape \code{theta} and rate
#'  \eqn{\theta \lambda_g / 2}, and each group's \eqn{\lambda_g} a gamma
#'  prior of \code{shape} and \code{rate}. By default the lags of order
#'  \eqn{l} form group \eqn{l}, with \code{theta = 0.1} and
#'  \code{shape = rate = 0.01}. \code{theta} is held fixed; the smaller it is,
#'  the more the prior pushes small coefficients to zero while leaving large
#'  ones alone.
#' }
#' Coefficients in group 0 -- by default the deterministic terms and the
#' exogenous variables -- keep the prior of \code{add_priors}. A \code{group} of
#' one's own has one element per coefficient, in the order of the columns of
#' \code{data$train$z}. The prior must be diagonal wherever it shrinks.
#'
#' Argument \code{stationary} redraws a draw of the coefficients that is not
#' stationary up to 100 times, and keeps the previous draw if none of them is,
#' which leaves the posterior restricted to the stationary region invariant.
#'
#' Argument \code{steady_state} puts the prior on the unconditional mean
#' \eqn{\mu} of the endogenous variables rather than on the intercept, which is
#' \eqn{(I - \sum_i A_i) \mu} in every draw, following Villani (2009). A prior
#' belief about the level a series returns to is more often available than one
#' about an intercept. The model must have lags and an intercept and nothing
#' else, and cannot be combined with \code{shrinkage}. The draws of \eqn{\mu}
#' come back as \code{posterior$mu}.
#'
#' The three are available for models with constant coefficients and
#' \code{error = "wishart"}, \code{"gamma"}, \code{"gamma+covar"}, \code{"sv"} or
#' \code{"sv+covar"}, without variable selection or a structural form.
#'
#' @return The object in \code{object} with \code{model$shrinkage},
#' \code{model$stationary} and \code{model$steady_state} set where they apply, and
#' the priors in \code{priors$a$shrinkage}, \code{priors$mu} and
#' \code{priors$constraints}.
#'
#' @examples
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' model <- create_bvarmodel(e1, p = 2, iterations = 20, burnin = 10)
#' model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' model <- add_prior_options(model, shrinkage = "minnesota", stationary = TRUE)
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' @references
#'
#' Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
#' estimator for sparse signals. \emph{Biometrika, 97}(2), 465--480.
#' \doi{10.1093/biomet/asq017}
#'
#' Griffin, J. E., & Brown, P. J. (2010). Inference with normal-gamma prior
#' distributions in regression problems. \emph{Bayesian Analysis, 5}(1),
#' 171--188. \doi{10.1214/10-BA507}
#'
#' Huber, F., & Feldkircher, M. (2019). Adaptive shrinkage in Bayesian vector
#' autoregressive models. \emph{Journal of Business & Economic Statistics,
#' 37}(1), 27--39. \doi{10.1080/07350015.2016.1256217}
#'
#' Chan, J. C. C. (2021). Minnesota-type adaptive hierarchical priors for large
#' Bayesian VARs. \emph{International Journal of Forecasting, 37}(3),
#' 1212--1226. \doi{10.1016/j.ijforecast.2021.01.002}
#'
#' Makalic, E., & Schmidt, D. F. (2016). A simple sampler for the horseshoe
#' estimator. \emph{IEEE Signal Processing Letters, 23}(1), 179--182.
#' \doi{10.1109/LSP.2015.2503725}
#'
#' Villani, M. (2009). Steady-state priors for vector autoregressions.
#' \emph{Journal of Applied Econometrics, 24}(4), 630--650.
#' \doi{10.1002/jae.1065}
#'
#' @family model set-up
#' @export
add_prior_options <- function(object, shrinkage = "none", stationary = FALSE,
                              steady_state = NULL, constraints = NULL) {

  if (inherits(object, "modellist")) {
    for (i in seq_along(object)) {
      object[[i]] <- add_prior_options(object[[i]], shrinkage = shrinkage, stationary = stationary,
                                       steady_state = steady_state, constraints = constraints)
    }
    return(object)
  }
  if (!inherits(object, "bvarmodel")) {
    stop("Argument 'object' must be of class 'bvarmodel' or 'modellist'.")
  }
  if (is.null(object[["priors"]][["a"]])) {
    stop("Argument 'object' has no prior on its coefficients. Call add_priors() first.")
  }

  model <- object[["model"]]
  k <- model[["k"]]
  p <- model[["p"]]

  # Shrinkage ----
  if (is.character(shrinkage)) {
    shrinkage <- list("type" = shrinkage)
  }
  if (!is.list(shrinkage) || length(shrinkage[["type"]]) != 1 ||
      !shrinkage[["type"]] %in% c("none", "minnesota", "horseshoe", "normal_gamma")) {
    stop("Argument 'shrinkage' must be 'none', 'minnesota', 'horseshoe' or 'normal_gamma', ",
         "or a list whose element 'type' is one of them.")
  }
  object[["priors"]][["a"]][["shrinkage"]] <- NULL
  model[["shrinkage"]] <- NULL
  if (shrinkage[["type"]] != "none") {
    n_par <- length(object[["priors"]][["a"]][["mu"]])
    group <- shrinkage[["group"]]
    if (is.null(group)) {
      # vec(A) of the k x n_x coefficient matrix: position (j - 1) k + i is
      # equation i and regressor j, and the first k p regressors are the lags.
      eq <- rep(seq_len(k), length.out = n_par)
      reg <- rep(seq_len(n_par / k), each = k)
      lag <- reg <= k * p
      own <- lag & ((reg - 1) %% k) + 1 == eq
      if (shrinkage[["type"]] == "minnesota") {
        group <- ifelse(own, 1, ifelse(lag, 2, 0))
        # A model of one variable has no other lags, and groups are numbered
        # without gaps.
        if (!any(group == 2) || !any(group == 1)) {
          group[group > 0] <- 1
        }
      } else if (shrinkage[["type"]] == "normal_gamma") {
        # One global rate per lag order, as in Huber and Feldkircher (2019).
        group <- ifelse(lag, (reg - 1) %/% k + 1, 0)
      } else {
        group <- ifelse(lag, 1, 0)
      }
    }
    if (length(group) != n_par) {
      stop("Element 'group' of argument 'shrinkage' must have one element per coefficient, ",
           n_par, ".")
    }
    prior <- list("group" = matrix(as.numeric(group), ncol = 1))
    if (shrinkage[["type"]] == "minnesota") {
      n_groups <- max(group)
      prior[["shape"]] <- .column(if (is.null(shrinkage[["shape"]])) 3 else shrinkage[["shape"]], n_groups)
      prior[["rate"]] <- .column(if (is.null(shrinkage[["rate"]])) 2 else shrinkage[["rate"]], n_groups)
    }
    if (shrinkage[["type"]] == "normal_gamma") {
      n_groups <- max(group)
      prior[["shape"]] <- .column(if (is.null(shrinkage[["shape"]])) 0.01 else shrinkage[["shape"]], n_groups)
      prior[["rate"]] <- .column(if (is.null(shrinkage[["rate"]])) 0.01 else shrinkage[["rate"]], n_groups)
      prior[["theta"]] <- .column(if (is.null(shrinkage[["theta"]])) 0.1 else shrinkage[["theta"]], n_groups)
    }
    object[["priors"]][["a"]][["shrinkage"]] <- prior
    model[["shrinkage"]] <- shrinkage[["type"]]
  }

  # Stationarity ----
  if (!is.logical(stationary) || length(stationary) != 1 || is.na(stationary)) {
    stop("Argument 'stationary' must be TRUE or FALSE.")
  }
  model[["stationary"]] <- if (stationary) TRUE else NULL

  # Steady state ----
  object[["priors"]][["mu"]] <- NULL
  model[["steady_state"]] <- NULL
  if (!is.null(steady_state)) {
    if (!is.list(steady_state) || is.null(steady_state[["mu"]]) || is.null(steady_state[["v_i"]])) {
      stop("Argument 'steady_state' must be a list with elements 'mu' and 'v_i'.")
    }
    object[["priors"]][["mu"]] <- list("mu" = .column(steady_state[["mu"]], k),
                                       "v_inv" = diag(rep_len(as.numeric(steady_state[["v_i"]]), k), k))
    model[["steady_state"]] <- TRUE
  }

  # Soft constraints ----
  object[["priors"]][["constraints"]] <- NULL
  n_soft <- .soft_groups(object)
  if (!is.null(constraints) && n_soft == 0) {
    stop("Argument 'constraints' is given, and the model has no soft constraints.")
  }
  if (n_soft > 0) {
    if (is.null(constraints)) {
      constraints <- list("shape" = 3, "rate" = 0.01)
    }
    object[["priors"]][["constraints"]] <- list("shape" = .column(constraints[["shape"]], n_soft),
                                                "rate" = .column(constraints[["rate"]], n_soft))
  }

  object[["model"]] <- model
  return(object)
}

# A value recycled to 'n' as a column matrix, which is how the priors keep a
# vector: the HDF5 writer stores it as the (1, n) row BayesTS reads.
.column <- function(value, n) {
  matrix(rep_len(as.numeric(value), n), ncol = 1)
}

# The number of groups of soft constraints in the sample.
.soft_groups <- function(object) {
  group <- object[["data"]][["train"]][["constraints"]][["group"]]
  if (length(group) == 0) 0 else max(group)
}

# The starting values the options of add_prior_options() and the constraints
# of create_bvarmodel() need, where the object does not carry them: filled in
# just before the sampler runs, so that add_initial_values() need not know
# about them and a model given starting values of its own keeps them. The
# unconditional mean starts at its prior mean, the precision of the soft
# constraints at the mean of its prior. The scales of a shrinkage prior are
# left to the sampler, which starts them at one.
.complete_initial_values <- function(object) {

  if (isTRUE(object[["model"]][["steady_state"]]) && is.null(object[["initial"]][["mu"]])) {
    object[["initial"]][["mu"]] <- object[["priors"]][["mu"]][["mu"]]
  }
  if (.soft_groups(object) > 0 && is.null(object[["initial"]][["constraints_inv"]])) {
    prior <- object[["priors"]][["constraints"]]
    if (is.null(prior)) {
      stop("The model has soft constraints and no prior on their error precision. ",
           "Call add_prior_options() after add_priors().")
    }
    object[["initial"]][["constraints_inv"]] <- prior[["shape"]] / prior[["rate"]]
  }

  object
}
