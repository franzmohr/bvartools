#' @include add_sign_restrictions.R
NULL

#' Sign Restrictions
#'
#' Identifies the shocks of an object of class 'bvarmodel' by the signs of the
#' impulse responses they produce.
#'
#' @param object an object of class 'bvarmodel', containing posterior draws of
#' the coefficients and the error covariance.
#' @param restrictions a data frame of the restrictions, with one row per
#' restriction. See 'Details'.
#' @param max_tries the number of rotations that are tried per posterior draw
#' before the draw is given up on. Defaults to 1000.
#' @param period integer. Index of the period, whose draws should be identified.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior
#' draws of the last time period are used.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The function identifies the structural shocks of the VAR model
#' \deqn{y_t = \sum_{i=1}^{p} A_i y_{t - i} + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}, by searching for rotations of its Choleski
#' factor whose impulse responses carry the signs that argument
#' \code{restrictions} asks for.
#'
#' Write \eqn{P} for the lower triangular Choleski factor of \eqn{\Sigma} and
#' \eqn{Q} for an orthogonal matrix. Then \eqn{P Q (P Q)^{\prime} = \Sigma} for
#' every such \eqn{Q}, so a model rotated by \eqn{Q} fits the data exactly as
#' well as the one the posterior draw describes. The likelihood is therefore
#' silent about which rotation is the right one, and the restrictions choose
#' among the models it cannot tell apart. The identification is in consequence
#' \emph{set valued}: what is obtained is not one impulse response per posterior
#' draw but the collection of those that the restrictions admit.
#'
#' For every posterior draw the function draws rotations uniformly over the
#' orthogonal group and keeps the first one whose responses satisfy every
#' restriction. Since a sign restriction cannot distinguish a shock from its own
#' negative, a column that fails is retried with its sign flipped before the
#' rotation is discarded. A draw for which no admissible rotation is found
#' within \code{max_tries} attempts is dropped from the identified sample, and
#' \code{\link{summary}} reports how many were. A low acceptance rate is a
#' statement about the restrictions, not a technicality: it says the model
#' rarely produces the pattern that was asked of it.
#'
#' Argument \code{restrictions} is a data frame with the columns
#' \describe{
#'  \item{\code{impulse}}{name of the endogenous variable the shock is named
#'  after.}
#'  \item{\code{response}}{name of the endogenous variable whose response is
#'  restricted.}
#'  \item{\code{sign}}{either \code{1} for a response that must be positive or
#'  \code{-1} for one that must be negative.}
#'  \item{\code{horizon}}{optional integer of the horizon the restriction
#'  applies to, counted from zero for the impact period. Defaults to \code{0}.}
#' }
#' Restrictions on several shocks may be combined in one data frame. A shock
#' that no row mentions is left as it is drawn, since nothing in the
#' restrictions distinguishes one rotation of it from another; only the
#' restricted shocks should be interpreted.
#'
#' Only sign restrictions are supported. Zero restrictions on the impact
#' responses cannot be imposed by rejection, because the set of rotations that
#' satisfies them has probability zero, and need the algorithm of Arias et al.
#' (2018) instead.
#'
#' The accepted rotations are added to the object as element \code{q} of its
#' posterior draws, from where \code{\link{irf}}, \code{\link{fevd}} and
#' \code{\link{spillover}} use them under \code{type = "sign"}.
#'
#' The result depends on the state of the random number generator, so
#' \code{\link{set.seed}} is needed to reproduce it.
#'
#' @return The object of class 'bvarmodel' with the accepted rotations in
#' element \code{q} of its \code{posterior}, one row per posterior draw, and the
#' specification of the restrictions in element \code{sign_restrictions} of its
#' \code{model}. Draws for which no admissible rotation was found are recorded
#' as \code{NA}.
#'
#' @examples
#'
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#'
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#'
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' model <- add_posterior_coefficients(model)
#'
#' # A demand shock raises investment and consumption on impact
#' restrictions <- data.frame(impulse = "invest",
#'                            response = c("invest", "cons"),
#'                            sign = c(1, 1),
#'                            horizon = 0)
#'
#' set.seed(1234)
#' model <- add_sign_restrictions(model, restrictions)
#'
#' # Obtain the identified impulse response
#' ir <- irf(model, impulse = "invest", response = "cons", type = "sign")
#'
#' @references
#'
#' Arias, J. E., Rubio-Ramirez, J. F., Waggoner, D. F. (2018). Inference based on structural vector
#' autoregressions identified with sign and zero restrictions: Theory and applications.
#' \emph{Econometrica, 86}(2), 685-720. \doi{10.3982/ECTA14468}
#'
#' Rubio-Ramirez, J. F., Waggoner, D. F., Zha, T. (2010). Structural vector autoregressions: Theory of
#' identification and algorithms for inference. \emph{The Review of Economic Studies, 77}(2), 665-696.
#' \doi{10.1111/j.1467-937X.2009.00578.x}
#'
#' Uhlig, H. (2005). What are the effects of monetary policy on output? Results from an agnostic
#' identification procedure. \emph{Journal of Monetary Economics, 52}(2), 381-419.
#' \doi{10.1016/j.jmoneco.2004.05.007}
#'
#' @export
#' @method add_sign_restrictions bvarmodel
add_sign_restrictions.bvarmodel <- function(object, restrictions, max_tries = 1000,
                                            period = NULL, ...) {

  if (is.null(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Argument 'object' must include draws of the variance-covariance matrix Sigma.")
  }

  # The rotation acts on the Choleski factor of the reduced form covariance. A
  # structural model has already spent its identification on A_0 and stores the
  # covariance of the structural errors in its place, so there is no reduced
  # form here to rotate. See irf.bvarmodel for the same restriction.
  if (object[["model"]][["structural"]]) {
    stop("Sign restrictions are not defined for a structural model: they would rotate the ",
         "covariance of the structural errors instead of the reduced form. Estimate the model ",
         "with 'structural = FALSE' to identify it by sign restrictions.")
  }

  # Rotating a covariance that was never estimated is rotating a diagonal
  # matrix, which produces responses that are not those of any shock the data
  # speak about.
  if (object[["model"]][["error"]] %in% c("gamma", "sv", "ald")) {
    stop("Sign restrictions need a model whose error covariances are estimated. Argument ",
         "'object' was estimated with error = \"", object[["model"]][["error"]], "\", which ",
         "leaves the off-diagonal elements at zero, so its Choleski factor is diagonal and a ",
         "rotation of it carries no information about the correlation of the errors.")
  }

  if (object[["model"]][["p"]] == 0) {
    stop("Sign restrictions are only supported for models with p > 0.")
  }

  if (length(max_tries) != 1 || !is.numeric(max_tries) || max_tries < 1) {
    stop("Argument 'max_tries' must be a single integer of at least 1.")
  }
  max_tries <- as.integer(max_tries)

  k <- object[["model"]][["k"]]
  restrictions <- .check_sign_restrictions(restrictions, object[["model"]][["endogen"]])

  # The draws in the shape the worker wants them, shared with irf() and fevd()
  # so that the three cannot disagree about which slice of a row is `period`.
  A <- .collect_draws(object, period = period, need_A0 = FALSE,
                      need_Sigma = TRUE)
  store <- length(A)

  positions <- as.matrix(restrictions[, c("impulse", "response", "sign", "horizon")])

  q <- matrix(NA_real_, store, k * k)
  for (i in seq_len(store)) {
    rotation <- .draw_sign_restricted_q(A[[i]], positions, max_tries)
    if (length(rotation) > 0) {
      q[i, ] <- as.numeric(rotation)
    }
  }

  accepted <- sum(!is.na(q[, 1]))
  if (accepted == 0) {
    stop("No rotation satisfying the restrictions was found for any of the ", store,
         " posterior draws within ", max_tries, " tries each. Either the restrictions ",
         "contradict each other, or the model does not produce the pattern they describe.")
  }

  # Carried with the start, end and thinning interval of the draws it belongs
  # to, so that thin() and window() move it along with them.
  mcpar <- coda::mcpar(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  object[["posterior"]][["q"]] <- list(
    "coeffs" = coda::mcmc(q, start = mcpar[1], end = mcpar[2], thin = mcpar[3])
  )

  # The number of draws that were identified is not recorded here: it is a
  # property of `q`, and thin() would leave a stored count behind.
  object[["model"]][["sign_restrictions"]] <- list(
    "restrictions" = restrictions,
    "max_tries" = max_tries,
    "period" = period
  )

  return(object)
}

# The restriction table, checked and translated into the positions the worker
# counts in.
#
# Variable names are the currency of the rest of the package and positions are
# what the recursion indexes, so the translation has to happen somewhere. Doing
# it once here means the worker never sees a name, and an error message never
# shows a number that the caller did not write.
.check_sign_restrictions <- function(restrictions, varnames) {

  if (!is.data.frame(restrictions)) {
    stop("Argument 'restrictions' must be a data frame.")
  }

  if (nrow(restrictions) == 0) {
    stop("Argument 'restrictions' does not contain any restriction.")
  }

  required <- c("impulse", "response", "sign")
  absent <- required[!required %in% names(restrictions)]
  if (length(absent) > 0) {
    stop("Argument 'restrictions' must contain the column",
         if (length(absent) > 1) "s " else " ",
         paste0("'", absent, "'", collapse = ", "), ".")
  }

  if (is.null(restrictions[["horizon"]])) {
    restrictions[["horizon"]] <- 0L
  }

  for (i in c("impulse", "response")) {
    name <- as.character(restrictions[[i]])
    unknown <- unique(name[!name %in% varnames])
    if (length(unknown) > 0) {
      stop("Column '", i, "' of argument 'restrictions' names ",
           if (length(unknown) > 1) "variables " else "a variable ",
           paste0("'", unknown, "'", collapse = ", "),
           ", which the model does not contain.")
    }
    restrictions[[i]] <- match(name, varnames)
  }

  if (!is.numeric(restrictions[["sign"]]) ||
      !all(restrictions[["sign"]] %in% c(-1, 1))) {
    stop("Column 'sign' of argument 'restrictions' must be either 1 or -1.")
  }

  horizon <- restrictions[["horizon"]]
  if (!is.numeric(horizon) || anyNA(horizon) || any(horizon < 0) ||
      any(horizon != round(horizon))) {
    stop("Column 'horizon' of argument 'restrictions' must contain non-negative integers.")
  }
  restrictions[["horizon"]] <- as.integer(horizon)

  # Two rows that ask opposite things of the same response would reject every
  # rotation, and would do it only after spending the whole budget of tries on
  # every draw. Saying so here costs nothing and names the cause.
  cell <- paste(restrictions[["impulse"]], restrictions[["response"]],
                restrictions[["horizon"]], sep = "-")
  signs_per_cell <- tapply(restrictions[["sign"]], cell,
                           function(s) length(unique(s)))
  if (any(signs_per_cell > 1)) {
    stop("Argument 'restrictions' asks for both signs of the same response at the same ",
         "horizon, which no rotation can satisfy.")
  }

  restrictions[, c("impulse", "response", "sign", "horizon")]
}
