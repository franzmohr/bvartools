#' @include add_sign_zero_restrictions.R
NULL

#' Sign and Zero Restrictions
#'
#' Identifies the shocks of an object of class 'bvarmodel' by the signs of the
#' impulse responses they produce together with responses that are restricted
#' to be exactly zero.
#'
#' @param object an object of class 'bvarmodel', containing posterior draws of
#' the coefficients and the error covariance.
#' @param restrictions a data frame of the restrictions, with one row per
#' restriction. See 'Details'.
#' @param draws integer. The number of draws to resample. Defaults to
#' \code{NULL}, which uses the effective sample size of the importance sampler,
#' the number of independent draws it actually produced. Asking for more than
#' that is allowed and duplicates draws.
#' @param one_sided logical. Should the numerical derivative behind the
#' importance weights be taken on one side only? It is about forty percent
#' faster and correspondingly less accurate. Defaults to \code{FALSE}.
#' @param period integer. Index of the period, whose draws should be identified.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior
#' draws of the last time period are used.
#' @param ... further arguments passed to or from other methods.
#'
#' @details A sign restriction can be imposed by trying rotations until one of
#' them carries the signs that were asked for, which is what
#' \code{\link{add_sign_restrictions}} does. A zero restriction cannot: the
#' rotations that satisfy one form a set of probability zero, so no number of
#' tries finds a single member of it. This function uses the algorithms of
#' Arias, Rubio-Ramirez and Waggoner (2018) instead, and needs at least one zero
#' restriction to be worth using -- with only signs, \code{add_sign_restrictions}
#' produces the same draws more cheaply, and this function says so rather than
#' running.
#'
#' The restrictions are imposed on \eqn{F(A_0, A_+)}, the responses stacked over
#' the horizons the restrictions mention. That function satisfies
#' \eqn{F(A_0 Q, A_+ Q) = F(A_0, A_+) Q} for every orthogonal \eqn{Q}, so a zero
#' restriction on the response to the \eqn{j}th shock is a \emph{linear}
#' restriction on the \eqn{j}th column of \eqn{Q} once the reduced form is
#' fixed. The columns are therefore built one at a time, each drawn uniformly
#' from the sphere in the subspace that the zero restrictions and the columns
#' already placed leave over. Every rotation the function draws satisfies the
#' zero restrictions exactly rather than approximately.
#'
#' That construction does not draw from the posterior conditional on the zero
#' restrictions but from a distribution that differs from it by a volume
#' element, which the function computes numerically and divides out. The result
#' is an importance sample: each draw carries a weight, and the draws are
#' resampled with replacement according to those weights so that what the
#' function returns is an ordinary, equally weighted sample. A draw whose
#' rotation fails the sign restrictions is given weight zero and so never
#' resampled; unlike the rejection sampler no second rotation is tried for it,
#' and a shock is not retried with its sign flipped, since the sphere its column
#' is drawn from already covers both of its signs.
#'
#' \strong{The draws that come back are a resample and no longer a chain.}
#' Their order carries no information, several of them may be copies of the same
#' original draw, and convergence diagnostics computed on them mean nothing. How
#' much independent information they carry is the effective sample size, which
#' \code{\link{summary}} reports.
#'
#' Argument \code{restrictions} is a data frame with the columns
#' \describe{
#'  \item{\code{impulse}}{name of the endogenous variable the shock is named
#'  after.}
#'  \item{\code{response}}{name of the endogenous variable whose response is
#'  restricted.}
#'  \item{\code{sign}}{\code{1} for a response that must be positive, \code{-1}
#'  for one that must be negative, and \code{0} for one restricted to be exactly
#'  zero.}
#'  \item{\code{horizon}}{optional horizon the restriction applies to, counted
#'  from zero for the impact period. \code{Inf} restricts the long-run response.
#'  Defaults to \code{0}.}
#' }
#'
#' \strong{The order of the endogenous variables matters here in a way it does
#' not for sign restrictions alone.} Columns of \eqn{Q} are built in the order
#' the variables appear, and the \eqn{j}th column has to keep at least one
#' dimension after the \eqn{j-1} columns before it and its own zero
#' restrictions have been taken out. A shock carrying many zero restrictions
#' therefore has to be named early. The function refuses an ordering that leaves
#' a column with nothing to draw and says which shock it was.
#'
#' The result depends on the state of the random number generator, so
#' \code{\link{set.seed}} is needed to reproduce it.
#'
#' Applied to a 'modellist' or an 'expandingwindow' the function identifies
#' each member on its own and returns the collection. Each therefore gets its
#' own effective sample size, and unless \code{draws} is given the members come
#' back with different numbers of draws -- which is what it means for one model,
#' or one window of the sample, to support the restrictions less well than
#' another. Pass \code{draws} to give them all the same number.
#'
#' @return The object of class 'bvarmodel' with its posterior draws resampled,
#' the accepted rotations in element \code{q} of its \code{posterior}, one row
#' per resampled draw, and the specification of the restrictions in element
#' \code{sign_zero_restrictions} of its \code{model}, together with the number
#' of draws that satisfied the sign restrictions and the effective sample size
#' of the importance sampler. Element \code{sign_restrictions} is set alongside
#' it, so that \code{\link{irf}}, \code{\link{fevd}} and \code{\link{spillover}}
#' read the rotations under \code{type = "sign"} as they do for
#' \code{\link{add_sign_restrictions}}.
#'
#' @examples
#'
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
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
#' # Consumption does not respond to the investment shock on impact, and income
#' # rises. Investment is the first variable of the data set, which is what
#' # leaves its shock a column to draw: see 'Details' on the ordering.
#' restrictions <- data.frame(impulse = "invest",
#'                            response = c("cons", "income"),
#'                            sign = c(0, 1),
#'                            horizon = 0)
#'
#' set.seed(1234)
#' model <- add_sign_zero_restrictions(model, restrictions)
#'
#' # Obtain the identified impulse response
#' ir <- irf(model, impulse = "invest", response = "income", type = "sign")
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
#' @family post-estimation analysis
#' @export
#' @method add_sign_zero_restrictions bvarmodel
add_sign_zero_restrictions.bvarmodel <- function(object, restrictions, draws = NULL,
                                                 one_sided = FALSE, period = NULL, ...) {

  .refuse_unrotatable(object, "Sign and zero restrictions")

  k <- object[["model"]][["k"]]
  restrictions <- .check_sign_zero_restrictions(restrictions, object[["model"]][["endogen"]], k)

  if (!is.logical(one_sided) || length(one_sided) != 1 || is.na(one_sided)) {
    stop("Argument 'one_sided' must be either TRUE or FALSE.")
  }

  A <- .collect_draws(object, period = period, need_A0 = FALSE, need_Sigma = TRUE,
                      all_regressors = TRUE)
  store <- length(A)
  m <- ncol(A[[1]][["A"]])

  setup <- .arw_setup(restrictions, k, m, object[["model"]][["p"]])

  q <- matrix(NA_real_, store, k * k)
  log_weight <- rep(-Inf, store)
  for (i in seq_len(store)) {
    identified <- .arw_draw_q(A[[i]], setup, weight = TRUE, one_sided = one_sided)
    if (length(identified[["q"]]) > 0) {
      q[i, ] <- as.numeric(identified[["q"]])
      log_weight[i] <- identified[["log_weight"]]
    }
  }

  accepted <- sum(!is.na(q[, 1]))
  if (accepted == 0) {
    stop("No rotation satisfying the sign restrictions was found for any of the ", store,
         " posterior draws. Every rotation the zero restrictions admit carries the wrong ",
         "signs, so either the restrictions contradict each other or the model does not ",
         "produce the pattern they describe.", call. = FALSE)
  }

  # A weight that could not be computed is not a draw that was rejected: it is
  # one whose volume element came back singular. Both are dropped, but only the
  # second is worth telling the user about.
  unweighted <- sum(!is.na(q[, 1]) & is.na(log_weight))
  if (unweighted > 0) {
    warning("The importance weight of ", unweighted, " of the ", accepted,
            " draws satisfying the sign restrictions could not be computed and they were ",
            "dropped. Their volume element was singular, which a near-singular draw of the ",
            "error covariance can cause.", call. = FALSE)
    log_weight[is.na(log_weight)] <- -Inf
  }

  weights <- .arw_weights(log_weight)
  effective <- max(1L, as.integer(floor(1 / sum(weights^2))))

  if (is.null(draws)) {
    draws <- effective
  } else {
    if (!is.numeric(draws) || length(draws) != 1 || is.na(draws) ||
        draws < 1 || draws != round(draws)) {
      stop("Argument 'draws' must be a single positive integer.")
    }
    draws <- as.integer(draws)
  }

  index <- sample.int(store, size = draws, replace = TRUE, prob = weights)

  object[["posterior"]] <- .resample_draws(object[["posterior"]], index, store)
  object[["posterior"]][["q"]] <- list(
    "coeffs" = coda::mcmc(q[index, , drop = FALSE], start = 1, end = draws, thin = 1)
  )

  # The resampled draws are one sample rather than several chains, whatever the
  # sampler produced: a resample mixes them, and leaving the old count behind
  # would have thin() split a sample that no longer divides.
  object[["model"]][["chains"]] <- 1L

  object[["model"]][["sign_zero_restrictions"]] <- list(
    "restrictions" = restrictions,
    "one_sided" = one_sided,
    "period" = period,
    "candidates" = store,
    "accepted" = accepted,
    "effective_sample_size" = effective
  )

  # What irf(), fevd() and spillover() read under type = "sign". Both
  # identifications leave the same thing behind -- a rotation per draw and the
  # period it was found in -- so neither has to know about the other.
  object[["model"]][["sign_restrictions"]] <- list(
    "restrictions" = restrictions,
    "period" = period
  )

  return(object)
}


# The restriction table of a sign and zero restricted identification, checked
# and translated into the positions the worker counts in.
#
# It differs from the table add_sign_restrictions() takes in two ways, and both
# are what the algorithm buys: a sign of zero is a restriction rather than an
# error, and a horizon may be infinite, which restricts the long-run response.
.check_sign_zero_restrictions <- function(restrictions, varnames, k) {

  restrictions <- .check_restriction_columns(restrictions, varnames)

  if (!is.numeric(restrictions[["sign"]]) ||
      !all(restrictions[["sign"]] %in% c(-1, 0, 1))) {
    stop("Column 'sign' of argument 'restrictions' must be -1, 0 or 1.")
  }

  horizon <- restrictions[["horizon"]]
  if (!is.numeric(horizon) || anyNA(horizon) || any(horizon < 0) ||
      any(is.finite(horizon) & horizon != round(horizon))) {
    stop("Column 'horizon' of argument 'restrictions' must contain non-negative integers ",
         "or Inf for the long run.")
  }

  if (!any(restrictions[["sign"]] == 0)) {
    stop("Argument 'restrictions' contains no zero restriction, so this function has ",
         "nothing to offer over add_sign_restrictions(), which imposes signs alone by ",
         "trying rotations. Use that one: it draws from the same distribution and does ",
         "not pay for the importance weights.", call. = FALSE)
  }

  cell <- paste(restrictions[["impulse"]], restrictions[["response"]],
                restrictions[["horizon"]], sep = "-")
  signs_per_cell <- tapply(restrictions[["sign"]], cell, function(s) length(unique(s)))
  if (any(signs_per_cell > 1)) {
    stop("Argument 'restrictions' asks for more than one thing of the same response at the ",
         "same horizon, which no rotation can satisfy.")
  }

  # The jth column of the rotation is drawn from what is left of the space once
  # the j - 1 columns before it and its own zero restrictions are taken out.
  # Nothing left means the ordering of the variables, not the restrictions, is
  # what has to change.
  zeros <- table(factor(restrictions[["impulse"]][restrictions[["sign"]] == 0],
                        levels = seq_len(k)))
  for (j in seq_len(k)) {
    if (k - ((j - 1) + zeros[[j]]) < 1) {
      stop("Shock ", j, ", named after '", varnames[j], "', carries ", zeros[[j]],
           " zero restrictions, which with the ", j - 1, " shocks before it leaves nothing ",
           "of its column to draw. Order the endogenous variables so that the shocks with ",
           "the most zero restrictions come first.", call. = FALSE)
    }
  }

  restrictions[, c("impulse", "response", "sign", "horizon")]
}


# The restriction table as the worker wants it: the horizons the responses are
# stacked over, and one block of rows per shock for the zeros and for the signs.
#
# Each block multiplies the stacked responses, so a restriction on the response
# of variable r at the ith horizon is a row with its sign in column
# (i - 1) * k + r.
.arw_setup <- function(restrictions, k, m, lags) {

  horizons <- sort(unique(restrictions[["horizon"]]))
  nh <- length(horizons)

  block <- function(rows) {
    out <- matrix(0, length(rows), nh * k)
    for (i in seq_along(rows)) {
      r <- restrictions[rows[i], ]
      column <- (match(r[["horizon"]], horizons) - 1) * k + r[["response"]]
      out[i, column] <- if (r[["sign"]] == 0) 1 else r[["sign"]]
    }
    out
  }

  z <- vector("list", k)
  s <- vector("list", k)
  w <- vector("list", k)
  for (j in seq_len(k)) {
    z[[j]] <- block(which(restrictions[["impulse"]] == j & restrictions[["sign"]] == 0))
    s[[j]] <- block(which(restrictions[["impulse"]] == j & restrictions[["sign"]] != 0))

    # The fixed matrix that completes the constraints on this column to a square
    # system. Any draw of it defines a valid algorithm -- see Appendix A.3 of
    # the paper -- but it has to be the same one for every posterior draw,
    # because the volume element differentiates the map it takes part in.
    rows <- k - ((j - 1) + nrow(z[[j]]))
    w[[j]] <- matrix(stats::rnorm(rows * k), rows, k)
  }

  list("k" = k, "m" = m, "lags" = lags, "horizons" = horizons,
       "z" = z, "s" = s, "w" = w)
}


# The normalised importance weights. The largest log weight is taken out before
# exponentiating, which is what keeps a draw with a large weight from
# overflowing and every draw from underflowing to zero together.
.arw_weights <- function(log_weight) {
  finite <- is.finite(log_weight)
  weights <- rep(0, length(log_weight))
  weights[finite] <- exp(log_weight[finite] - max(log_weight[finite]))
  weights / sum(weights)
}


# Every element of a posterior that holds one row per draw, resampled by
# `index`. The walk is the one .thin_draws() makes, and for the same reason:
# anything with as many rows as there are draws is a block of draws, and taking
# some blocks and not others would leave row i belonging to different draws in
# different blocks.
#
# The labels start afresh at one. A resample is not a chain and there is no
# iteration number to carry: draws appear more than once and their order says
# nothing.
.resample_draws <- function(posterior, index, draws) {

  for (i in names(posterior)) {
    element <- posterior[[i]]
    if (is.null(element)) {
      next
    }
    if (is.list(element) && !inherits(element, "mcmc")) {
      posterior[[i]] <- .resample_draws(element, index, draws)
    } else if (NROW(element) == draws) {
      posterior[[i]] <- coda::mcmc(.draws_matrix(element)[index, , drop = FALSE],
                                   start = 1, end = length(index), thin = 1)
    }
  }

  return(posterior)
}
