#' @include multipliers.R
NULL

#' Dynamic Multipliers of a VAR Model with Exogenous Variables
#'
#' Computes the response of an endogenous variable to a change in a weakly
#' exogenous variable of an object of class 'bvarmodel'.
#'
#' @param x an object of class 'bvarmodel' with at least one exogenous
#' variable.
#' @param impulse name of the exogenous variable that is moved.
#' @param response name of the endogenous variable whose response is returned.
#' @param n_ahead number of steps ahead. Zero is allowed and returns the impact
#' response alone.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered
#' by the credible intervals. Defaults to 0.95.
#' @param shock size of the change in the exogenous variable, in its own units.
#' Defaults to one.
#' @param type \code{"permanent"} (default) for a change that is held from
#' period zero on, or \code{"transitory"} for a change in period zero alone.
#' See 'Details'.
#' @param cumulative logical specifying whether the responses should be
#' cumulated.
#' @param keep_draws logical specifying whether the function should return all
#' draws of the posterior multipliers. Defaults to \code{FALSE}, so that the
#' median and the credible intervals of the posterior draws are returned.
#' @param period integer. Index of the period the coefficients are taken from.
#' Only used for models with time varying parameters. Defaults to \code{NULL},
#' so that the draws of the last period are used.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' For the model
#' \deqn{y_t = \sum_{l = 1}^{p} A_{l} y_{t-l} + \sum_{j = 0}^{s} B_{j} x_{t-j} +
#' C d_t + u_t,}
#' where \eqn{x_t} is weakly exogenous, the dynamic multipliers \eqn{M_h} are
#' obtained from the recursion
#' \deqn{M_h = \sum_{l = 1}^{\min(h, p)} A_{l} M_{h-l} + \sum_{j = 0}^{\min(h, s)} B_{j} e,}
#' with \eqn{M_h = 0} for \eqn{h < 0} and \eqn{e} the unit vector that selects
#' the exogenous variable in \code{impulse}. The sum over \eqn{B_j} is what
#' makes the change permanent: the variable is one unit higher from period zero
#' on, so every lag of it that has come into range contributes. With
#' \code{type = "transitory"} the variable is one unit higher in period zero
#' alone, the sum collapses to \eqn{B_h}, and the response returns to zero
#' unless the endogenous block has a unit root.
#'
#' The deterministic terms and the errors play no part: a multiplier is a
#' difference between two paths of the same model that differ only in the
#' exogenous variable, so everything both paths share cancels. The responses are
#' in the units of the endogenous variables per unit of the exogenous one, and
#' they are computed draw by draw, so what is returned is a posterior
#' distribution of multipliers rather than a point estimate of them.
#'
#' Whether the multipliers settle down is a property of the endogenous block
#' alone. If the companion matrix of the \eqn{A_l} has a root on or outside the
#' unit circle -- which is the normal case for a model in levels whose variables
#' are integrated -- a permanent change moves the endogenous variables
#' permanently, and the multipliers converge to a level rather than to zero.
#'
#' @return A time-series object of class 'bvarirf', which is what
#' \code{\link{irf}} returns, so that the same \code{plot} method applies.
#'
#' @examples
#'
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' # Investment and income as endogenous, consumption as weakly exogenous
#' model <- create_bvarmodel(data = e1[, c("invest", "income")],
#'                           exogen = e1[, "cons", drop = FALSE],
#'                           p = 2, s = 1, deterministic = "const",
#'                           iterations = 100, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' model <- add_posterior_coefficients(add_initial_values(model))
#'
#' dm <- multipliers(model, impulse = "cons", response = "income", n_ahead = 8)
#' plot(dm)
#'
#' @references
#'
#' Pesaran, M. H., Schuermann, T., Weiner, S. M. (2004). Modeling regional
#' interdependencies using a global error-correcting macroeconometric model.
#' \emph{Journal of Business & Economic Statistics, 22}(2), 129-162.
#'
#' @family post-estimation analysis
#' @export
#' @method multipliers bvarmodel
multipliers.bvarmodel <- function(x, impulse = NULL, response = NULL,
                                  n_ahead = 5, ci = .95, shock = 1,
                                  type = "permanent", cumulative = FALSE,
                                  keep_draws = FALSE, period = NULL, ...) {

  if (!type %in% c("permanent", "transitory")) {
    stop("Argument 'type' must be \"permanent\" or \"transitory\".")
  }
  if (length(n_ahead) != 1 || !is.numeric(n_ahead) || is.na(n_ahead) || n_ahead < 0) {
    stop("Argument 'n_ahead' must be a single integer of at least 0.")
  }
  if (length(shock) != 1 || !is.numeric(shock) || is.na(shock)) {
    stop("Argument 'shock' must be a single number.")
  }

  specs <- x[["model"]]
  k <- specs[["k"]]
  p <- specs[["p"]]
  m <- specs[["m"]]
  s <- specs[["s"]]

  if (is.null(m) || m == 0) {
    stop("Dynamic multipliers need a model with weakly exogenous variables, ",
         "and this model has none.")
  }
  # The coefficients of a structural model are the structural ones, and the
  # recursion below is the reduced form's. Reading the first for the second
  # would give multipliers that belong to no model.
  if (isTRUE(specs[["structural"]])) {
    stop("Dynamic multipliers are not available for a structural model.")
  }
  if (is.null(x[["posterior"]][["a"]][["coeffs"]])) {
    stop("The model holds no posterior draws of its coefficients.")
  }

  exogen <- specs[["exogen"]]
  endogen <- specs[["endogen"]]
  impulse <- which(exogen == impulse)
  if (length(impulse) != 1) {
    stop("Impulse variable not available among the exogenous variables.")
  }
  response <- which(endogen == response)
  if (length(response) != 1) {
    stop("Response variable not available among the endogenous variables.")
  }

  draws <- x[["posterior"]][["a"]][["coeffs"]]
  store <- nrow(draws)

  # With time varying parameters a row holds the whole path, so the block of
  # one period is cut out of it, as the impulse responses do.
  offset <- 0
  if (isTRUE(specs[["tvp"]])) {
    nparams <- ncol(x[["data"]][["train"]][["z"]])
    tt <- nrow(x[["data"]][["train"]][["y"]])
    if (is.null(tt) || tt == 0) {
      tt <- length(x[["data"]][["train"]][["y"]]) / k
    }
    period <- if (is.null(period)) tt else .check_period(period, tt)
    offset <- (period - 1) * nparams
  }

  n_a <- k * k * p
  n_b <- k * m * (s + 1)

  result <- matrix(NA_real_, store, n_ahead + 1)

  for (i in seq_len(store)) {

    a_i <- if (p > 0) {
      matrix(draws[i, offset + seq_len(n_a)], k)
    } else {
      matrix(0, k, 0)
    }
    b_i <- matrix(draws[i, offset + n_a + seq_len(n_b)], k)

    # The column of each B block that belongs to the exogenous variable moved.
    b <- lapply(0:s, function(j) b_i[, j * m + impulse])

    path <- matrix(0, n_ahead + 1, k)
    for (h in 0:n_ahead) {
      value <- rep(0, k)
      if (p > 0) {
        for (l in seq_len(min(h, p))) {
          value <- value + a_i[, (l - 1) * k + 1:k] %*% path[h - l + 1, ]
        }
      }
      if (type == "permanent") {
        for (j in 0:min(h, s)) {
          value <- value + b[[j + 1]]
        }
      } else {
        if (h <= s) {
          value <- value + b[[h + 1]]
        }
      }
      path[h + 1, ] <- value
    }

    result[i, ] <- path[, response] * shock
  }

  if (cumulative) {
    result <- matrix(apply(result, 1, cumsum), nrow = nrow(result), byrow = TRUE)
  }

  if (!keep_draws) {
    result <- .summarise_irf_draws(result, ci)
  }

  class(result) <- append("bvarirf", class(result))
  result
}
