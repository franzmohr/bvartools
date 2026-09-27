#' @include persistence_profiles.R
NULL

#' Persistence Profiles of a VEC Model
#'
#' Calculates the persistence profiles of the cointegrating relations of an
#' object of class 'bvecmodel'.
#'
#' @param object an object of class 'bvecmodel', usually, the result of a call
#' to \code{\link{add_posterior_coefficients}}.
#' @param n_ahead the number of horizons the profile is followed over.
#' @param ci a numeric between 0 and 1 specifying the probability of the
#' credible band. Defaults to 0.95.
#' @param keep_draws logical. If \code{FALSE} (default) the draws are summarised
#' by their median and the band; if \code{TRUE} they are returned as they are.
#' @param period the period of a model with time varying coefficients the
#' profile is calculated at. Defaults to \code{NULL}, the last period. Ignored
#' for a model whose coefficients are constant.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The persistence profile of Pesaran and Shin (1996) is the response of a
#' cointegrating relation to a system-wide shock, scaled to one on impact:
#' \deqn{PP_j(h) = \frac{\beta_j' \Psi_h \Sigma \Psi_h' \beta_j}{\beta_j' \Sigma \beta_j},}
#' where \eqn{\Psi_h} are the moving average coefficients of the level VAR the
#' model implies, \eqn{\Psi_0 = I}, and \eqn{\Sigma} is the covariance of its
#' errors. It starts at one by construction.
#'
#' \strong{A profile that does not fall to zero says the relation is not
#' cointegrating, whatever the rank says.} That is what the statistic is for: a
#' rank chosen by a test or a criterion is an assertion about how many
#' stationary combinations exist, and the profile is the check on it. How
#' quickly the profile falls is the speed of convergence to equilibrium.
#'
#' Unlike \code{\link{irf}} and \code{\link{fevd}}, which a VEC model reaches
#' through \code{\link{vec_to_var}}, this statistic cannot be taken from the
#' level VAR alone: \code{vec_to_var} drops \code{beta}, and the profile is a
#' statement about the cointegrating vectors themselves. The level form is used
#' for the moving average coefficients and the model's own \code{beta} for the
#' relations.
#'
#' \strong{For a model with weakly exogenous variables the profile is partial.}
#' A cointegrating vector of a VECX model spans the domestic variables and the
#' weakly exogenous ones, and only the domestic block responds to a shock to
#' this model. The relation is therefore followed over the part of it this model
#' governs, which is the whole relation only where the model has no
#' \code{exogen}. Where the weakly exogenous variables are endogenous to a
#' larger system -- a global model assembled from sub-models -- take the profile
#' there instead, and the function warns to that effect.
#'
#' @return A list with one element per cointegrating relation, each a matrix
#' with one row per horizon and, unless \code{keep_draws} is \code{TRUE},
#' the columns \code{median}, \code{lower} and \code{upper}.
#'
#' @examples
#'
#' # Load data
#' data("e6")
#'
#' # Create model
#' model <- create_bvecmodel(e6, p = 2, r = 1, const = "unrestricted",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#'
#' model <- add_priors(model,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 0.0001))
#'
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#'
#' profiles <- persistence_profiles(model, n_ahead = 12)
#' round(profiles[[1]][1:5, ], 3)
#'
#' @references
#'
#' Pesaran, M. H., & Shin, Y. (1996). Cointegration and speed of convergence to
#' equilibrium. \emph{Journal of Econometrics, 71}(1-2), 117--143.
#' \doi{10.1016/0304-4076(94)01697-6}
#'
#' @seealso \code{\link{persistence_profiles}} for the generic.
#'
#' @export
#' @method persistence_profiles bvecmodel
persistence_profiles.bvecmodel <- function(object, n_ahead = 20, ci = 0.95,
                                           keep_draws = FALSE, period = NULL, ...) {

  if (ci <= 0 || ci >= 1) {
    stop("Argument 'ci' is not within the permitted range of 0 and 1.")
  }
  if (n_ahead < 1) {
    stop("Argument 'n_ahead' must be at least 1.")
  }

  rank <- object[["model"]][["rank"]]
  if (is.null(rank) || rank < 1) {
    stop("A model of rank ", if (is.null(rank)) "NULL" else rank,
         " has no cointegrating relations, so it has no persistence profiles. ",
         "Estimate the model at a rank of at least one, with argument 'r' of ",
         "create_bvecmodel().", call. = FALSE)
  }
  if (is.null(object[["posterior"]][["beta"]][["coeffs"]])) {
    stop("The model carries no posterior draws of 'beta'. Run ",
         "add_posterior_coefficients() first.", call. = FALSE)
  }

  k <- object[["model"]][["k"]]
  k_w <- ncol(object[["data"]][["train"]][["w"]])
  if (k_w > k) {
    warning("The cointegrating vectors of this model span ", k_w, " variables ",
            "while only ", k, " of them are endogenous to it, so the profiles ",
            "follow the part of each relation this model governs. For the whole ",
            "relation, take the profiles of the system in which the remaining ",
            "variables are endogenous.", call. = FALSE)
  }

  # The level VAR of the model supplies the moving average coefficients, and the
  # model itself the cointegrating vectors, which vec_to_var() drops.
  level <- vec_to_var(object, period = period)
  a <- .draws_matrix(level[["posterior"]][["a"]][["coeffs"]])
  sigma <- .draws_matrix(level[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  beta <- .draws_matrix(object[["posterior"]][["beta"]][["coeffs"]])

  # A time varying model has one beta per period, and the level VAR was taken at
  # one period, so the same one is taken here.
  n_beta <- k_w * rank
  if (ncol(beta) > n_beta) {
    tt <- ncol(beta) / n_beta
    at <- if (is.null(period)) tt else period
    if (at < 1 || at > tt) {
      stop("Argument 'period' must be between 1 and ", tt, ".", call. = FALSE)
    }
    beta <- beta[, (at - 1) * n_beta + seq_len(n_beta), drop = FALSE]
  }

  p <- level[["model"]][["p"]]
  draws <- nrow(a)
  profiles <- array(NA_real_, c(draws, n_ahead + 1, rank))

  for (d in seq_len(draws)) {

    coefficients <- matrix(a[d, ], k)
    covariance <- solve(matrix(sigma[d, ], k))
    vectors <- matrix(beta[d, ], k_w, rank)[seq_len(k), , drop = FALSE]

    # Psi_0 = I and Psi_h = sum_j A_j Psi_{h-j}. Written out rather than taken
    # from .ir(), which returns one impulse and response pair while the profile
    # needs the whole matrix at every horizon.
    psi <- vector("list", n_ahead + 1)
    psi[[1]] <- diag(k)
    for (h in seq_len(n_ahead)) {
      total <- matrix(0, k, k)
      for (j in seq_len(min(h, p))) {
        total <- total + coefficients[, (j - 1) * k + seq_len(k), drop = FALSE] %*% psi[[h - j + 1]]
      }
      psi[[h + 1]] <- total
    }

    for (j in seq_len(rank)) {
      b <- vectors[, j, drop = FALSE]
      scale <- as.numeric(crossprod(b, covariance %*% b))
      for (h in seq_len(n_ahead + 1)) {
        m <- crossprod(b, psi[[h]])
        profiles[d, h, j] <- as.numeric(m %*% covariance %*% t(m)) / scale
      }
    }
  }

  names_of <- colnames(object[["data"]][["train"]][["w"]])
  result <- lapply(seq_len(rank), function(j) {
    x <- profiles[, , j, drop = FALSE][, , 1, drop = TRUE]
    if (is.null(dim(x))) {
      x <- matrix(x, nrow = 1)
    }
    if (keep_draws) {
      colnames(x) <- 0:n_ahead
      return(x)
    }
    low <- (1 - ci) / 2
    out <- cbind(median = apply(x, 2, stats::median),
                 lower = apply(x, 2, stats::quantile, probs = low),
                 upper = apply(x, 2, stats::quantile, probs = 1 - low))
    rownames(out) <- 0:n_ahead
    out
  })
  names(result) <- paste0("relation_", seq_len(rank))
  attr(result, "variables") <- names_of
  class(result) <- append("bvecpp", class(result))
  result
}
