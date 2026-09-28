#' Sign and Zero Restricted Rotations of Posterior Draws
#'
#' The identification step of \code{\link{add_sign_zero_restrictions}} on its
#' own: rotations of the reduced-form draws of a VAR that satisfy sign and zero
#' restrictions on the impulse responses, with the importance weights of Arias,
#' Rubio-Ramirez and Waggoner (2018). It is exported for packages that estimate
#' a VAR as part of something larger -- the state equation of a factor
#' augmented VAR, for example -- and want to identify it the same way.
#'
#' @param draws a list with one element per posterior draw, each a list with
#' elements \code{A}, the \eqn{K \times M} coefficient matrix with the
#' \eqn{p} lag blocks first and any deterministic terms after them, and
#' \code{Sigma}, the \eqn{K \times K} error covariance.
#' @param restrictions a data frame of the restrictions, as in
#' \code{\link{add_sign_zero_restrictions}}: columns \code{impulse},
#' \code{response}, \code{sign} and optionally \code{horizon}, the first two
#' naming elements of \code{variables}. A table without a zero restriction is
#' accepted here; the rotations are then drawn uniformly and the weights are
#' equal.
#' @param variables the names of the \eqn{K} variables, in the order of the
#' rows of \code{A}. Shocks are named after them and their columns are drawn in
#' this order.
#' @param lags the lag order \eqn{p}.
#' @param max_tries integer. The largest number of rotations drawn for each
#' posterior draw. See 'Details' of \code{\link{add_sign_zero_restrictions}}.
#' @param smooth logical. Should the importance weights be Pareto smoothed?
#' @param one_sided logical. Should the numerical derivative behind the
#' importance weights be taken on one side only?
#'
#' @details The rotation \eqn{Q} of a draw maps into the impact matrix
#' \eqn{L Q}, with \eqn{L} the lower Cholesky factor of \code{Sigma}: the
#' response of the variables to the shocks on impact.
#'
#' The function draws from the random number generator in the order
#' \code{\link{add_sign_zero_restrictions}} does, so the same seed gives the
#' same rotations through either. It warns when the importance sample is not
#' fit to summarise, and stops when no draw has an admissible rotation.
#'
#' @return A list with
#' \item{q}{a matrix with one row per draw holding its rotation, by column, and
#' \code{NA} for a draw with no admissible rotation;}
#' \item{weights}{the normalised importance weights, zero for such a draw;}
#' \item{restrictions}{the restrictions, with positions in place of names;}
#' \item{tries}{the number of rotations drawn in all;}
#' \item{accepted}{the number of draws with an admissible rotation;}
#' \item{effective_sample_size, max_weight_share, pareto_k}{the diagnostics of
#' the importance sample, as \code{\link{add_sign_zero_restrictions}} records
#' them.}
#'
#' @references
#'
#' Arias, J. E., Rubio-Ramirez, J. F., Waggoner, D. F. (2018). Inference based on structural vector
#' autoregressions identified with sign and zero restrictions: Theory and applications.
#' \emph{Econometrica, 86}(2), 685-720. \doi{10.3982/ECTA14468}
#'
#' @seealso \code{\link{add_sign_zero_restrictions}}, which applies it to a
#' 'bvarmodel' and resamples the draws.
#'
#' @export
arias_rubio_ramirez_waggoner_2018 <- function(draws, restrictions, variables, lags,
                                              max_tries = 1, smooth = TRUE,
                                              one_sided = FALSE) {

  if (!is.list(draws) || length(draws) == 0) {
    stop("Argument 'draws' must be a non-empty list of draws.")
  }
  k <- length(variables)
  for (i in c(1, length(draws))) {
    d <- draws[[i]]
    if (!is.matrix(d[["A"]]) || !is.matrix(d[["Sigma"]]) || nrow(d[["A"]]) != k ||
        !identical(dim(d[["Sigma"]]), c(k, k))) {
      stop("Each element of 'draws' must hold a ", k, " x M matrix 'A' and a ", k,
           " x ", k, " matrix 'Sigma', one row and column per element of 'variables'.")
    }
  }
  if (!is.numeric(lags) || length(lags) != 1 || is.na(lags) || lags < 0 ||
      lags != round(lags) || ncol(draws[[1]][["A"]]) < k * lags) {
    stop("Argument 'lags' must be a non-negative integer, and 'A' must hold that many ",
         "lag blocks.")
  }
  if (!is.logical(one_sided) || length(one_sided) != 1 || is.na(one_sided)) {
    stop("Argument 'one_sided' must be either TRUE or FALSE.")
  }
  if (!is.logical(smooth) || length(smooth) != 1 || is.na(smooth)) {
    stop("Argument 'smooth' must be either TRUE or FALSE.")
  }
  if (!is.numeric(max_tries) || length(max_tries) != 1 || is.na(max_tries) ||
      max_tries < 1 || max_tries != round(max_tries) || max_tries > .Machine$integer.max) {
    stop("Argument 'max_tries' must be a single positive integer.")
  }
  max_tries <- as.integer(max_tries)

  restrictions <- .check_sign_zero_restrictions(restrictions, variables, k,
                                                require_zero = FALSE)

  store <- length(draws)
  m <- ncol(draws[[1]][["A"]])
  setup <- .arw_setup(restrictions, k, m, lags)

  q <- matrix(NA_real_, store, k * k)
  log_weight <- rep(-Inf, store)
  tries <- 0
  for (i in seq_len(store)) {
    identified <- .arw_draw_q(draws[[i]], setup, weight = TRUE, one_sided = one_sided,
                              max_tries = max_tries)
    tries <- tries + identified[["tries"]]
    if (length(identified[["q"]]) > 0) {
      q[i, ] <- as.numeric(identified[["q"]])
      log_weight[i] <- identified[["log_weight"]]
    }
  }

  accepted <- sum(!is.na(q[, 1]))
  if (accepted == 0) {
    stop("No rotation satisfying the sign restrictions was found for any of the ", store,
         " posterior draws, with ", max_tries, if (max_tries == 1) " try" else " tries",
         " each. ",
         if (max_tries == 1) {
           paste0("With many sign restrictions a single try rarely satisfies them all; ",
                  "'max_tries' draws more rotations per draw. ")
         },
         "Otherwise every rotation the zero restrictions admit carries the wrong ",
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

  smoothed <- .arw_weights(log_weight, smooth = smooth)
  weights <- smoothed[["weights"]]
  pareto_k <- smoothed[["pareto_k"]]
  effective <- max(1L, as.integer(floor(1 / sum(weights^2))))
  largest <- max(weights)

  # An importance sampler can fail quietly. When one draw carries most of the
  # weight the effective sample size collapses, and because a caller's resample
  # defaults to it the posterior handed back would be a handful of rows that
  # irf() and fevd() would summarise without complaint. Both halves of that are
  # worth saying out loud, and they are separate symptoms: a small effective
  # sample can come from many mildly unequal weights, while a tail too heavy for
  # the estimator to have a finite variance is a different problem with a
  # different fix.
  .warn_importance_sample(effective, largest, accepted, pareto_k)

  list("q" = q,
       "weights" = weights,
       "restrictions" = restrictions,
       "tries" = tries,
       "accepted" = accepted,
       "effective_sample_size" = effective,
       "max_weight_share" = largest,
       "pareto_k" = pareto_k)
}
