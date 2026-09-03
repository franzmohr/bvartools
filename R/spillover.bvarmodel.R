#' Spillover Index
#'
#' Produces the connectedness measures of Diebold and Yilmaz (2012) for an
#' object of class 'bvarmodel'.
#'
#' @param object an object of class 'bvarmodel'.
#' @param n_ahead number of steps ahead. Defaults to 10, the horizon of Diebold
#' and Yilmaz (2012).
#' @param type type of the impulse responses the decomposition is based on.
#' Possible choices are generalised \code{gir} (default) and orthogonalised
#' \code{oir}. See 'Details'.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param keep_draws logical specifying whether the function should return all draws of
#' the indices. Defaults to \code{FALSE}, so that the credible intervals are returned.
#' @param period integer. Index of the period, for which the measures should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of
#' the last time period are used.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The function produces the connectedness measures of Diebold and
#' Yilmaz (2012) for the VAR model
#' \deqn{y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}.
#'
#' Let \eqn{\theta_{jk}(H)} be the share of the \eqn{H} step forecast error
#' variance of variable \eqn{j} that is attributed to a shock to variable
#' \eqn{k}. Under \code{type = "gir"} this is the generalised decomposition of
#' Pesaran and Shin (1998),
#' \deqn{\theta_{jk}(H) = \frac{\sigma^{-1}_{kk} \sum_{i = 0}^{H - 1} (e_j^{\prime} \Phi_i \Sigma e_k)^2}{\sum_{i = 0}^{H - 1} (e_j^{\prime} \Phi_i \Sigma \Phi_i^{\prime} e_j)},}
#' where \eqn{\Phi_i} is the forecast error impulse response of period \eqn{i}
#' and \eqn{\sigma_{kk}} the variance of the \emph{shock} variable. These shares
#' do not add up over \eqn{k}, so each row is normalised,
#' \deqn{\tilde{\theta}_{jk}(H) = \theta_{jk}(H) / \sum_{k} \theta_{jk}(H).}
#' Under \code{type = "oir"} the decomposition uses the Choleski factor of
#' \eqn{\Sigma}, adds up by construction and depends on the ordering of the
#' variables, which is what the generalised version of Diebold and Yilmaz (2012)
#' avoids.
#'
#' From the normalised table the measures are
#' \describe{
#'  \item{\code{total}}{the total spillover index,
#'  \eqn{100 \sum_{j \neq k} \tilde{\theta}_{jk} / k}, the share of forecast
#'  error variance across the system that comes from other variables.}
#'  \item{\code{from}}{spillovers received by each variable from all others,
#'  \eqn{100 (1 - \tilde{\theta}_{jj}) / k}.}
#'  \item{\code{to}}{spillovers transmitted by each variable to all others,
#'  \eqn{100 \sum_{j \neq k} \tilde{\theta}_{jk} / k} summed down column \eqn{k}.}
#'  \item{\code{net}}{\code{to} less \code{from}. Positive values mark net
#'  transmitters of shocks and negative values net receivers.}
#'  \item{\code{pairwise}}{net pairwise spillovers,
#'  \eqn{100 (\tilde{\theta}_{kj} - \tilde{\theta}_{jk}) / k}.}
#' }
#'
#' Every measure is computed once per posterior draw and only then summarised.
#' That is not interchangeable with computing it from the posterior mean table:
#' the measures are ratios, so the index of the mean is not the mean of the
#' index. It is also what gives the index a credible interval, which the
#' original point estimate based version does not have.
#'
#' Objects of class \code{'bvecmodel'} have to be transformed with
#' \code{\link{vec_to_var}} first, as for \code{\link{fevd}}.
#'
#' @return An object of class 'bvarspillover', a list containing
#' \item{total}{the total spillover index. The median and the bounds of the
#' credible interval, or all draws if \code{keep_draws = TRUE}.}
#' \item{to, from, net}{the directional measures, one column per variable, in
#' the same form.}
#' \item{table}{the posterior mean of the normalised decomposition table, in
#' percent, with the \code{from} column, the \code{to} row and the total index
#' in the corner attached, so it is \eqn{(k + 1) \times (k + 1)}.}
#' \item{pairwise}{the posterior mean of the net pairwise spillovers.}
#' \item{specification}{a list recording \code{n_ahead}, \code{type}, \code{ci},
#' \code{period}, \code{k} and the number of draws.}
#'
#' @examples
#'
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' e1 <- window(e1, end = c(1978, 4))
#'
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#'
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#'
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' object <- add_posterior_coefficients(model)
#'
#' # Obtain the connectedness measures
#' sp <- spillover(object)
#' sp
#'
#' # Net directional spillovers
#' plot(sp)
#'
#' @references
#'
#' Diebold, F. X., & Yilmaz, K. (2012). Better to give than to receive: Predictive
#' directional measurement of volatility spillovers. \emph{International Journal of
#' Forecasting, 28}(1), 57--66.
#'
#' Diebold, F. X., & Yilmaz, K. (2014). On the network topology of variance
#' decompositions: Measuring the connectedness of financial firms.
#' \emph{Journal of Econometrics, 182}(1), 119--134.
#'
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in
#' linear multivariate models. \emph{Economics Letters, 58}, 17--29.
#'
#' @export
spillover.bvarmodel <- function(object, n_ahead = 10, type = "gir", ci = .95,
                                keep_draws = FALSE, period = NULL, ...) {

  if (is.null(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Argument 'object' must include draws of the variance-covariance matrix Sigma.")
  }

  if (!type %in% c("gir", "oir")) {
    stop("Argument 'type' must be either 'gir' or 'oir'.")
  }

  if (length(n_ahead) != 1 || n_ahead < 1) {
    stop("Argument 'n_ahead' must be a single integer of at least 1.")
  }

  if (length(ci) != 1 || ci <= 0 || ci >= 1) {
    stop("Argument 'ci' must be a single value between 0 and 1.")
  }

  if (object[["model"]][["p"]] == 0) {
    stop("Spillover measures are only supported for models with p > 0.")
  }

  k <- object[["model"]][["k"]]
  if (k < 2) {
    stop("Spillover measures need at least two endogenous variables.")
  }
  varnames <- object[["model"]][["endogen"]]

  # Shared with fevd, so the two agree on which slice of a row is `period`.
  A <- .collect_draws(object, period = period, need_A0 = FALSE)
  store <- length(A)

  tables <- lapply(A, .spillover_table, h = n_ahead, type = type)

  # One set of measures per draw. The table is normalised, so it sums to k and
  # that is the denominator of every one of them.
  total <- numeric(store)
  from <- matrix(NA_real_, store, k)
  to <- matrix(NA_real_, store, k)
  mean_table <- matrix(0, k, k)
  mean_pairwise <- matrix(0, k, k)

  for (i in seq_len(store)) {
    theta <- tables[[i]]
    own <- diag(theta)

    from[i, ] <- 100 * (1 - own) / k
    to[i, ] <- 100 * (colSums(theta) - own) / k
    # Computed off the table rather than as sum(from), so that the accounting
    # identity sum(from) == sum(to) == total stays something a test can check.
    total[i] <- 100 * (sum(theta) - sum(own)) / k

    mean_table <- mean_table + theta
    mean_pairwise <- mean_pairwise + 100 * (t(theta) - theta) / k
  }

  net <- to - from

  mean_table <- mean_table / store
  mean_pairwise <- mean_pairwise / store
  dimnames(mean_table) <- list(varnames, varnames)
  dimnames(mean_pairwise) <- list(varnames, varnames)

  # The table as this literature prints it: the shares, a column of what each
  # variable receives, a row of what it transmits, and the total index in the
  # corner. In percent throughout -- the body is scaled up to meet the margins,
  # which are percentages already, so that no part of this object is in
  # different units from another.
  printed <- cbind(100 * mean_table, "from" = colMeans(from))
  printed <- rbind(printed, "to" = c(colMeans(to), mean(total)))

  result <- list(
    total = .spillover_summary(matrix(total, ncol = 1), ci, keep_draws, "total"),
    from = .spillover_summary(from, ci, keep_draws, varnames),
    to = .spillover_summary(to, ci, keep_draws, varnames),
    net = .spillover_summary(net, ci, keep_draws, varnames),
    table = printed,
    pairwise = mean_pairwise,
    specification = list(n_ahead = n_ahead, type = type, ci = ci,
                         period = period, k = k, draws = store,
                         endogen = varnames)
  )

  class(result) <- append("bvarspillover", class(result))

  return(result)
}

# Draws of one measure summarised the way irf.bvarmodel summarises its own: the
# median between the two bounds of the credible interval, unless the draws
# themselves were asked for.
.spillover_summary <- function(draws, ci, keep_draws, varnames) {
  colnames(draws) <- varnames
  if (keep_draws) {
    return(coda::as.mcmc(draws))
  }
  ci_low <- (1 - ci) / 2
  pr <- c(ci_low, .5, 1 - ci_low)
  # apply() drops to a vector for a single column, so the shape is restored:
  # every element of the result carries the same three rows.
  out <- matrix(apply(draws, 2, stats::quantile, probs = pr), nrow = length(pr))
  dimnames(out) <- list(paste0(pr * 100, "%"), varnames)
  out
}
