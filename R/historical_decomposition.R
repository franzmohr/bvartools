#' Historical Decomposition
#'
#' Decomposes the path of a variable over the estimation sample into the
#' contributions of the identified structural shocks.
#'
#' @param x an object of class 'bvarmodel' with posterior draws.
#' @param ... further arguments passed to or from other methods.
#'
#' @return See the method, \code{\link{historical_decomposition.bvarmodel}}.
#'
#' @family post-estimation analysis
#' @export
historical_decomposition <- function(x, ...) {
  UseMethod("historical_decomposition")
}

#' Historical Decomposition of a Vector Autoregressive Model
#'
#' Decomposes the path of a variable over the estimation sample into the
#' contributions of the identified structural shocks and a baseline.
#'
#' @param x an object of class 'bvarmodel' with posterior draws of the
#' coefficients and the error covariance.
#' @param response name of the variable whose path is decomposed.
#' @param type the identification: \code{"oir"} (default) for orthogonalised
#' shocks from the Choleski factor of each draw's covariance, \code{"sign"} for
#' the rotations \code{\link{add_sign_restrictions}} or
#' \code{\link{add_sign_zero_restrictions}} stored in the model, and
#' \code{"custom"} for impact matrices given in \code{impact}.
#' @param impact for \code{type = "custom"}, the impact matrix: one matrix, a
#' list of one matrix per draw, or a function of a draw, as in
#' \code{\link{irf.bvarmodel}}. Column \eqn{j} is the impact response to shock
#' \eqn{j}; its column names, where it has any, name the shocks.
#' @param statistic the posterior summary of each contribution,
#' \code{"mean"} (default) or \code{"median"}. Only the means add up to the
#' data exactly.
#' @param ci an optional probability, the coverage of the credible interval
#' returned alongside, such as \code{0.68}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details With \eqn{u_t = P \epsilon_t} the reduced form errors and
#' \eqn{\epsilon_t} the structural shocks, the model's moving average form splits
#' each period's value into
#' \deqn{y_t = b_t + \sum_{j} \sum_{s = 0}^{t - 1} \Phi_s P_{\cdot j} \epsilon_{j, t - s},}
#' the contribution of each shock \eqn{j} over the sample so far plus a
#' baseline \eqn{b_t}: what the lags before the sample, the deterministic terms
#' and the exogenous variables would have produced without any shock. For every
#' draw the shocks are recovered from that draw's residuals and its impact
#' matrix, and the contributions are propagated by the draw's lag coefficients,
#' so that the decomposition carries the posterior uncertainty of all three.
#'
#' Where the panel was not observed whole (see \code{\link{create_bvarmodel}}),
#' each draw is decomposed over its own completed panel, so that the periods
#' nobody observed are decomposed as the model filled them.
#'
#' Available for models with constant coefficients, a lag order of at least one
#' and a covariance that does not change over the sample; that excludes
#' time-varying parameters and stochastic volatility, whose impact matrix would
#' differ from period to period, and structural models.
#'
#' @return A time series of class \code{"bvarhd"}, one row per period of the
#' estimation sample and one column per shock plus a last column
#' \code{baseline}. With \code{ci}, attributes \code{lower} and \code{upper}
#' hold the bounds of the credible interval in the same shape. The response
#' itself is attribute \code{data}; for a completed panel it is the sum of the
#' columns, the posterior mean of the completed series.
#'
#' @examples
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' model <- create_bvarmodel(e1, p = 2, iterations = 100, burnin = 50)
#' model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' model <- add_posterior_coefficients(add_initial_values(model))
#' hd <- historical_decomposition(model, response = "cons")
#' plot(hd)
#'
#' @family post-estimation analysis
#' @export
historical_decomposition.bvarmodel <- function(x, response = NULL, type = "oir", impact = NULL,
                                               statistic = "mean", ci = NULL, ...) {

  if (!type %in% c("oir", "sign", "custom")) {
    stop("Argument 'type' must be \"oir\", \"sign\" or \"custom\".")
  }
  if (type == "custom" && is.null(impact)) {
    stop("A historical decomposition of type \"custom\" needs an impact matrix in argument 'impact'.")
  }
  if (!statistic %in% c("mean", "median")) {
    stop("Argument 'statistic' must be \"mean\" or \"median\".")
  }
  .check_ci(ci)
  .refuse_quantile_covariance(x, "Historical decompositions")

  if (is.null(x[["posterior"]][["a"]][["coeffs"]]) ||
      is.null(x[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Argument 'x' must include posterior draws of the coefficients and the error covariance.")
  }
  if (isTRUE(x[["model"]][["tvp"]])) {
    stop("Historical decompositions are available for models with constant coefficients only.")
  }
  if (isTRUE(x[["model"]][["structural"]])) {
    stop("Historical decompositions are not available for structural models.")
  }
  k <- x[["model"]][["k"]]
  p <- x[["model"]][["p"]]
  tt <- .train_periods(x, k)
  if (p < 1) {
    stop("Historical decompositions need a model with at least one lag.")
  }
  if (.u_sigma_is_path(x, k, tt)) {
    stop("Historical decompositions need an error covariance that is the same in every period; ",
         "the impact matrix of a model with stochastic volatility differs from period to period.")
  }

  varnames <- x[["model"]][["endogen"]]
  if (is.null(response) || !response %in% varnames) {
    stop("Argument 'response' must name an endogenous variable.")
  }
  r <- which(varnames == response)

  if (type == "sign") {
    impact <- .sign_impact(x, "Historical decompositions")
  } else if (type == "oir") {
    impact <- function(draw, i) t(chol(draw[["Sigma"]]))
  }
  shock_names <- varnames
  if (is.matrix(impact) && !is.null(colnames(impact))) {
    shock_names <- colnames(impact)
  }

  A <- .collect_draws(x, need_Sigma = type != "custom", impact = impact, all_regressors = TRUE)
  if (length(A) == 0) {
    stop("No draw has an impact matrix to decompose with.")
  }

  y_data <- as.matrix(x[["data"]][["train"]][["y"]])
  x_data <- as.matrix(x[["data"]][["train"]][["x"]])
  completed <- x[["posterior"]][["y"]][["coeffs"]]
  kept <- if (is.null(completed)) NULL else .kept_draws(x, A)

  result <- array(NA_real_, c(length(A), tt, k + 1))
  for (d in seq_along(A)) {
    y_d <- y_data
    x_d <- x_data
    if (!is.null(completed)) {
      # This draw's completed panel, and the lags rebuilt from it where they
      # fall inside the sample; before it, the regressors keep what was given.
      y_d <- matrix(completed[kept[d], ], tt, k, byrow = TRUE)
      for (l in seq_len(p)) {
        if (tt > l) {
          x_d[(l + 1):tt, (l - 1) * k + 1:k] <- y_d[1:(tt - l), ]
        }
      }
    }
    result[d, , ] <- .hd_draw(y_d, x_d, A[[d]][["A"]], A[[d]][["P"]], k, p, r)
  }

  summarise <- function(f) {
    out <- apply(result, c(2, 3), f)
    out <- stats::ts(out, start = stats::start(x[["data"]][["train"]][["y"]]),
                     frequency = stats::frequency(x[["data"]][["train"]][["y"]]))
    colnames(out) <- c(shock_names, "baseline")
    out
  }
  out <- summarise(if (statistic == "mean") mean else stats::median)
  if (!is.null(ci)) {
    attr(out, "lower") <- summarise(function(v) stats::quantile(v, (1 - ci) / 2, names = FALSE))
    attr(out, "upper") <- summarise(function(v) stats::quantile(v, 1 - (1 - ci) / 2, names = FALSE))
  }
  attr(out, "data") <- stats::ts(y_data[, r], start = stats::start(x[["data"]][["train"]][["y"]]),
                                 frequency = stats::frequency(x[["data"]][["train"]][["y"]]))
  if (!is.null(completed)) {
    attr(out, "data") <- stats::ts(rowSums(out), start = stats::start(x[["data"]][["train"]][["y"]]),
                                   frequency = stats::frequency(x[["data"]][["train"]][["y"]]))
  }
  attr(out, "response") <- response
  class(out) <- append("bvarhd", class(out))
  out
}

# The rows of the posterior draws each element of .collect_draws() came from:
# it drops the draws a sign restricted identification could not identify, and
# a completed panel has to be read at the same rows as the coefficients.
.kept_draws <- function(x, A) {
  store <- nrow(x[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  if (length(A) == store) {
    return(seq_len(store))
  }
  rotations <- x[["posterior"]][["q"]][["coeffs"]]
  if (is.null(rotations)) {
    stop("An impact function that drops draws cannot be combined with a completed panel: ",
         "which draw of the panel belongs to which impact matrix is lost.")
  }
  which(!apply(rotations, 1, anyNA))
}

# One draw's decomposition of response `r`: a tt x (k + 1) matrix, one column
# per shock and the baseline last. `a_full` is the k x n_x coefficient matrix
# whose first k p columns are the lags, `impact` the k x k matrix P with
# u_t = P e_t. The contributions of all shocks are propagated at once: column j
# of M_t is what shock j has done to y_t so far,
#
#     M_t = sum_l A_l M_{t-l} + P diag(e_t),
#
# and the baseline is whatever of the data they leave.
.hd_draw <- function(y, x, a_full, impact, k, p, r) {
  tt <- nrow(y)
  u <- y - x %*% t(a_full)
  e <- u %*% t(solve(impact))
  a_lag <- a_full[, seq_len(k * p), drop = FALSE]
  history <- matrix(0, k * p, k) # M_{t-1} stacked over M_{t-p}
  out <- matrix(0, tt, k + 1)
  for (t in seq_len(tt)) {
    m <- a_lag %*% history + impact %*% diag(e[t, ], k)
    out[t, seq_len(k)] <- m[r, ]
    history <- rbind(m, history[seq_len(k * (p - 1)), , drop = FALSE])
  }
  out[, k + 1] <- y[, r] - rowSums(out[, seq_len(k), drop = FALSE])
  out
}

# A coverage probability, where one is given.
.check_ci <- function(ci) {
  if (!is.null(ci) && (!is.numeric(ci) || length(ci) != 1 || is.na(ci) || ci <= 0 || ci >= 1)) {
    stop("Argument 'ci' must be a single probability between 0 and 1.")
  }
}

#' Plot a Historical Decomposition
#'
#' Stacked bars of the contributions of the shocks in each period, and the
#' part of the response they account for as a line.
#'
#' @param x an object of class \code{"bvarhd"}, the result of
#' \code{\link{historical_decomposition}}.
#' @param baseline logical: should the baseline be stacked with the shocks?
#' Defaults to \code{FALSE}, which plots the response net of the baseline --
#' the part the shocks explain.
#' @param ... further arguments passed to \code{\link[graphics]{barplot}}.
#'
#' @return \code{x}, invisibly.
#'
#' @export
plot.bvarhd <- function(x, baseline = FALSE, ...) {
  orig_par <- graphics::par(mar = graphics::par("mar"))
  on.exit(graphics::par(orig_par))

  parts <- unclass(x)
  attributes(parts) <- list(dim = dim(x), dimnames = dimnames(x))
  if (!baseline) {
    parts <- parts[, colnames(parts) != "baseline", drop = FALSE]
  }
  line <- rowSums(parts)
  positive <- pmax(parts, 0)
  negative <- pmin(parts, 0)
  colours <- grDevices::gray.colors(ncol(parts))

  dots <- list(...)
  args <- list(ylab = attr(x, "response"), border = NA, space = 0,
               ylim = range(c(rowSums(positive), rowSums(negative), line)))
  args <- args[!(names(args) %in% names(dots))]
  mids <- do.call(graphics::barplot, c(list(t(positive), col = colours), args, dots))
  graphics::barplot(t(negative), col = colours, border = NA, space = 0, add = TRUE,
                    axes = FALSE)
  graphics::lines(mids, line, lwd = 2)

  at <- pretty(seq_along(line))
  at <- at[at >= 1 & at <= length(line)]
  graphics::axis(1, at = mids[at], labels = round(stats::time(x)[at], 2))
  graphics::legend("topleft", legend = colnames(parts), fill = colours, bty = "n", cex = .8)

  invisible(x)
}
