#' Test for Time Variation
#'
#' Generic function used to compute Bayes factors for time variation in the
#' coefficients and volatilities of a model with time varying parameters.
#'
#' @param object an object with suitable posterior draws passed forward to method.
#' @param ... arguments passed forward to method.
#'
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{time_variation_test.bvarmodel}}.
#'
#' @export
time_variation_test <- function(object, ...) {
  UseMethod("time_variation_test")
}


#' @export
time_variation_test.default <- function(object, ...) {
  stop("time_variation_test() is available for 'bvarmodel' objects with time varying ",
       "parameters and stochastic volatility, estimated with the prior 'omega_v' in ",
       "add_priors().", call. = FALSE)
}


#' Test for Time Variation in a VAR Model
#'
#' Computes the Bayes factors of Chan (2018) for time variation in each
#' coefficient, each covariance coefficient and each log-volatility of a VAR
#' model with time varying parameters and stochastic volatility, from the draws
#' of a single estimation of that model.
#'
#' @param object an object of class \code{"bvarmodel"} with \code{tvp = TRUE} and
#' \code{error = "sv"} or \code{"sv+covar"}, whose priors were set with
#' \code{coef$omega_v} or \code{sigma$omega_v} in \code{\link{add_priors}} and
#' whose posterior was drawn by \code{\link{add_posterior_coefficients}}.
#' @param joint logical. Should the Bayes factor for the time variation of every
#' state of a block at once be reported as well? Default is \code{TRUE}.
#' @param batches integer. The number of batches the draws are split into for
#' the numerical standard errors. Default is 20.
#' @param ... further arguments, which are ignored.
#'
#' @details Under the prior \code{omega_v} a random walk
#' \eqn{x_t = x_{t-1} + v_t}, \eqn{v_t \sim N(0, \omega^2)}, is estimated in the
#' non-centred form \eqn{x_t = x_0 + \omega \tilde{x}_t} of Frühwirth-Schnatter and
#' Wagner (2010), with \eqn{\tilde{x}_t} a standard random walk and the prior
#' \eqn{\omega \sim N(0, V_\omega)}. The state does not move exactly when
#' \eqn{\omega = 0}, which is a point inside the prior, so the Bayes factor of the
#' model in which the state varies against the one in which it is constant is the
#' Savage-Dickey density ratio
#' \deqn{BF = \frac{p(\omega = 0)}{p(\omega = 0 | y)}.}
#' The numerator is the density of the prior at zero. The denominator is
#' estimated by the average over the draws of the density at zero of the
#' conditional posterior of \eqn{\omega}, which the sampler stores for every draw
#' in \code{omega_log_zero} (Chan 2018). A positive log Bayes factor favours
#' time variation. On the scale of Kass and Raftery (1995), values of the log
#' Bayes factor between 1 and 3 are positive evidence, between 3 and 5 strong
#' evidence and above 5 very strong evidence, and the same values with a
#' negative sign are evidence for a constant state.
#'
#' The joint Bayes factor of a block compares the model in which every state of
#' the block varies with the model in which none does. It is not a test of
#' whether at least one state varies: every state that is constant costs the
#' joint model about the same as one of its own Bayes factors against time
#' variation, so a block in which one state moves and many do not can come out
#' against time variation as a whole. Conversely, the joint Bayes factor can be
#' far larger than the individual ones: each state may have posterior mass near
#' zero on its own while little of the joint posterior sits where all of them
#' are near zero at once. Chan (2018, Table 2 and footnote 8) reports log Bayes
#' factors of about 3 for each of two volatilities of Italian inflation and of
#' 235 for the two together.
#'
#' The estimate is least precise where the Bayes factor is large, since the
#' ordinates are then an average of small numbers (Chan 2018, section 2.3). The
#' numerical standard error of each log Bayes factor is obtained by the delta
#' method from batch means, and says how much of the value is owed to the length
#' of the chain. Where it is large, more draws are needed before the value is
#' read as more than its sign.
#'
#' Each block chooses its prior separately in \code{\link{add_priors}}, and only
#' the blocks estimated under \code{omega_v} are reported.
#'
#' @return A data frame of class \code{"bvartimevar"} with one row per state and,
#' if \code{joint = TRUE}, one per block, and the columns
#' \describe{
#'   \item{\code{block}}{\code{"coefficients"}, \code{"covariances"} or
#'   \code{"volatilities"}.}
#'   \item{\code{equation}}{the endogenous variable of the equation the state
#'   belongs to; \code{NA} for the joint rows.}
#'   \item{\code{term}}{the regressor of a coefficient, the variable whose error
#'   a covariance coefficient loads on, \code{"log-volatility"}, or
#'   \code{"(joint)"}.}
#'   \item{\code{log_bf}}{the log Bayes factor in favour of time variation.}
#'   \item{\code{nse}}{its numerical standard error.}
#' }
#'
#' @references
#'
#' Chan, J. C. C. (2018). Specification tests for time-varying parameter models with stochastic
#' volatility. \emph{Econometric Reviews, 37}(8), 807--823. \doi{10.1080/07474938.2016.1167948}
#'
#' Frühwirth-Schnatter, S., & Wagner, H. (2010). Stochastic model specification search for
#' Gaussian and partial non-Gaussian state space models. \emph{Journal of Econometrics, 154}(1),
#' 85--100. \doi{10.1016/j.jeconom.2009.07.003}
#'
#' Kass, R. E., & Raftery, A. E. (1995). Bayes factors. \emph{Journal of the American
#' Statistical Association, 90}(430), 773--795. \doi{10.1080/01621459.1995.10476572}
#'
#' @examples
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' # Generate model
#' model <- create_bvarmodel(e1, p = 1, deterministic = "const", tvp = TRUE,
#'                           error = "sv", iterations = 100, burnin = 50)
#'
#' # The non-centred prior on the coefficients and the log-volatilities
#' model <- add_priors(model,
#'                     coef = list(v_i = 1 / 10, omega_v = 0.001),
#'                     sigma = list(mu = 0, v_i = 1 / 100, omega_v = 0.1,
#'                                  state_variance = 0.05, offset = 1e-4))
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' model <- add_posterior_coefficients(model)
#'
#' # Bayes factors for time variation
#' time_variation_test(model)
#'
#' @export
#' @rdname time_variation_test.bvarmodel
time_variation_test.bvarmodel <- function(object, joint = TRUE, batches = 20, ...) {

  if (!is.logical(joint) || length(joint) != 1 || is.na(joint)) {
    stop("Argument 'joint' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(batches) || length(batches) != 1 || is.na(batches) ||
      batches < 2 || batches != round(batches)) {
    stop("Argument 'batches' must be a single integer of at least 2.", call. = FALSE)
  }

  posterior <- object[["posterior"]]
  if (is.null(posterior)) {
    stop("The model has no posterior draws. Use add_posterior_coefficients() first.",
         call. = FALSE)
  }

  k <- object[["model"]][["k"]]
  y_names <- object[["model"]][["endogen"]]
  if (is.null(y_names)) {
    y_names <- paste0("y", seq_len(k))
  }

  # The posterior block, the prior group its omega_v is in, and what the rows
  # of the block are called.
  blocks <- list(
    list(posterior = "a", prior = "a", title = "coefficients",
         labels = .time_variation_labels_a(object, y_names)),
    list(posterior = "psi", prior = "psi", title = "covariances",
         labels = .time_variation_labels_psi(k, y_names)),
    list(posterior = "u_sigma_inv", prior = "u_sigma", title = "volatilities",
         labels = list(equation = y_names, term = rep("log-volatility", k))))

  result <- NULL
  for (b in blocks) {
    draws <- posterior[[b[["posterior"]]]]
    omega_v <- object[["priors"]][[b[["prior"]]]][["omega_v"]]
    if (is.null(draws[["omega_log_zero"]]) || is.null(omega_v)) {
      next
    }

    log_zero <- .draws_matrix(draws[["omega_log_zero"]])
    omega_v <- rep(as.numeric(omega_v), length.out = ncol(log_zero))
    prior_zero <- stats::dnorm(0, 0, sqrt(omega_v), log = TRUE)

    equation <- b[["labels"]][["equation"]]
    term <- b[["labels"]][["term"]]
    if (length(equation) != ncol(log_zero)) {
      equation <- rep(NA_character_, ncol(log_zero))
      term <- paste0(b[["posterior"]], "[", seq_len(ncol(log_zero)), "]")
    }

    rows <- data.frame(
      block = b[["title"]],
      equation = equation,
      term = term,
      log_bf = prior_zero - apply(log_zero, 2, .log_mean_exp),
      nse = apply(log_zero, 2, .nse_log_mean_exp, batches = batches),
      stringsAsFactors = FALSE)

    if (joint && !is.null(draws[["omega_log_zero_joint"]])) {
      log_zero_joint <- as.numeric(.draws_matrix(draws[["omega_log_zero_joint"]]))
      rows <- rbind(rows, data.frame(
        block = b[["title"]],
        equation = NA_character_,
        term = "(joint)",
        log_bf = sum(prior_zero) - .log_mean_exp(log_zero_joint),
        nse = .nse_log_mean_exp(log_zero_joint, batches = batches),
        stringsAsFactors = FALSE))
    }

    result <- rbind(result, rows)
  }

  if (is.null(result)) {
    stop("No block of the model was estimated under the non-centred prior. Set ",
         "'coef$omega_v' or 'sigma$omega_v' in add_priors() for a model with ",
         "tvp = TRUE and error = \"sv\" or \"sv+covar\", and draw the posterior again.",
         call. = FALSE)
  }

  rownames(result) <- NULL
  class(result) <- c("bvartimevar", "data.frame")
  return(result)
}


# The equation and the regressor of every coefficient, in the order of the
# draws: vec(A) with A of dimension k x n_x, then for a structural model the free
# elements of A_0 column by column, as summary.bvarmodel() reads them.
.time_variation_labels_a <- function(object, y_names) {

  k <- object[["model"]][["k"]]
  n_a <- ncol(object[["data"]][["train"]][["z"]])
  if (is.null(n_a)) {
    return(list(equation = character(0), term = character(0)))
  }

  n_structural <- if (isTRUE(object[["model"]][["structural"]])) k * (k - 1) / 2 else 0
  n_non_structural <- n_a - n_structural
  x_names <- .get_regressor_names_bvarmodel(object, add_block = FALSE)
  n_x <- n_non_structural / k
  if (is.null(x_names) || length(x_names) < n_x) {
    x_names <- paste0("x", seq_len(n_x))
  }

  pos <- seq_len(n_non_structural) - 1
  equation <- y_names[pos %% k + 1]
  term <- x_names[pos %/% k + 1]

  if (n_structural > 0) {
    struct <- which(lower.tri(matrix(0, k, k)), arr.ind = TRUE)
    equation <- c(equation, y_names[struct[, 1]])
    term <- c(term, paste0(y_names[struct[, 2]], " (contemporaneous)"))
  }

  list(equation = equation, term = term)
}


# The row and the column of every free element of Psi, row by row -- the order
# of the covariance block in the vendored core, (2,1), (3,1), (3,2), ... -- which
# is not the column by column order of lower.tri().
.time_variation_labels_psi <- function(k, y_names) {
  if (k < 2) {
    return(list(equation = character(0), term = character(0)))
  }
  rows <- unlist(lapply(2:k, function(i) rep(i, i - 1)))
  cols <- unlist(lapply(2:k, function(i) seq_len(i - 1)))
  list(equation = y_names[rows], term = y_names[cols])
}


#' @export
print.bvartimevar <- function(x, digits = 2, ...) {

  cat("Bayes factors for time variation (Savage-Dickey, Chan 2018)\n")
  cat("log BF > 0 favours time variation; |log BF| > 3 is strong evidence\n\n")

  for (b in unique(x[["block"]])) {
    rows <- x[x[["block"]] == b, , drop = FALSE]
    out <- data.frame(equation = ifelse(is.na(rows[["equation"]]), "", rows[["equation"]]),
                      term = rows[["term"]],
                      `log BF` = formatC(rows[["log_bf"]], format = "f", digits = digits),
                      NSE = formatC(rows[["nse"]], format = "f", digits = digits),
                      check.names = FALSE, stringsAsFactors = FALSE)
    cat(toupper(substring(b, 1, 1)), substring(b, 2), ":\n", sep = "")
    print(out, row.names = FALSE, right = FALSE)
    cat("\n")
  }

  invisible(x)
}
