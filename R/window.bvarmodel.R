#' Time Series Windows
#' 
#' Restricts the observations in the training sample of a model of class
#' 'bvarmodel' to the specified time window.
#'
#' @param x an object of class 'bvarmodel'.
#' @param start the start time of the period of interest.
#' @param end the end time of the period of interest.
#' @param ... further arguments passed to or from other methods.
#'
#' @details Posterior draws that form a path over the periods of the training
#' sample are cut to the periods that remain, so that they still refer to the
#' same observations as the data. These are the coefficients of a model with
#' time varying parameters, the draws of the error term that vary by period --
#' as under stochastic volatility -- the panel completed by a model with
#' \code{missing = "estimate"}, and the pointwise log-likelihood. Draws that
#' do not vary by period are left as they are. The posterior of a discounted
#' model, which holds the moments of every period rather than draws, is cut to
#' the same periods. In either case the draws or moments that remain were
#' estimated from the whole sample; re-estimate the model to condition on the
#' window alone.
#'
#' \strong{A window that ends before the sample drops the forecast input, the
#' forecasts and the realised values they were scored against}, with a warning:
#' they start after the original end, from its last observations. Call
#' \code{add_forecast_input()} again to forecast from the end of the window.
#'
#' @return An object of class 'bvarmodel'.
#'
#' @examples
#'
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' model <- add_priors(model,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     sigma = list(df = "k", scale = 1))
#' model <- add_posterior_coefficients(add_initial_values(model))
#'
#' # The last ten years of the sample
#' recent <- window(model, start = c(1972, 1))
#' start(recent$data$train$y)
#'
#' @export
#' @method window bvarmodel
window.bvarmodel <- function(x, start = NULL, end = NULL, ...) {

  k <- x[["model"]][["k"]]
  orig_time <- stats::time(x[["data"]][["train"]][["y"]])

  x[["data"]][["train"]][["y"]] <- stats::window(x[["data"]][["train"]][["y"]],
                                                         start = start, end = end, ...)

  if (!is.null(x[["data"]][["train"]][["x"]])) {
    x[["data"]][["train"]][["x"]] <- stats::window(x[["data"]][["train"]][["x"]],
                                                        start = start, end = end, ...)
  }

  periods <- which(orig_time %in% stats::time(x[["data"]][["train"]][["y"]]))

  if (!is.null(x[["data"]][["train"]][["z"]])) {
    pos <- rep(periods * k, each = k) - k + 1 + rep(0:(k - 1), length(periods))
    x[["data"]][["train"]][["z"]] <- x[["data"]][["train"]][["z"]][pos,]
  }

  # The observations of a panel not observed whole that the window leaves out,
  # and the starting values of what it did not observe, follow the cut.
  x <- .window_panel(x, periods, orig_time)

  .warn_empty_dummy_variables(x)

  if (!is.null(x[["posterior"]])) {
    x[["posterior"]] <- .window_posterior(x[["posterior"]], periods, length(orig_time),
                                          .path_widths(x, k), .is_discount(x))
  }
  x <- .window_forecast(x, periods, length(orig_time))

  # A sign restricted identification holds for the period it was found in,
  # which is an index into the sample and so moves with it. The last period is
  # what no stored period means, and it is not the last period of the window.
  sign <- x[["model"]][["sign_restrictions"]]
  if (!is.null(sign)) {
    stored <- if (is.null(sign[["period"]])) length(orig_time) else sign[["period"]]
    moved <- match(stored, periods)
    if (is.na(moved)) {
      x[["model"]][["sign_restrictions"]] <- NULL
      x[["model"]][["sign_zero_restrictions"]] <- NULL
      x[["posterior"]][["q"]] <- NULL
      warning("The sign restrictions were imposed in a period the window leaves out, so ",
              "the identification was dropped. Run add_sign_restrictions() again.",
              call. = FALSE)
    } else {
      x[["model"]][["sign_restrictions"]][["period"]] <- as.integer(moved)
      # The record of add_sign_zero_restrictions() names the same period, and
      # was left pointing into the original sample.
      if (!is.null(x[["model"]][["sign_zero_restrictions"]])) {
        x[["model"]][["sign_zero_restrictions"]][["period"]] <- as.integer(moved)
      }
    }
  }

  return(x)
}


# How many columns each posterior block holds per period if it is a path.
#
# A block that does not vary by period holds exactly this many columns, and a
# path this many times the number of periods, so the width alone tells the two
# apart; no block needs to be known in advance to be time varying. Only these
# blocks are considered, because something else -- a forecast, say -- can have
# a number of columns that is a multiple of the periods without being a path.
.path_widths <- function(x, k) {
  widths <- list("a$coeffs" = NCOL(x[["data"]][["train"]][["z"]]),
                 # The samplers store the whole lower triangular Psi per
                 # period, not its k(k - 1)/2 free elements, which is only
                 # what initial$psi holds.
                 "psi$coeffs" = k * k,
                 "u_sigma_inv$coeffs" = k * k,
                 "u_omega_inv$coeffs" = k,
                 "u_scale$coeffs" = k,
                 # The panel completed by a model with missing = "estimate",
                 # one draw per row, k values per period.
                 "y$coeffs" = k,
                 # The mixing scale of Student t errors, one per period.
                 "student_t$lambda" = 1,
                 "loglik" = 1)
  if (isTRUE(x[["model"]][["rank"]] > 0) && !is.null(x[["data"]][["train"]][["w"]])) {
    widths[["beta$coeffs"]] <- NCOL(x[["data"]][["train"]][["w"]]) * x[["model"]][["rank"]]
  }
  return(widths)
}


# Cuts the posterior paths among 'widths' to the retained 'periods' of a sample
# of 'tt' periods, keeping the attributes of the mcmc objects. Used to be
# missing entirely: window() cut the data and left every path at the length of
# the original sample, so that a summary of the last period of the window
# silently reported a period of the original sample instead.
.window_posterior <- function(posterior, periods, tt, widths, discount = FALSE) {

  if (tt < 2) {
    return(posterior)
  }

  # A discounted model has no draws but its posterior itself, one row per
  # period: the smoothed moments of every period of the sample. They are cut to
  # the periods that remain, as a path of draws is below.
  if (discount) {
    for (path in list(c("a", "mean"), c("a", "scale"), c("a", "cov"),
                      c("u_sigma", "scale"), "df")) {
      block <- tryCatch(posterior[[path]], error = function(e) NULL)
      if (!is.null(block) && NROW(block) == tt) {
        posterior[[path]] <- block[periods, , drop = FALSE]
      }
    }
    # The log predictive density of each period, one column per period. Left
    # whole, the marginal likelihood summed from it covered the original sample.
    loglik <- posterior[["loglik"]]
    if (!is.null(loglik) && NCOL(loglik) == tt) {
      posterior[["loglik"]] <- loglik[, periods, drop = FALSE]
    }
    return(posterior)
  }

  for (name in names(widths)) {
    path <- strsplit(name, "$", fixed = TRUE)[[1]]
    draws <- tryCatch(posterior[[path]], error = function(e) NULL)
    width <- widths[[name]]
    if (is.null(draws) || width < 1 || NCOL(draws) != width * tt) {
      next
    }
    columns <- rep((periods - 1) * width, each = width) + seq_len(width)
    mcpar <- attr(draws, "mcpar")
    cut <- .draws_matrix(draws)[, columns, drop = FALSE]
    if (!is.null(mcpar)) {
      cut <- coda::mcmc(cut, start = mcpar[1], end = mcpar[2], thin = mcpar[3])
    }
    posterior[[path]] <- cut
  }

  return(posterior)
}

# Drops what belongs to the end of the original sample when the window moves
# it: the forecast input, whose lags are the last observations, the horizon,
# the realised values, and the forecasts and their scores. Kept, they described
# forecasts made from periods the window no longer ends in -- predict() dated
# them after the window and they started from the old last observations.
.window_forecast <- function(x, periods, tt) {
  if (length(periods) == 0 || max(periods) == tt) {
    return(x)
  }
  had <- !is.null(x[["data"]][["forecast"]]) || !is.null(x[["data"]][["test"]]) ||
    !is.null(x[["model"]][["h"]]) || !is.null(x[["posterior"]][["forecast"]])
  if (!had) {
    return(x)
  }
  x[["data"]][["forecast"]] <- NULL
  x[["data"]][["test"]] <- NULL
  x[["model"]][["h"]] <- NULL
  x[["model"]][["forecast_start"]] <- NULL
  x[["posterior"]][["forecast"]] <- NULL
  warning("The window ends before the original sample, so the forecast input, the ",
          "forecasts and their test values, which start after the original end, were ",
          "dropped. Use add_forecast_input() again.", call. = FALSE)
  x
}
