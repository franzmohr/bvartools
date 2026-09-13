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
#' as under stochastic volatility -- and the pointwise log-likelihood. Draws that
#' do not vary by period are left as they are.
#'
#' @return An object of class 'bvarmodel'.
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

  if (!is.null(x[["posterior"]])) {
    x[["posterior"]] <- .window_posterior(x[["posterior"]], periods, length(orig_time),
                                          .path_widths(x, k))
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
                 "psi$coeffs" = k * (k - 1) / 2,
                 "u_sigma_inv$coeffs" = k * k,
                 "u_omega_inv$coeffs" = k,
                 "u_scale$coeffs" = k,
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
.window_posterior <- function(posterior, periods, tt, widths) {

  if (tt < 2) {
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