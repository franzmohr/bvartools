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
#' do not vary by period are left as they are. The posterior of a discounted
#' model, which holds the moments of every period rather than draws, is cut to
#' the same periods. In either case the draws or moments that remain were
#' estimated from the whole sample; re-estimate the model to condition on the
#' window alone.
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
                                          .path_widths(x, k), .is_discount(x))
  }

  # A sign restricted identification holds for the period it was found in,
  # which is an index into the sample and so moves with it. The last period is
  # what no stored period means, and it is not the last period of the window.
  sign <- x[["model"]][["sign_restrictions"]]
  if (!is.null(sign)) {
    stored <- if (is.null(sign[["period"]])) length(orig_time) else sign[["period"]]
    moved <- match(stored, periods)
    if (is.na(moved)) {
      x[["model"]][["sign_restrictions"]] <- NULL
      x[["posterior"]][["q"]] <- NULL
      warning("The sign restrictions were imposed in a period the window leaves out, so ",
              "the identification was dropped. Run add_sign_restrictions() again.",
              call. = FALSE)
    } else {
      x[["model"]][["sign_restrictions"]][["period"]] <- as.integer(moved)
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