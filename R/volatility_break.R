#' Check a volatility break and place it in the sample
#'
#' Validates the volatility break of Lenza and Primiceri (2022) and returns the
#' two fields a model carries for it, or \code{NULL} when no break was asked
#' for.
#'
#' @param volatility_break the first period of the break as a date of the data,
#' \code{c(year, period)} such as \code{c(2020, 3)}, or a single time such as
#' \code{2020}. \code{NULL} for no break.
#' @param volatility_scales a whole number of at least one, how many periods
#' from the break on get a scale of their own before the excess decays.
#' @param y the endogenous series of the estimation sample, whose
#' \code{\link[stats]{time}} the date is looked up in.
#' @param algorithm the sampler the model would be estimated with.
#'
#' @details
#' \code{\link{create_bvarmodel}} and \code{\link{create_bvecmodel}} call
#' this, and it is exported so that a package building its own models -- the
#' sub-models of a global model in \pkg{bgvars} among them -- validates a
#' break the same way rather than repeating the arithmetic. The break is given
#' as a date and kept as the period of \code{y} it falls on, counted from one,
#' which is what \code{/model/volatility_break} of a model file holds and what
#' the sampler reads.
#'
#' \strong{It is refused for every algorithm but the two constant coefficient
#' Wishart samplers}, \code{VarNormalWishart} and \code{VecNormalWishart},
#' which are the ones that scale the error covariance from a known period on. A
#' break outside the estimation sample is refused as well: the periods after it
#' are what estimate its size.
#'
#' @return A named list with the integer elements \code{volatility_break} and
#' \code{volatility_scales}, or \code{NULL}.
#'
#' @references
#' Lenza, M. and Primiceri, G. E. (2022). How to estimate a vector
#' autoregression after March 2020. \emph{Journal of Applied Econometrics},
#' 37(4), 688--699.
#'
#' @seealso \code{\link{volatility_break_prior}}
#'
#' @export
check_volatility_break <- function(volatility_break, volatility_scales, y, algorithm) {
  if (is.null(volatility_break)) {
    return(NULL)
  }
  if (!algorithm %in% c("VarNormalWishart", "VecNormalWishart")) {
    stop("A volatility break is available for a model with constant coefficients and ",
         "error = \"wishart\", whose sampler scales the error covariance from a known period ",
         "on; this model would be estimated with ", algorithm, ".")
  }
  if (!is.numeric(volatility_break) || !length(volatility_break) %in% 1:2 ||
      anyNA(volatility_break)) {
    stop("Argument 'volatility_break' must be the first period of the break as a date of ",
         "the data, c(year, period) such as c(2020, 3) for March 2020, or a single time ",
         "such as 2020.",
         call. = FALSE)
  }
  if (!is.numeric(volatility_scales) || length(volatility_scales) != 1 ||
      is.na(volatility_scales) || volatility_scales != round(volatility_scales) ||
      volatility_scales < 1) {
    stop("Argument 'volatility_scales' must be a whole number of at least one: the periods ",
         "from the break on that get a scale of their own before it decays.", call. = FALSE)
  }

  frequency <- stats::frequency(y)
  when <- volatility_break[1]
  if (length(volatility_break) == 2) {
    when <- when + (volatility_break[2] - 1) / frequency
  }
  period <- which(abs(stats::time(y) - when) < getOption("ts.eps"))
  if (length(period) != 1) {
    span <- stats::tsp(y)[1:2]
    stop("Argument 'volatility_break' is not a period of the estimation sample, which runs ",
         "from ", format(span[1]), " to ", format(span[2]), " once the lags are taken: ",
         "the break has to fall inside it, since the periods after it are what estimate ",
         "its size.", call. = FALSE)
  }

  list("volatility_break" = as.integer(period),
       "volatility_scales" = as.integer(volatility_scales))
}

#' The prior of a volatility break
#'
#' Checks \code{sigma$volatility} against the model and returns the prior a
#' volatility break is drawn under, or \code{NULL} for a model without one.
#'
#' @param object a model object whose \code{model$volatility_break} says
#' whether it has a break.
#' @param sigma the list given to \code{\link{add_priors}}, whose element
#' \code{volatility} is the prior: \code{shape}, the Pareto shape of every
#' free scale, and \code{rho_shape1} and \code{rho_shape2}, the two shapes of
#' the beta prior on the rate the excess decays at.
#'
#' @details
#' Exported for the same reason as \code{\link{check_volatility_break}}: a
#' package that builds its own models checks the prior the same way. There is
#' no default, as there is none for the rest of \code{sigma}, and the message
#' names the values of Lenza and Primiceri (2022). A prior given to a model
#' without a break is refused rather than ignored, and a break without a prior
#' is refused as well.
#'
#' @return A named list of the prior, or \code{NULL}.
#'
#' @seealso \code{\link{check_volatility_break}}
#'
#' @export
volatility_break_prior <- function(object, sigma) {
  prior <- sigma[["volatility"]]
  if (is.null(object[["model"]][["volatility_break"]])) {
    if (!is.null(prior)) {
      stop("Element 'volatility' of argument 'sigma' is the prior of a volatility break, and ",
           "this model has none. Create it with argument 'volatility_break'.")
    }
    return(NULL)
  }
  if (is.null(prior)) {
    stop("This model has a volatility break, whose prior is element 'volatility' of ",
         "argument 'sigma': list(shape = 1, rho_shape1 = 3.04, rho_shape2 = 1.51) is the ",
         "one of Lenza and Primiceri (2022), a Pareto with shape one on every scale and a ",
         "beta with mode 0.8 and standard deviation 0.2 on the rate the break decays at.")
  }
  names_needed <- c("shape", "rho_shape1", "rho_shape2")
  if (!is.list(prior) || !all(names_needed %in% names(prior))) {
    stop("Argument 'sigma$volatility' must be a list with elements 'shape', 'rho_shape1' ",
         "and 'rho_shape2'.")
  }
  for (i in names_needed) {
    value <- prior[[i]]
    if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || value <= 0) {
      stop("Argument 'sigma$volatility$", i, "' must be a single positive number.")
    }
  }
  list("shape" = prior[["shape"]],
       "rho_shape1" = prior[["rho_shape1"]],
       "rho_shape2" = prior[["rho_shape2"]])
}

# Names the columns of posterior$volatility$coeffs: the free scales s1, s2, ...
# from the period of the break on, then rho, the rate the last one's excess
# over one decays at. NULL passes through, as it does for a model without one.
.label_volatility_draws <- function(volatility) {
  if (is.null(volatility[["coeffs"]])) {
    return(volatility)
  }
  m <- ncol(volatility[["coeffs"]]) - 1
  colnames(volatility[["coeffs"]]) <- c(paste0("s", seq_len(m)), "rho")
  volatility
}

# The period of a volatility break moves with a window of the sample, which
# window() keeps as 'periods', indices into the original sample. A window that
# ends before the break leaves a sample without one, and the break is dropped
# with its prior and draws, with a warning; one that starts after it would need
# the break before the first period, which the sampler cannot represent.
.window_volatility_break <- function(x, periods) {
  period <- x[["model"]][["volatility_break"]]
  if (is.null(period)) {
    return(x)
  }
  moved <- period - periods[1] + 1
  if (moved < 1) {
    stop("The window starts after the volatility break, and a break has to fall inside ",
         "the sample of the model. Create the model again on the shorter sample, without ",
         "'volatility_break' if the window is past the break's effect.", call. = FALSE)
  }
  if (moved > length(periods)) {
    x[["model"]][["volatility_break"]] <- NULL
    x[["model"]][["volatility_scales"]] <- NULL
    if (!is.null(x[["priors"]])) {
      x[["priors"]][["volatility"]] <- NULL
    }
    if (!is.null(x[["posterior"]])) {
      x[["posterior"]][["volatility"]] <- NULL
    }
    warning("The window ends before the volatility break, so the break was dropped ",
            "together with its prior and its draws.", call. = FALSE)
    return(x)
  }
  x[["model"]][["volatility_break"]] <- as.integer(moved)
  x
}

# The free parameters a volatility break adds, for the counts AIC, BIC and HQ
# penalise by: the scales of the periods that have one of their own, and the
# rate the last one's excess over one decays at.
#
# Counted against the sample rather than from 'volatility_scales' alone, because
# a scale the sample does not reach is not a parameter the data paid for: a
# break in the last period of the sample has one scale and no decay, whatever
# the model asked for. In the ordinary case, a break with room after it, this is
# volatility_scales + 1.
.volatility_break_nparams <- function(object) {
  period <- object[["model"]][["volatility_break"]]
  if (is.null(period)) {
    return(0)
  }
  tt <- nrow(object[["data"]][["train"]][["y"]])
  scales <- as.numeric(object[["model"]][["volatility_scales"]])
  after <- tt - period + 1
  min(scales, after) + (after > scales)
}
