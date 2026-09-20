#' Selection Criteria
#'
#' Calculates the criteria that a pointwise log-likelihood and a scored forecast
#' are enough for, whatever the model is.
#'
#' @param object any object with \code{posterior$loglik}, \code{posterior$forecast$loglik}
#' or both. A model of a class this package has no method for -- a dynamic factor
#' model of \pkg{dfmtools}, say -- reaches this one.
#' @param ci the width of the credible bands, a value between 0 and 1. Defaults
#' to 0.95.
#' @param ... additional arguments.
#'
#' @details
#' The criteria of \code{\link{selection_criteria.bvarmodel}} fall into two
#' groups. \code{AIC}, \code{BIC} and \code{HQ} charge a model for its size, so
#' they need a count of its free parameters, which is a property of the model and
#' cannot be read off its draws. \code{FE}, \code{AFE} and \code{RSFE} need the
#' forecast errors, which need the variables to be named and paired up. Neither
#' group is available here.
#'
#' What is left needs nothing but the draws themselves:
#' \itemize{
#'  \item \code{LL}, the log-likelihood of the whole sample per draw, summed over
#' the periods of \code{posterior$loglik};
#'  \item \code{WAIC} and \code{LOOIC}, which penalise by the flexibility the fit
#' used rather than by a count of parameters -- see \code{\link{selection_criteria}}
#' for why that is what makes constant, time varying and stochastic volatility
#' specifications comparable at all;
#'  \item \code{LPL}, the log predictive likelihood of
#' \code{posterior$forecast$loglik}, the score of a forecast against what its
#' horizon realised.
#' }
#'
#' The periods of the \code{"terms"} attribute of \code{LPL} are numbered from
#' one rather than dated: this method knows nothing of where a model keeps its
#' sample, so it cannot say when the scored periods were. A class that can say
#' should write its own method, as \code{\link{selection_criteria.bvarmodel}}
#' does.
#'
#' @return A list of class 'selcrit', which also inherits the class of
#' \code{object}, with the element \code{model} and one data frame per criterion
#' that the draws supported, each with the columns \code{mean}, \code{median},
#' \code{qlower} and \code{qupper}. \code{\link{choose_best_model}} ranks a list
#' of them and \code{print} shows them.
#'
#' @examples
#'
#' # Any object with the draws will do, whatever produced them.
#' set.seed(7)
#' model <- list(model = list(k = 2),
#'               posterior = list(loglik = matrix(stats::rnorm(200, -3), 50, 4),
#'                                forecast = list(loglik = matrix(stats::rnorm(150, -3), 50, 3))))
#'
#' criteria <- selection_criteria(model)
#' names(criteria)
#' criteria[["LPL"]][["mean"]]
#'
#' @family model comparison
#' @export
#' @method selection_criteria default
selection_criteria.default <- function(object, ci = 0.95, ...) {

  if (ci < 0 | ci > 1) {
    stop("Argument 'ci' must be between 0 and 1.")
  }
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low

  loglik <- object[["posterior"]][["loglik"]]
  score <- object[["posterior"]][["forecast"]][["loglik"]]

  if (is.null(loglik) & is.null(score)) {
    stop("Object must contain posterior draws of the log-likelihood in posterior$loglik, ",
         "a scored forecast in posterior$forecast$loglik, or both.")
  }

  result <- NULL
  result[["model"]] <- object[["model"]]

  if (!is.null(loglik)) {

    total <- rowSums(as.matrix(loglik))
    ll <- data.frame("mean" = mean(total),
                     "median" = stats::median(total),
                     "qlower" = stats::quantile(total, probs = ci_low),
                     "qupper" = stats::quantile(total, probs = ci_high))
    row.names(ll) <- NULL
    result[["LL"]] <- ll

    waic <- .waic_data_frame(loglik, ci_low, ci_high)
    if (!is.null(waic)) {
      row.names(waic) <- NULL
      result[["WAIC"]] <- waic
    }

    looic <- .loo_data_frame(loglik, ci_low, ci_high)
    if (!is.null(looic)) {
      row.names(looic) <- NULL
      result[["LOOIC"]] <- looic
    }
  }

  if (!is.null(score)) {
    score <- as.matrix(score)
    result[["LPL"]] <- .lpl_entry(lapply(seq_len(ncol(score)), function(i) score[, i]),
                                  seq_len(ncol(score)), ci_low, ci_high)
  }

  attr(result, "ci") <- c(paste0(ci_low * 100, "%"), paste0(ci_high * 100, "%"))
  class(result) <- c("selcrit", class(object))

  return(result)
}
