#' Choose Best Model
#' 
#' Chooses the best model according the selection criteria in an object of class
#' 'selcritlist'.
#' 
#' @param object object of class 'selcritlist', usually, a result of a call
#' to \code{\link{selection_criteria}}.
#' @param criterion the selection criterion that should be used. Available choices
#' are \code{"LL"}, \code{"AIC"}, \code{"BIC"}, \code{"HQ"},
#' \code{"WAIC"} (default), \code{"LOOIC"}, \code{"LPL"} and, for the
#' discounted models, \code{"LML"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' If argument \code{criterion} is "LL", "LPL" or "LML", the model with the
#' maximum value is chosen, otherwise, the model with the minimum value.
#'
#' Which criterion to use depends on the models that are compared. See
#' \code{\link{selection_criteria}}, whose details set out what each of them
#' penalises and when a count of parameters ceases to describe a model. The
#' default is \code{"WAIC"}, because it is the criterion that stays defined
#' across the models this package estimates: it penalises by the flexibility
#' the fit used rather than by a count of parameters, which describes neither a
#' model whose coefficients or variances follow a state equation nor one whose
#' prior shrinks them. \code{"AIC"}, \code{"BIC"} and \code{"HQ"} remain the
#' criteria for the lag order of a model with constant coefficients and a weak
#' prior, which is the case they are derived for. \code{"LPL"}, the log
#' predictive likelihood of an expanding window exercise with
#' \code{\link{add_predictive_loglik}}, is the criterion for the rank of a VEC
#' model whose coefficients or variances follow a state equation.
#'
#' \code{"LML"} is the one criterion of the discounted models and the only one
#' they carry. It is the exact log marginal likelihood of the sample, read off
#' the file rather than estimated, so a grid over the rank, the cointegration
#' matrix, the lag order or the two discounts is compared by it directly. It is
#' not comparable with \code{"LL"}, which conditions on the parameters where
#' this integrates them out.
#'
#' @return An integer giving the position of the best model in the list provided
#' in argument \code{object}. Where several models attain the best value
#' exactly, the positions of all of them, in increasing order.
#' 
#' @examples
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 1:3, deterministic = "const",
#'                           iterations = 50, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' model <- add_posterior_coefficients(model)
#' 
#' # Obtain posterior log-likelihoods
#' model <- add_posterior_loglik(model)
#' 
#' # Calculate selection criteria
#' sel <- selection_criteria(model)
#' 
#' # View results
#' sel
#' 
#' # Choose best model according to WAIC, the default
#' choose_best_model(sel)
#' 
#' # Choose best model according to AIC
#' choose_best_model(sel, criterion = "AIC")
#' 
#' @family model comparison
#' @export
#' @method choose_best_model selcritlist
choose_best_model.selcritlist <- function(object, criterion = "WAIC", ...) {

  # One number per model. The forecast error criteria -- FE, AFE, RSFE -- hold
  # one per variable and horizon, which stacked into a matrix made which()
  # below return a position in that matrix rather than a model.
  allowed <- c("LL", "AIC", "BIC", "HQ", "WAIC", "LOOIC", "LPL", "LML")
  if (!is.character(criterion) || length(criterion) != 1 || !criterion %in% allowed) {
    stop("Argument 'criterion' must be one of ", paste0("\"", allowed, "\"", collapse = ", "),
         ". Forecast error criteria hold a value per variable and horizon and name no ",
         "single best model.")
  }

  # Models, which do not contain the criterion, cannot be chosen, but their
  # positions in 'object' are maintained
  res <- vapply(object, function(y) {
    if (is.null(y[[criterion]])) {
      return(NA_real_)
    }
    value <- y[[criterion]][, "mean"]
    if (length(value) != 1) {
      stop("Criterion '", criterion, "' holds more than one value per model.", call. = FALSE)
    }
    as.numeric(value)
  }, numeric(1), USE.NAMES = FALSE)

  if (all(is.na(res))) {
    stop("None of the models in argument 'object' contains criterion '", criterion, "'.")
  }

  # LML, the exact log marginal likelihood of a discounted model, is maximised
  # like the other two densities rather than minimised like the deviance-based
  # criteria.
  if (criterion %in% c("LL", "LPL", "LML")) {
    pos <- which(res == max(res, na.rm = TRUE))
  } else {
    pos <- which(res == min(res, na.rm = TRUE))
  }
  
  return(pos)
}
