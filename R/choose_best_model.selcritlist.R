#' Choose Best Model
#' 
#' Chooses the best model according the selection criteria in an object of class
#' 'selcritlist'.
#' 
#' @param object object of class 'selcritlist', usually, a result of a call
#' to \code{\link{selection_criteria}}.
#' @param criterion the selection criterion that should be used. Available choices
#' are \code{"LL"}, \code{"AIC"}, \code{"BIC"}, \code{"HQ"},
#' \code{"WAIC"} (default) and \code{"LOOIC"}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details
#' If argument \code{criterion} is "LL", the model with the maximum value is chosen,
#' otherwise, the model with the minimum value.
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
#' prior, which is the case they are derived for.
#' 
#' @returns An integer giving the position of the best model in the list provided
#' in argument \code{object}.
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
#' @export
#' @method choose_best_model selcritlist
choose_best_model.selcritlist <- function(object, criterion = "WAIC", ...) {
  
  # Models, which do not contain the criterion, cannot be chosen, but their
  # positions in 'object' are maintained
  res <- lapply(object, function(y, criterion) {
    if (is.null(y[[criterion]])) {
      return(NA_real_)
    }
    y[[criterion]][, "mean"]
  }, criterion = criterion)
  res <- do.call("rbind", res)

  if (all(is.na(res))) {
    stop("None of the models in argument 'object' contains criterion '", criterion, "'.")
  }

  if (criterion == "LL") {
    pos <- which(res == max(res, na.rm = TRUE))
  } else {
    pos <- which(res == min(res, na.rm = TRUE))
  }
  
  return(pos)
}
