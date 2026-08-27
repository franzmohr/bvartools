#' Choose Best Model
#' 
#' Chooses the best model according the selection criteria in an object of class
#' 'selcritlist'.
#' 
#' @param object object of class 'selcritlist', usually, a result of a call
#' to \code{\link{selection_criteria}}.
#' @param criterion the selection criterion that should be plotted. Available choices
#' are \code{"LL"}, \code{"AIC"}, \code{"BIC"} (default), \code{"HQ"}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details
#' If argument \code{criterion} is "LL", the model with the maximum value is chosen,
#' otherwise, the model with the minimum value.
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
#' # Choose best model according to AIC
#' choose_best_model(sel, criterion = "AIC")
#' 
#' # Choose best model according to BIC
#' choose_best_model(sel, criterion = "BIC")
#' 
#' @export
#' @method choose_best_model selcritlist
choose_best_model.selcritlist <- function(object, criterion = "BIC", ...) {
  
  res <- lapply(object, function(y, criterion) {y[[criterion]][, "mean"]}, criterion = criterion)
  res <- do.call("rbind", res)
  
  if (criterion == "LL") {
    pos <- which(res == max(res))
  } else {
    pos <- which(res == min(res)) 
  }
  
  return(pos)
}
