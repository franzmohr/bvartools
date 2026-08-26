#' Align Observations Across Models
#'  
#' Generic function used to restrict each model in a list to the set of
#' observations common to all models, ensuring that comparisons are computed
#' on the same underlying sample.
#' 
#' @param object a list of models.
#' @param ... arguments passed forward to method.
#' 
#' @export
align_model_obs <- function (object, ...) {
 UseMethod("align_model_obs")
}