#' Align Observations Across Models
#' 
#' Restricts each model in an object of class 'modellist' to the set of
#' observations common to all models, ensuring that comparisons are computed
#' on the same underlying sample.
#' 
#' @param object an object of class 'modellist'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'modellist'.
#' 
#' @examples
#' 
#' # Load data
#' data(e6)
#' 
#' # Create multiple VAR models
#' model_1 <- create_bvarmodel(diff(e6), p = 0:2, deterministic = "const")
#' 
#' # Create multiple VEC models
#' model_2 <- create_bvecmodel(e6, p = 1, r = 0:1, const = "unrestricted")
#' 
#' # Combine the models in one object
#' model <- combine_models(model_1, model_2)
#' 
#' model <- align_model_obs(model)
#' 
#' @export
align_model_obs.modellist <- function(object, ...) {
  
  # Get sample sizes
  avail <- lapply(object, function(x) {stats::tsp(x[["data"]][["train"]][["y"]])})
  avail <- do.call("rbind", avail)
  
  min_date <- max(avail[, 1]) # First period overall is latest in start
  max_date <- min(avail[, 2]) # First period overall is earlierst in end
  
  object <- stats::window(x = object, start = min_date, end = max_date, ...)
  
}