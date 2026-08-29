#' Combine Models
#' 
#' Convenience function that combines multiple models into one object.
#' 
#' @param ... One or multiple lists that contain either model objects or
#' posterior objects.
#' 
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
#' @export
combine_models <- function(...) {
  
  input <- list(...)
  
  classes <- unlist(lapply(input, class))
  classes <- classes[-which(classes == "list")]
  
  supported_models <- c("bvarmodel", "bvecmodel", "klgs2010", "modellist",
                        "expandingwindow")
  supported_total <- supported_models
  if (!all(classes %in% supported_total)) {
    stop("The class for at least one input object could not be recognised.")
  }
  
  result <- NULL
  
  for (i in 1:length(input)) {
    if (any(c("bvarmodel", "bvecmodel", "expandingwindow") %in% class(input[[i]]))) {
      result <- c(result, list(input[[i]]))
    }
    if (any(c("modellist") %in% class(input[[i]]))) {
      result <- c(result, input[[i]])
    }
  }
  class(result) <- append("modellist", class(result))
  
  return(result)
  
}