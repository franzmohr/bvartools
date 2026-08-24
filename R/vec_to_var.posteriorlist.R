#' Transform a VEC Model to a VAR in Levels
#' 
#' The elements of an object of class 'posteriorlist', which are VEC model,
#' are transformed to a VAR in level representation.
#' 
#' @param object an object of class 'posteriorlist'.
#' @param ... arguments passed forward to method.
#' 
#' @return A list of class 'posteriorlist'.
#' 
#' 
#' @export
#' @method vec_to_var posteriorlist
vec_to_var.posteriorlist <- function(object, ...){
  
  object <- lapply(object, vec_to_var, ...)
  
  class(object) <- append("posteriorlist", class(object))
  
  return(object)
}
