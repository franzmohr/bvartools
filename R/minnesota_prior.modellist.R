#' Minnesota Prior
#'
#' Forwards the elements of an object of class 'modellist' forward to its
#' respective method.
#'
#' @param object a list of class 'modellist'.
#' @param ... further arguments passed to or from other methods.
#' 
#' The same for an object of class 'expandingwindow', each window getting the
#' prior of its own sample.
#'
#' @return A list of class 'modellist'.
#' 
#' @export
#' @method minnesota_prior modellist
minnesota_prior.modellist <- function(object, ...){
  
  for (i in 1:length(object)) {
    object[[i]] <- minnesota_prior(object[[i]], ...)
  }
  
  return(object)
}