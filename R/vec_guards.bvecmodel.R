# Methods of 'bvecmodel' for the analysis functions that go through its VAR
# representation. They used to live in add_forecast_input.bvecmodel.R.

#' Impulse Response Function
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param x an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since impulse responses of a VEC model are
#' obtained from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{\link[bayests]{irf}} to the resulting 'bvarmodel'.
#'
#' @family post-estimation analysis
#' @export
#' @method irf bvecmodel
irf.bvecmodel <- function(x, ...){
  .use_vec_to_var("irf")
}



#' Forecast Error Variance Decomposition
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param x an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since variance decompositions of a VEC model
#' are obtained from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{\link[bayests]{fevd}} to the resulting 'bvarmodel'.
#'
#' @family post-estimation analysis
#' @export
#' @method fevd bvecmodel
fevd.bvecmodel <- function(x, ...){
  .use_vec_to_var("fevd")
}



#' Spillover Index
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param object an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since spillover indices of a VEC model are
#' obtained from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{\link{spillover}} to the resulting 'bvarmodel'.
#'
#' @family post-estimation analysis
#' @export
#' @method spillover bvecmodel
spillover.bvecmodel <- function(object, ...){
  .use_vec_to_var("spillover")
}



# The error the guards above raise. Impulse responses, variance decompositions
# and spillovers of a VEC model go through its VAR representation, so those
# functions have no method of their own for it. Without these guards they
# stopped with R's "no applicable method", which says nothing about what to do
# instead. Forecasts are the exception: they are simulated from the VEC model's
# own draws, in levels, by the methods above.
.use_vec_to_var <- function(fun) {
  stop("'", fun, "' does not work directly on a 'bvecmodel' object.\n",
       "Use 'vec_to_var()' first and then use '", fun, "' on the\n",
       "resulting 'bvarmodel' object.",
       call. = FALSE)
}
