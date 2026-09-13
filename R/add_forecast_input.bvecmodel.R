#' Add Forecast Input Data
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param object an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since forecasts of a VEC model are obtained
#' from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{\link{add_forecast_input}} to the resulting 'bvarmodel'.
#'
#' @family posterior simulation
#' @export
#' @method add_forecast_input bvecmodel
add_forecast_input.bvecmodel <- function(object, ...){
  .use_vec_to_var("add_forecast_input")
}



#' Add Forecasts
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param object an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since forecasts of a VEC model are obtained
#' from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{\link{add_posterior_forecasts}} to the resulting 'bvarmodel'.
#'
#' @family posterior simulation
#' @export
#' @method add_posterior_forecasts bvecmodel
add_posterior_forecasts.bvecmodel <- function(object, ...){
  .use_vec_to_var("add_posterior_forecasts")
}



#' Predict Method for Objects of Class bvecmodel
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param object an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since forecasts of a VEC model are obtained
#' from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{predict} to the resulting 'bvarmodel'.
#'
#' @family posterior simulation
#' @export
#' @method predict bvecmodel
predict.bvecmodel <- function(object, ...){
  .use_vec_to_var("predict")
}



#' Impulse Response Function
#'
#' Guard method for objects of class 'bvecmodel'.
#'
#' @param x an object of class 'bvecmodel'.
#' @param ... additional arguments.
#'
#' @return Nothing. The method raises an error, since impulse responses of a VEC model are
#' obtained from its VAR representation: apply \code{\link{vec_to_var}} first and then
#' \code{\link{irf}} to the resulting 'bvarmodel'.
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
#' \code{\link{fevd}} to the resulting 'bvarmodel'.
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



# The error the functions above raise. The analysis of a VEC model goes through
# its VAR representation, so the functions that analyse or forecast a fitted
# model have no method of their own for it. Without these guards all but
# add_forecast_input() stopped with R's "no applicable method", which says
# nothing about what to do instead.
.use_vec_to_var <- function(fun) {
  stop("'", fun, "' does not work directly on a 'bvecmodel' object.\n",
       "Use 'vec_to_var()' first and then use '", fun, "' on the\n",
       "resulting 'bvarmodel' object.",
       call. = FALSE)
}
