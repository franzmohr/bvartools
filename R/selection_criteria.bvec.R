#' Model Selection Criteria
#'
#' Calculates model selection criteria for an object of class \code{"bvec"}.
#'
#' @param object an object of class \code{"bvec"}, usually, the
#' output of a call to \code{\link{draw_posterior}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list.
#' 
#' @export
selection_criteria.bvec <- function(object, ...){
  
  # Save r for later
  r <- object[["specifications"]][["rank"]]
  
  object <- bvec_to_bvar(object)
  
  x <- selection_criteria(object)
  
  # Add r to summary table
  pos_specs <- 1:(ncol(x[["summary"]]) - 4)
  names_smry <- names(x[["summary"]])[pos_specs]
  smry <- x[["summary"]][, pos_specs]
  smry <- cbind(smry, data.frame(r = r))
  names(smry) <- c(names_smry, "r")
  x[["summary"]] <- cbind(smry, x[["summary"]][, c(ncol(x[["summary"]]) - 4 + 1:4)])
  
  return(x)
}