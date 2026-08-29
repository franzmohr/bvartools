#' Printing Model Information
#' 
#' print method for objects of class 'bvarmodel'.
#' 
#' @param x an object of class 'bvarmodel'.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
print.bvarmodel <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  result <- get_model_specifications(x)
  names(result)[names(result) == "type"] <- "Type"
  names(result)[names(result) == "varsel"] <- "Variable selection"
  
  if (all(result[, "m"] == 0)) {
    result <- result[, !names(result) %in% c("m", "s")]
  }
  
  print(result, digits = digits, ...)
  invisible(result)
}