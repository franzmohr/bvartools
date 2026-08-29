#' Printing Model Information
#' 
#' print method for objects of class 'modellist'.
#' 
#' @param x an object of class 'modellist'.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
print.modellist <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  result <- do.call("rbind", lapply(x, get_model_specifications))
  names(result)[which(names(result) == "type")] <- "Type"
  names(result)[which(names(result) == "varsel")] <- "Variable selection"
  
  # Drop columns, where all entries of a column are the same
  pos <- NULL
  for (i in names(result)) {
    if (!i %in% c("Type", "Variable selection", "k", "p")) {
      if (all(result[, i] == result[1, i])) {
        pos <- append(i, pos)
      } 
    }
  }
  if (!is.null(pos)) {
    result <- result[, -which(names(result) %in% pos)]
  }
  
  print(result, digits = digits, ...)
  invisible(result)
}