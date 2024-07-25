#' @include selection_criteria.R
#'
#' @export
#' @rdname selection_criteria
print.selcritlist <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  res <- lapply(x, function(y) {y[["summary"]]})
  
  # # Harmonise table formats so that "rbind" works
  # unique_names <- unique(unlist(lapply(object, names)))
  # for (i in 1:length(object)) {
  #   if (all(unique_names %in% names(object[[i]]))) {
  #     next
  #   } else {
  #     add_columns <- unique_names[which(!unique_names %in% names(object[[i]]))]
  #     for (j in add_columns) {
  #       temp <- matrix(NA, 1, length(add_columns))
  #       dimnames(temp) <- list(NULL, add_columns)
  #       object[[i]] <- cbind(object[[i]], as.data.frame(temp))
  #     }
  #   }
  # }
  
  res <- do.call("rbind", res)
  print(res, digits = digits, ...)
  invisible(res)
} 

#' @include selection_criteria.R
#'
#' @export
#' @rdname selection_criteria
print.selcrit <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  print(x[["summary"]], digits = digits, ...)
  invisible(x)
} 