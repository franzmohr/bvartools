#' Inclusion Probabilities
#' 
#' Calculates the inclusion probablities of variable in a country model, where BVS or SSVS were used.
#' 
#' @param data a list containing the solved GVAR model as produced by \code{\link{solve_gvar}}.
#' @param all a logical. If \code{TRUE}, the contemporaneous effects of all estimated models will be printed. Otherwise, only
#' those models will be considered, which enter the global model.
#' 
#' @return A list of matrices.
#' 
#' @export
get_inclusion_probabilities <- function(data, all = FALSE, decimal = NULL){
  result <- lapply(data$country.models, function(x) {return(x$coefs$A_p_include)})
  if (all(is.null(unlist(result)))) {
    stop("No additional shinkage methods were used.")
  }
  if (all) {
    pos_final <- 1:length(data$country.models)
  } else {
    pos_final <- which(data$GVAR$specs$criteria[, "Maximum"])
  }
  final <- c()
  for (i in pos_final) {
    if (is.null(decimal)) {
      final <- c(final, list(result[[i]])) 
    } else {
      final <- c(final, list(round(result[[i]], decimal)))
    }
  }
  names(final) <- names(result)[pos_final]
  result <- final
  
  return(result)
}