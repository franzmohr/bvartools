#' Printing Model Information
#' 
#' print method for objects of class 'bvecmodel'.
#' 
#' @param x an object of class 'bvecmodel'.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
print.bvecmodel <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  result <- get_model_specifications(x)
  names(result)[names(result) == "type"] <- "Type"
  names(result)[names(result) == "varsel"] <- "Variable selection"
  
  if (all(result[, "m"] == 0)) {
    result <- result[, !names(result) %in% c("m", "s")]
  }
  
  print(result, digits = digits, ...)

  # The autocorrelation of the cointegration state equation belongs to the
  # specification of a time varying model, and whether it is fixed or drawn
  # showed up nowhere before the posterior was looked at directly.
  beta_prior <- x[["priors"]][["beta"]]
  if (isTRUE(x[["model"]][["tvp"]]) && !is.null(beta_prior[["rho"]])) {
    if (!is.null(beta_prior[["rho_min"]])) {
      cat("\nAutocorrelation of the cointegration state equation (rho): drawn from a ",
          "uniform prior on (", format(beta_prior[["rho_min"]]), ", ",
          format(beta_prior[["rho_max"]]), "), starting at ",
          format(beta_prior[["rho"]]), "\n", sep = "")
    } else {
      cat("\nAutocorrelation of the cointegration state equation (rho): fixed at ",
          format(beta_prior[["rho"]]), "\n", sep = "")
    }
  }

  invisible(result)
}