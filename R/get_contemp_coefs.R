#' Contemporaneous Coefficients
#' 
#' Calculates statistics for draws of contemporaneous coefficients of foreign variables on their domestic counterpart
#' for each country model, which was used to construct the global model.
#' 
#' @param data a list containing the solved GVAR model as produced by \code{\link{solve_gvar}}.
#' @param type a character specifying the type of statistic that should be calculated.
#' Possible option are "mean" (default), "median" and "sd".
#' @param t a numeric specifying the time index of the coefficients. Only applicable if the country models were
#' estimated with time-varying parameters or stochastic volatility.
#' If \code{NULL} the latest or only available period will be used.
#' @param all a logical. If \code{TRUE}, the contemporaneous effects of all estimated models will be printed. Otherwise, only
#' those models will be considered, which enter the global model.
#' 
#' @return A data frame.
#' 
#' @export
get_contemp_coefs <- function(data, type = "mean", t = NULL, all = FALSE){
  d.vars <- unique(unlist(lapply(data$Country.models, function(x){return(x$specs$domestic.variables)})))
  s.vars <- unique(unlist(lapply(data$Country.models, function(x){return(x$specs$foreign.variables)})))
  vars <- s.vars[which(is.element(s.vars, d.vars))]
  
  result <- as.data.frame(matrix(NA, length(data$Country.models), length(vars) + 1))
  names(result) <- c("Country", vars)
  result[, "Country"] <- names(data$Country.models)
  
  if (is.null(t)){
    t <- dim(data$Country.models[[1]]$coefs$A.s.0)[2]
  }
  
  for (j in 1:nrow(result)){
    if (type == "mean") {
      A.s.0 <- apply(data$Country.models[[j]]$coefs$A.s.0[, t,], 1, mean)
    }
    if (type == "median") {
      A.s.0 <- apply(data$Country.models[[j]]$coefs$A.s.0[, t,], 1, median)
    }
    if (type == "sd") {
      A.s.0 <- apply(data$Country.models[[j]]$coefs$A.s.0[, t,], 1, sd)
    }
    names.A0 <- names(A.s.0)
    v <- strsplit(names.A0, split = "_", fixed = TRUE)
    v1 <- unique(unlist(lapply(v, function(x){return(x[1])})))
    v2 <- unique(unlist(lapply(v, function(x){return(x[2])})))
    v2 <- strsplit(v2, split = ".", fixed = TRUE)
    v2 <- unlist(lapply(v2, function(x){return(x[length(x)])}))
    
    A.s.0 <- matrix(A.s.0, length(v1))
    
    pos.v1 <- which(is.element(v1, v2))
    pos.v2 <- which(is.element(v2, v1))
    vars.i <- v2[pos.v2]
    
    for (i in vars.i){
      result[j, which(is.element(vars, i)) + 1]  <- round(A.s.0[which(is.element(v1, i)), which(is.element(v2, i))], 2)
    }
  }
  if (!all) {
    result <- result[which(data$GVAR$specs$criteria[, "Maximum"]),]
  }
  return(result)
}