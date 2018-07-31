#' Generate a VARX model
#' 
#' Produces the model matrix of a VAR model with exogenous variables.
#' 
#' @param data containing the data and specification of a country model.
#' 
#' @return List containing the final model data and specifications.
#' @export
gen.varx <- function(data){
  n <- ncol(data$x)
  n.star <- ncol(data$x.star)
  
  x <- data$x
  x.star <- data$x.star
  
  global <- !is.na(data$specs$global.variables)
  n.global <- 0
  p.global <- NA
  if (global){
    p.global <- data$specs$lags$lags$global
    if (!is.na(s)){
      n.global <- length(data$specs$global.variables)
      x.global <- data$x.g
    }
  }
  
  p <- data$specs$lags$lags$domestic
  p.star <- data$specs$lags$lags$foreign
  
  total <- x
  names.total <- data$specs$domestic.variables
  # Lags of domestic variables
  n.d <- 0
  if (p > 0){
    for (i in 1:p){
      total <- cbind(total, stats::lag(x, -i))
      names.total <- c(names.total, paste(data$specs$domestic.variables, ".l", i, sep = ""))
      n.d <- n.d + n
    }
  }
  
  # Lags of star variables
  total <- cbind(total, x.star)
  names.total <- c(names.total, paste("s.", data$specs$foreign.variables, sep = ""))
  n.s <- n.star
  if (p.star > 0){
    for (i in 1:p.star){
      total <- cbind(total, stats::lag(x.star, -i))
      names.total <- c(names.total, paste("s.", data$specs$foreign.variables, ".l", i, sep = ""))
      n.s <- n.s + n.star
    }
  }
  
  # Lags of global variables
  n.g <- 0
  if (global){
    if (!is.na(p.global)){
      total <- cbind(total, x.global)
      names.total <- c(names.total, data$specs$global.variables)
      n.g <- n.global
      if ((p.global - 1) > 0){
        for (i in 1:p.global){
          total <- cbind(total, stats::lag(x.global, -i))
          names.total <- c(names.total, paste(data$specs$global.variables, ".l", i, sep = ""))
          n.g <- n.g + n.global
        }
      }
    }
  }
  
  # Add deterministic terms
  case <- data$specs$deterministic.terms
  n.det <- 0
  if (!is.null(case)){
    if (case == "II" | case == "III"){
      total <- cbind(total, 1)
      names.total <- c(names.total, "const")
      n.det <- 1
    }
    if (case == "IV" | case == "V"){
      total <- cbind(total, 1, 1:nrow(total))
      names.total <- c(names.total, "const", "trend")
      n.det <- 2
    }
  }
  
  total <- stats::na.omit(total)
  dimnames(total)[[2]] <- names.total
  
  # Final data preparations
  used.t <- seq(to = nrow(total), length.out = nrow(data$x) - data$specs$lags$maximum.lag)
  y <- t(total[used.t, 1:n])
  x <- t(total[used.t, -(1:n)])
  t <- ncol(y)
  
  result <- list(y = y, x = x,
                 domestic = c("dim" = n, "lag" = as.numeric(p), "total" = n.d),
                 foreign=c("dim" = n.star, "lag" = as.numeric(p.star), "total" = n.s),
                 global=c("dim" = n.global, "lag" = as.numeric(p.global), "total" = n.g),
                 deterministic = n.det)
  return(result)
}