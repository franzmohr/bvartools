#' Generate a VECX model
#' 
#' Produces the model matrix of a VEC model with exogenous variables.
#' 
#' @param data containing the data and specification of a country model.
#' 
#' @return List containing the final model data and specifications.
#' 
#' @export
gen.vecx <- function(data){
  n <- ncol(data$x)
  n.star <- ncol(data$x.star)
  r <- data$specs$rank$rank
  
  level.domestic <- data$x
  level.star <- data$x.star
  ect <- cbind(level.domestic, level.star)
  names.ect <- c(paste(data$specs$domestic.variables, ".l1", sep = ""),
                 paste("s.", data$specs$foreign.variables, ".l1", sep = ""))
  
  global <- !is.na(data$specs$global.variables)
  n.global <- 0
  s <- NA
  if (global){
    s <- data$specs$lags$lags$global
    if (!is.na(s)){
      n.global <- length(data$specs$global.variables)
      level.global <- data$x.g
      ect <- cbind(ect, level.global)
      names.ect <- c(names.ect, paste(data$specs$global.variables, ".l1", sep = "")) 
    }
  }
  
  # Add restricted deterministic terms to error correction term
  det <- data$specs$deterministic.terms
  n.det.r <- 0
  if (!is.null(det)){
    case <- det
    if (class(case) == "character"){
      if (case == "II"){
        ect <- cbind(ect, 1)
        names.ect <- c(names.ect, "const")
        n.det.r <- 1
      }
      if (case=="IV"){
        ect <- cbind(ect, 1:nrow(ect))
        names.ect <- c(names.ect, "trend")
        n.det.r <- 1
      }
    }
  }
  n.ect <- dim(ect)[2]
  ect <- stats::lag(ect, -1)
  
  p <- data$specs$lags$lags$domestic - 1
  q <- data$specs$lags$lags$foreign - 1
  
  diff.domestic <- diff(level.domestic)
  total <- cbind(diff.domestic, ect)
  names.total <- c(paste("d.", data$specs$domestic.variables, sep = ""), names.ect)
  
  # Lags of domestic variables
  n.d <- 0
  if (p > 0){
    for (i in 1:p){
      total <- cbind(total, stats::lag(diff.domestic, -i))
      names.total <- c(names.total, paste("d.", data$specs$domestic.variables, ".l", i, sep = ""))
      n.d <- n.d + n
    }
  }
  
  # Lags of star variables
  diff.star <- diff(data$x.star)
  total <- cbind(total, diff.star)
  names.total <- c(names.total, paste("d.s.", data$specs$foreign.variables, sep = ""))
  n.s <- n.star
  if (q > 0){
    for (i in 1:q){
      total <- cbind(total, stats::lag(diff.star, -i))
      names.total <- c(names.total, paste("d.s.", data$specs$foreign.variables, ".l", i, sep = ""))
      n.s <- n.s + n.star
    }
  }
  
  # Lags of global variables
  n.g <- 0
  if (global){
    if (!is.na(s)){
      diff.global <- diff(data$x.g)
      total <- cbind(total, diff.global)
      names.total <- c(names.total, paste("d.", data$specs$global.variables, sep = ""))
      n.g <- n.global
      if ((s - 1) > 0){
        for (i in 1:(s - 1)){
          total <- cbind(total, stats::lag(diff.global, -i))
          names.total <- c(names.total, paste("d.", data$specs$global.variables, ".l", i, sep = ""))
          n.g <- n.g + n.global
        }
      }
    }
  }
  
  # Add unrestrited deterministic terms
  n.det.ur <- 0
  if (!is.null(det)){
    if (case == "III"){
      total <- cbind(total, 1)
      names.total <- c(names.total, "const")
      n.det.ur <- 1
    }
    if (case == "IV"){
      total <- cbind(total, 1)
      names.total <- c(names.total, "const")
      n.det.ur <- 1
    }
    if (case == "V"){
      total <- cbind(total, 1, 1:nrow(total))
      names.total <- c(names.total, "const", "trend")
      n.det.ur <- 2
    }
  }
  
  total <- stats::na.omit(total)
  names(total) <- names.total
  
  # Final data preparations
  used.t <- seq(to = nrow(total), length.out = nrow(data$x) - data$specs$lags$maximum.lag)
  y <- t(total[used.t, 1:n])
  ect <- t(total[used.t, n + 1:(n.ect)])
  x <- t(total[used.t, -(1:(n + n.ect))])
  
  result <- list(y = y, ect = ect, x = x,
                 domestic = c("dim" = n, "lag" = as.numeric(p), "total" = n.d),
                 foreign = c("dim" = n.star, "lag" = as.numeric(q), "total" = n.s),
                 global = c("dim" = n.global, "lag" = as.numeric(s), "total" = n.g),
                 deterministic = c("restricted" = n.det.r, "unrestricted" = n.det.ur),
                 rank = r)
  return(result)
}