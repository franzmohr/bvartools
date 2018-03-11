#' Weight Matrices
#' 
#' Produces a weight matrices for each country model.
#' 
#' @param weight.data a named matrix of trade weights.
#' @param index a data frame produced by the function \code{\link{global_series}}.
#' @param country.specs a list of model specifications as produced by \code{\link{country_specifications}}.
#' 
#' @return A list containing the weight matrices for each country model.
#' 
#' @export
weight_matrices <- function(weight.data, index, country.specs){
  tv <- !is.na(dim(weight.data)[3])
  
  # Check data availability for each country
  countries <- names(country.specs)
  variables <- unique(index[, 2])
  n.countries <- length(countries)
  n.variables <- length(variables)
  var.exists <- as.data.frame(matrix(NA, n.countries, n.variables))
  row.names(var.exists) <- countries
  names(var.exists) <- variables
  for (i in countries){
    for (j in variables){
      var.exists[i, j] <- is.element(j, index[index[,1] == i, 2])
    }
  }
  k <- dim(index)[1]
  
  W <- c()
  W.names <- c()
  for (i in countries){
    n.endo <- length(index[index[,1] == i, 1])
    star.vars <- country.specs[[i]]$foreign.variables
    n.star <- length(star.vars)
    if (tv){
      W.i <- array(0, dim = c(n.endo + n.star, k, dim(weight.data)[3]))
      dimnames(W.i) <- list(c(index[index[, 1] == i, 2], paste("s.", star.vars, sep = "")),
                            index[, 2], dimnames(weight.data)[[3]])
      for (l in 1:dim(weight.data)[3]) {
        W.i[1:n.endo, which(index[, 1] == i), l] <- diag(1, n.endo)
        for (j in star.vars){
          w.temp <- weight.data[i,,l]
          # Set all values to zero, where variable is not available
          w.temp[!var.exists[, j]] <- 0
          # Extract weights for available values
          exist.weights <- weight.data[i, var.exists[, j], l]
          # Make sum over weights of not available values
          non.exist.weights <- sum(weight.data[i, !var.exists[, j], l])
          # Recalculate weights
          w.temp[var.exists[, j]] <- exist.weights/(1 - non.exist.weights)
          # Add weights to country's weight matrix
          pos.x <- n.endo + which(star.vars == j)
          pos.y <- which(index[, 2] == j)
          W.i[pos.x, pos.y, l] <- w.temp[var.exists[, j]]
        } 
      }
    } else {
      W.i <- matrix(0, n.endo + n.star, k)
      dimnames(W.i) <- list(c(index[index[, 1] == i, 2], paste("s.", star.vars, sep = "")),
                            index[, 2])
      W.i[1:n.endo, which(index[, 1] == i)] <- diag(1, n.endo)
      for (j in star.vars){
        w.temp <- weight.data[i,]
        # Set all values to zero, where variable is not available
        w.temp[!var.exists[, j]] <- 0
        # Extract weights for available values
        exist.weights <- weight.data[i, var.exists[, j]]
        # Make sum over weights of not available values
        non.exist.weights <- sum(weight.data[i, !var.exists[, j]])
        # Recalculate weights
        w.temp[var.exists[, j]] <- exist.weights/(1 - non.exist.weights)
        # Add weights to country's weight matrix
        pos.x <- n.endo + which(star.vars == j)
        pos.y <- which(index[, 2] == j)
        W.i[pos.x, pos.y] <- w.temp[var.exists[, j]]
      }
    }

    W <- c(W, list(W.i))
    W.names <- c(W.names, i)
  }
  names(W) <- W.names
  return(W)
}