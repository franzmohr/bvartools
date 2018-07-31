#' Create Regional Series
#' 
#' Combines multiple country-specific time series to regional series.
#' 
#' @param country.data a named list of "zoo" objects containing country data.
#' @param regions a named list of character vectors containing the names of the countries, which belong to a region.
#' The name of the list element will be the name of the region.
#' @param region.weights an object of class "zoo" containing the data used to weight the observations to construct regions.
#' @param period a character vector specifiying the periods in \code{region.weights}, which should be used to construct weights.
#' @param weight.data a matrix or an array of weight data used to generate weight matrices in the GVAR model.
#' 
#' @return A list containing the following components:
#' \item{country.data}{a named list containing the updated country models.}
#' \item{weight.data}{a matrix or an array containing the updated weights.}
#' 
#' @examples
#' \dontrun{
#' data(gvar2016)
#' 
#' country.data <- gvar2016$country.data
#' weight.data <- gvar2016$weight.data
#' region.weights <- gvar2016$region.weights
#' 
#' period <- c("2008", "2009", "2010", "2011", "2012")
#' regions <- list("EuroArea" = c("Austria", "Belgium", "Finland",
#' "France", "Germany", "Italy", "Netherlands", "Spain"))
#' regional.data <- create_regions(country.data, regions, period, region.weights, weight.data)
#' 
#' country.data <- regional.data$country.data
#' weight.data <- regional.data$weight.data
#' }
#' 
#' @export
create_regions <- function(country.data, regions, period, region.weights, weight.data){
  m <- FALSE
  if (class(weight.data) == "matrix") {
    w.names <- dimnames(weight.data)
    weight.data <- array(weight.data, dim = c(dim(weight.data)[1], dim(weight.data)[2], 1))
    dimnames(weight.data) <- list(w.names[[1]], w.names[[2]], "1")
    m <- TRUE
  }
  
  tt <- unique(unlist(lapply(country.data, NROW)))
  if (length(tt) > 1) {stop("Country data must have the same numbers of observations.")}
  
  if ((class(regions) != "list") | is.null(names(regions))) {stop("Object 'regions' must be a named list.")}
  
  if (length(unique(unlist(regions))) < length(unlist(regions))) {
    stop("The same country is not allowed to be in more than one region.") 
  }
  
  vars <- unique(unlist(lapply(country.data, function(x){return(dimnames(x)[[2]])})))
  
  var_exist <- matrix(FALSE, length(country.data), length(vars))
  dimnames(var_exist) <- list(names(country.data), vars)
  for (i in names(country.data)) {
    var_exist[i, dimnames(country.data[[i]])[[2]]] <- TRUE
  }
  
  if (length(period) == 1) {
    t.temp <- as.numeric(time(country.data[[1]]))
    t.avail <- as.numeric(time(region.weights))
    
    t <- matrix(NA, tt, period)
    for (i in 1:tt) {
      if (t.temp[i] <= t.avail[period]) {
        t[i,] <- t.avail[1:period]
      }
      if (t.temp[i] >= t.avail[period]) {
        pos_t <- which(floor(t.temp[i]) == t.avail)
        pos_t <- (pos_t - period + 1):pos_t
        t[i,] <- t.avail[pos_t]
      }
    }
  }
  
  r_names <- names(regions)
  all_r_countries <- unlist(regions)
  names(all_r_countries) <- NULL
  r_tsp <- tsp(country.data[[1]])
  
  n.new <- dim(weight.data)[1] - length(unlist(regions)) + length(regions)
  w.temp <- array(0, dim = c(n.new, n.new, dim(weight.data)[3]))
  w.names <- dimnames(weight.data)[[1]][-which(dimnames(weight.data)[[2]] %in% unlist(regions))]
  w.temp[1:length(w.names), 1:length(w.names) , ] <- weight.data[w.names, w.names, ]
  w.names <- c(w.names, names(regions))
  dimnames(w.temp) <- list(w.names, w.names, dimnames(weight.data)[[3]])
  
  r.temp <- c()
  for (i in 1:length(regions)) {
    vars_r <- apply(var_exist[regions[[i]], ], 2, any)
    vars_r <- dimnames(var_exist)[[2]][vars_r]
    
    r_temp <- ts(matrix(NA, tt, length(vars_r)), start = r_tsp[1], frequency = r_tsp[3])
    dimnames(r_temp)[[2]] <- vars_r
    
    for (j in vars_r) {
      c_temp <- matrix(NA, tt, length(regions[[i]]))
      dimnames(c_temp)[[2]] <- regions[[i]]
      for (k in regions[[i]]) {
        if (var_exist[k , j]) {
          c_temp[, k] <- country.data[[k]][, j] 
        }
      }
      c_temp <- c_temp[, var_exist[regions[[i]], j]]
      
      if (NCOL(c_temp) > 1) {
        if (length(period) == 1) {
          for (k in 1:tt) {
            temp <- colSums(region.weights[which(dimnames(region.weights)[[1]] %in% t[k,]), dimnames(c_temp)[[2]]])
            temp <- temp / sum(temp)
            r_temp[k, j] <- sum(c_temp[k, ] * temp)
          }
        } else {
          temp <- colSums(region.weights[which(dimnames(region.weights)[[1]] %in% as.character(period)), dimnames(c_temp)[[2]]])
          temp <- temp / sum(temp)
          r_temp[, j] <- c_temp %*% matrix(temp)
        }
      } else {
        r_temp[, j] <- c_temp
      }
    }
    r.temp <- c(r.temp, list(r_temp))
    
    if (m) {
      w.temp[-which(dimnames(w.temp)[[1]] %in% r_names), r_names[i], ] <- rowSums(weight.data[-which(dimnames(weight.data)[[1]] %in% all_r_countries), regions[[i]],])
      w.temp[r_names[i] , -which(dimnames(w.temp)[[2]]%in%r_names), ] <- colSums(weight.data[regions[[i]], -which(dimnames(weight.data)[[2]] %in% all_r_countries),])
    } else {
      w.temp[-which(dimnames(w.temp)[[1]] %in% r_names), r_names[i], ] <- apply(weight.data[-which(dimnames(weight.data)[[1]] %in% all_r_countries), regions[[i]],], 3, rowSums)
      w.temp[r_names[i] , -which(dimnames(w.temp)[[2]]%in%r_names), ] <- apply(weight.data[regions[[i]], -which(dimnames(weight.data)[[2]] %in% all_r_countries), ], 3, colSums)
    }
    for (j in r_names) {
      if (r_names[i] != j) {
        w.temp[r_names[i], j,] <- apply(weight.data[regions[[i]], regions[[j]],], 3, sum)
      } else {
        next
      }
    }
  }
  names(r.temp) <- names(regions)
  
  w.temp[,,] <- apply(w.temp, 3, function(x) {x / rowSums(x)})
  
  data <- c()
  data.names <- c()
  for (i in names(country.data)) {
    if (!i %in% all_r_countries) {
      data <- c(data , list(country.data[[i]]))
      data.names <- c(data.names, i)
    }
  }
  for (i in names(r.temp)) {
    data <- c(data, list(r.temp[[i]]))
    data.names <- c(data.names, i)
  }
  names(data) <- data.names
  
  country.data <- data
  weight.data <- w.temp
  
  if (m) {
    w.names <- dimnames(weight.data)
    weight.data <- matrix(weight.data, dim(weight.data)[1])
    dimnames(weight.data) <- list(w.names[[1]], w.names[[2]])
  }
  
  result <- list("country.data" = country.data,
                 "weight.data" = weight.data)
  return(result)
}