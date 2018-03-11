#' Create Region Series
#' 
#' Combines country-specific time series to regional series and updates the country and weight data.
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
#' regions <- list("EuroArea" = c("Austria", "Belgium", "Finland", "France", "Germany", "Italy", "Netherlands", "Spain"))
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
  
  if ((class(regions) != "list") | is.null(names(regions))) {stop("Object 'regions' must be a named list.")}
  
  if (length(unique(unlist(regions))) < length(unlist(regions))) {
    stop("Different regions may not contain the same country.") 
  }
  
  vars <- unique(unlist(lapply(country.data, names)))
  weights <- c()
  for (i in 1:length(regions)){
    s <- colSums(region.weights[period, regions[[i]]])
    share <- s / sum(s)
    weights <- c(weights, list(share))
    names(weights)[i] <- names(regions)[i]
  }
  r_names <- names(regions)
  all_r_countries <- unlist(lapply(weights, function(x){return(names(x))}))
  names(all_r_countries) <- NULL
  
  n.new <- dim(weight.data)[1] - length(all_r_countries) + length(weights)
  w.temp <- array(0, dim = c(n.new, n.new, dim(weight.data)[3]))
  w.names <- dimnames(weight.data)[[1]][-which(dimnames(weight.data)[[2]]%in%all_r_countries)]
  w.temp[1:length(w.names), 1:length(w.names) , ] <- weight.data[w.names, w.names, ]
  w.names <- c(w.names, names(weights))
  dimnames(w.temp) <- list(w.names, w.names, dimnames(weight.data)[[3]])
  
  r.temp <- c()
  r.names <- c()
  for (i in r_names){
    c.temp <- zoo::zoo(NA, order.by = zoo::index(country.data[[1]]))
    c.names <- c()
    for (j in vars) {
      v.temp <- zoo::zoo(NA, order.by = zoo::index(country.data[[1]]))
      v.names <- c()
      for (k in regions[[i]]) {
        if (j %in% names(country.data[[k]])) {
          v.temp <- cbind(v.temp, country.data[[k]][, j])
          v.names <- c(v.names, k)
        }
      }
      v.temp <- v.temp[, -1]
      names(v.temp) <- v.names
      non.w <- 1 - sum(weights[[i]][v.names])
      if (NCOL(v.temp) > 1){
        c.temp <- cbind(c.temp, v.temp %*% weights[[i]][v.names]/(1-non.w))
      } else {
        c.temp <- cbind(c.temp, v.temp * weights[[i]][v.names]/(1-non.w))
      }
      c.names <- c(c.names, j)
    }
    c.temp <- c.temp[, -1]
    names(c.temp) <- c.names
    r.temp <- c(r.temp, list(c.temp))
    r.names <- c(r.names, i)
    
    w.temp[-which(dimnames(w.temp)[[1]]%in%r_names), i, ] <- apply(weight.data[-which(dimnames(weight.data)[[1]]%in%all_r_countries), regions[[i]],], 3, rowSums)
    w.temp[i , -which(dimnames(w.temp)[[2]]%in%r_names), ] <- apply(weight.data[regions[[i]], -which(dimnames(weight.data)[[2]]%in%all_r_countries), ], 3, colSums)
    
    for (j in r_names) {
      if (i != j) {
        w.temp[i, j,] <- apply(weight.data[regions[[i]], regions[[j]],], 3, sum)
      } else {
        next
      }
    }
  }
  names(r.temp) <- r.names
  
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
