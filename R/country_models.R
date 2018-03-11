#' Country Models
#'
#' Produces a list of country models, which can be estimated with the \code{\link{gvar_fit}} function.
#'
#' @param country.data a named list of "zoo" objects containing country data.
#' @param weight.data a matrix or an array of weights used to generate the weight matrices.
#' @param global.data a "zoo" object with global data.
#' @param country.specs a list of model specifications as produced by \code{\link{country_specifications}}.
#' @param prior a list of priors as produced by \code{\link{standard_priors}}.
#' 
#' @return A list containing the following components:
#' \item{Country.data}{a list of data, model specifications and priors for each estimated country model.}
#' \item{Global.data}{a list containing the data of the global model, the country-specific link matrices and a variable index.}
#' 
#' @export
country_models <- function(country.data, weight.data, global.data = NULL,
                           country.specs = NULL, prior = NULL){
  
  if(!requireNamespace("zoo")) {stop("Function 'country_models' requires the package 'zoo'.")}
  if(any(unlist(lapply(country.data, class)) != "zoo")) {stop("Country data must be of class 'zoo'.")}
  if (!is.null(global.data)) {
    if(class(global.data) != "zoo") {stop("Global data must be of class 'zoo'.")}
  }
  
  # Check if country data is a list
  if (class(country.data) != "list") {
    stop("Country data must have class 'list'.")
  }
  
  # Check if time series data has class zoo
  if(!all(unlist(lapply(country.data, class)) == "zoo")) {
    stop("Country time series data must be class 'zoo'.")
  }
  if (!is.null(global.data)){
    if(!all(unlist(lapply(global.data, class)) == "zoo")) {
      stop("Global time series data must be class 'zoo'.")
    }
  }
  
  # Check if data is named
  if(is.null(names(country.data))){
    stop("Country data must be named. Please, provide a name for each country.")
  }
  
  if (is.null(dimnames(weight.data)[[1]]) | is.null(dimnames(weight.data)[[2]])) {
    stop("Weight matrix rows and columns must both be named.")
  }
  
  if(is.null(names(global.data))){
    stop("Global data must be named. Please, provide a name for each global variable.")
  }
  
  const.weight <- class(weight.data) == "matrix"
  
  # Limit country data and weight matrix to the selected countries
  countries <- names(country.specs)
  temp.spec <- NULL
  temp.c <- NULL
  if (const.weight) {
    temp.w <- matrix(0, length(countries), length(countries))
    dimnames(temp.w) <- list(countries, countries)
  } else {
    temp.w <- array(0, dim = c(length(countries), length(countries), dim(weight.data)[[3]]))
    dimnames(temp.w) <- list(countries, countries, dimnames(weight.data)[[3]])
  }

  for (i in countries){
    temp.spec <- c(temp.spec, list(country.specs[[i]]))
    temp.c <- c(temp.c, list(country.data[[i]]))
    if (const.weight){
      temp.w[i,] <- weight.data[i, countries]/(1 - sum(weight.data[i, !is.element(dimnames(weight.data)[[2]], countries)]))
    } else {
      for (j in 1:dim(weight.data)[3]){
        temp.w[i,,j] <- weight.data[i, countries, j]/(1 - sum(weight.data[i, !is.element(dimnames(weight.data)[[2]], countries), j])) 
      }
    }
  }
  country.specs <- temp.spec; rm(temp.spec)
  country.data <- temp.c; rm(temp.c)
  names(country.specs) <- countries
  names(country.data) <- countries
  weight.data <- temp.w; rm(temp.w)
  
  # Check if all country time series have the same length
  n.c <- length(country.data)
  t.max <- c()
  for (i in 1:n.c){
    if (any(apply(country.data[[i]], 2, function(x){any(is.na(x))}))){
      vars <- dimnames(country.data[[i]])[[2]][which(apply(country.data[[i]], 2, function(x){any(is.na(x))}))]
      country <- names(country.data)[i]
      stop(paste("In country ", country, " the variable(s) ", paste(vars, sep = ", "), " contain NA values.", sep = ""))
    }
    t.max <- append(t.max, dim(stats::na.omit(country.data[[i]]))[1])
  }
  rm(n.c)
  
  if (min(t.max) != max(t.max)) {
    stop("Numbers of observations differ across countries. For now, this package requires the same number for each country.")
  } else {
    t <- t.max[1]
    rm(t.max)
  }
  
  # Check if weight matrix has the same order of names in the rows and columns
  if (any(dimnames(weight.data)[[1]] != dimnames(weight.data)[[2]])){
    stop("The order of country names in the rows and columns of the weight matrix is not the same.")
  }
  
  # Check if weight matrix contains data on all required countries
  if (any(!is.element(names(country.data), dimnames(weight.data)[[1]]))){
    stop("The weight matrix does not contain data for at least one country in the country sample or is named differently.")
  }
  
  # If weight matrix is time varying, check if the number of periods is the same as in the country series
  if (!is.na(dim(weight.data)[3])){
    if (dim(weight.data)[3] != t) {
      warning("The weight array does not contain as much periods as the country sample. Trying to correct that.")
      w.t <- zoo::as.yearqtr(ts(as.numeric(dimnames(weight.data)[[3]]),
                                start = as.numeric(dimnames(weight.data)[[3]])[1],
                                frequency = frequency(country.data[[1]])))
      weight.data <- weight.data[,,-which(!w.t%in%zoo::index(country.data[[1]]))]
    }
  }

  # Trim global variables if neccessary
  used.global.vars <- unique(unlist(lapply(country.specs, function(x) {return(x$global.variables)})))
  if (!all(is.na(used.global.vars))){
    global.data <- global.data[, used.global.vars]
  }
  
  # Get endogenous variables
  endo.data <- global_series(country.data, country.specs, global.data)
  X <- endo.data$endogenous.variables # Data for global model (without deterministic terms and global data)
  global.data <- endo.data$global.data # Updated global variables
  index <- endo.data$index # Index with country and variable names
  country.specs <- endo.data$country.specs # Updated specifications of country models
  rm(endo.data)
  
  # Make global and country data have the same length
  global <- !is.null(global.data) # If a global variable is used
  if (global){
    k <- dim(X)[2]
    names.X <- names(X)
    names.g <- names(global.data)
    temp <- stats::na.omit(cbind(X, global.data))
    X <- temp[, 1:k]
    X.global <- temp[, -(1:k)]
    names(X) <- names.X
    names(X.global) <- names.g
    rm(list=c("temp", "names.X", "names.g"))
  } else {
    X.global <- NULL
  }
  rm(global.data)
  
  # Add the number of the maximum lag order of the country models to the country specifications.
  # This is necessary for the calculation of comparable information criteria for model selection.
  lag.max <- max(unlist(lapply(country.specs, function(x) {return(max(unlist(x$lags$lags), na.rm = TRUE))})))
  for (i in 1:length(country.specs)){
    country.specs[[i]]$lags <- c(country.specs[[i]]$lags, maximum.lag = lag.max)
  }
  
  # Generate weight matrices for each country
  W <- weight_matrices(weight.data, index, country.specs)
  
  # Generate vector z = (x, x.star)' for each country
  Z <- lapply(W, get_z, X)
  
  X.index <- zoo::index(X)
  data <- c()
  for (i in 1:length(Z)){
    # Split z into domestic and foreign
    x <- zoo::zoo(t(matrix(Z[[i]][country.specs[[i]]$domestic.variables,], length(country.specs[[i]]$domestic.variables))), order.by = X.index)
    names(x) <- country.specs[[i]]$domestic.variables
    x.star <- zoo::zoo(t(Z[[i]][-(1:length(country.specs[[i]]$domestic.variables)),]), order.by = X.index)
    # Check for global variables
    if (!any(is.na(country.specs[[i]][[3]]))){
      g <- X.global
    } else {
      g <- NULL
    }
    # Collect data
    data <- c(data, list(list(x = x, x.star = x.star, x.g = g, weights = W[[i]], specs = country.specs[[i]])))
    rm(list=c("x","x.star","g"))
  }
  names(data) <- names(Z)
  
  # Generate country model for each lag specification
  data.temp <- NULL
  names.temp <- NULL
  for (i in 1:length(data)){
    for (j in 1:length(data[[i]]$specs$lags$lags$domestic)){
      for (k in 1:length(data[[i]]$specs$lags$lags$foreign)){
        for (l in 1:length(data[[i]]$specs$rank$rank)){ # For multiple rank specifications of country VEC models.
          temp <- data[[i]]
          temp$specs$lags$lags$domestic <- temp$specs$lags$lags$domestic[j]
          temp$specs$lags$lags$foreign <- temp$specs$lags$lags$foreign[k]
          temp$specs$rank$rank <- temp$specs$rank$rank[l]
          data.temp <- c(data.temp, list(temp))
          names.temp <- c(names.temp, names(data)[i])  
        }
      }
    }
  }
  names(data.temp) <- names.temp
  data <- data.temp
  
  # Add priors to specifications
  for (i in 1:length(data)){
    temp <- data[[i]]
    temp.type <- temp$specs$type
    n <- length(temp$specs$domestic.variables)
    n.s <- length(temp$specs$foreign.variables)
    if (is.na(temp$specs$global.variables)){
      n.g <- 0
    } else {
      n.g <- length(temp$specs$global.variables)
    }
    det <- temp$specs$deterministic.terms
    n.det <- 0
    if (temp.type == "VAR"){
      if (det == "II" || det == "III"){
        n.det <- 1
      }
      if (det == "IV" || det == "V"){
        n.det <- 2
      }
      p <- temp$specs$lags$lags$domestic
      q <- temp$specs$lags$lags$foreign
      s <- temp$specs$lags$lags$global
    }
    if (temp.type == "VEC"){
      if (det == "III" || det == "IV"){
        n.det <- 1
      }
      if (det == "V"){
        n.det <- 2
      } 
      p <- temp$specs$lags$lags$domestic - 1
      q <- temp$specs$lags$lags$foreign - 1
      s <- temp$specs$lags$lags$global - 1
    }
    if (is.na(s)){
      s <- 0
    }
    
    r <- temp$specs$rank$rank
    n.A <- n*(r + n*p + n.s*(1 + q) + n.g*(1 + s) + n.det)
    n.A0 <- n*(n-1)/2
    
    if (is.null(prior)){
      prior <- standard_priors()
    }
    
    A.mu.prior <- matrix(prior$A$constant[1], n.A)
    A.V.prior.i <- diag(prior$A$constant[2], n.A)
    A.H.df.prior <- prior$A$tvp[1]
    A.H.V.prior <- diag(prior$A$tvp[2],n.A)
    if (r > 0){
      pos.r <- 1:(n*r)
      A.mu.prior[pos.r, 1] <- prior$Pi$alpha$constant[1]
      diag(A.V.prior.i)[pos.r] <- prior$Pi$alpha$constant[2]
      diag(A.H.V.prior)[pos.r] <- prior$Pi$alpha$tvp[2]
    }
    if (n.det > 0){
      pos.det <- n.A - (n.det*n) + 1:(n.det*n)
      A.mu.prior[pos.det, 1] <- prior$Deterministic$constant[1]
      diag(A.V.prior.i)[pos.det] <- prior$Deterministic$constant[2]
    }
    
    A0.mu.prior <- matrix(prior$A0$constant[1], n.A0)
    A0.V.prior.i <- diag(prior$A0$constant[2], n.A0)
    A0.H.df.prior <- prior$A0$tvp[1]
    A0.H.V.prior <- diag(prior$A0$tvp[2], n.A0)
    
    Sigma.df.prior <- prior$Sigma$constant[1]
    Sigma.V.prior <- diag(prior$Sigma$constant[2], n)
    
    log.V <- NULL
    for (j in 1:n){
      log.V <- c(log.V, list(list(para = list(mu = prior$Sigma$sv$para$mu,
                                              phi = prior$Sigma$sv$para$phi,
                                              sigma = prior$Sigma$sv$para$sigma),
                                  latent = prior$Sigma$sv$latent,
                                  priors = list(priormu = prior$Sigma$sv$priors$priormu,
                                                priorphi = prior$Sigma$sv$priors$priorphi,
                                                priorsigma = prior$Sigma$sv$priors$priorsigma))))
    }
    
    if (prior$Shrinkage$type != "none"){
      shrinkage.type <- prior$Shrinkage$type
      if (shrinkage.type == "BVS"){
        if (prior$A$constant[2] == 0 | prior$A0$constant[2] == 0) {
          stop("BVS requires a mildly informative prior.")
        }
        if (prior$Shrinkage$exclude.deterministic){
          A.restricted.variables <- matrix(1:(n.A - n.det*n), n.A - n.det*n)
        } else {
          A.restricted.variables <- matrix(1:n.A, n.A) 
        }
        if (r > 0){
          A.restricted.variables <- matrix(A.restricted.variables[-(1:(n*r)),])
        }
        A0.restricted.variables <- matrix(1:n.A0, n.A0)
        A.lpr.include.prior <- log(matrix(prior$Shrinkage$spec, n.A))
        A.lpr.exclude.prior <- log(matrix(1 - prior$Shrinkage$spec, n.A))
        A0.lpr.include.prior <- log(matrix(prior$Shrinkage$spec, n.A0))
        A0.lpr.exclude.prior <- log(matrix(1 - prior$Shrinkage$spec, n.A0))
        shrinkage <- list("type" = shrinkage.type,
                          "spec" = list("A" = list(A.restricted.variables, A.lpr.include.prior, A.lpr.exclude.prior),
                                        "A0" = list(A0.restricted.variables, A0.lpr.include.prior, A0.lpr.exclude.prior)))
        if (temp.type == "VEC"){
          names(shrinkage$spec)[1] <- "B"
        }
      }
    } else {
      shrinkage <- list("type" = "none")
    }
    
    data[[i]]$priors <- list(A = list("constant" = list(A.mu.prior, A.V.prior.i),
                                      "tvp" = list(A.H.df.prior, A.H.V.prior)),
                             A0 = list("constant" = list(A0.mu.prior, A0.V.prior.i),
                                       "tvp" = list(A0.H.df.prior, A0.H.V.prior)),
                             Pi = list("rho" = prior$Pi$rho),
                             Sigma = list("constant" = list(Sigma.df.prior, Sigma.V.prior),
                                          "sv" = log.V),
                             Shrinkage = shrinkage) 
    if (temp.type == "VEC"){
      names(data[[i]]$priors)[1] <- "B"
    }
  }
  data <- list(Country.data = data, Global.data = list(X = X, X.global = X.global, W = W, index = index))
  return(data)
}