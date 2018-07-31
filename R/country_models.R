#' Country Models
#'
#' Produces a list of country models, which can be estimated with the \code{\link{gvar_fit}} function.
#'
#' @param country.data a named list of "ts" objects containing country data.
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
  # Check if country data is a list
  if (class(country.data) != "list") {
    stop("Country data must have class 'list'.")
  }
  if(sum(unlist(lapply(country.data, class)) == "ts") != length(country.data)) {stop("country.data must be of class 'ts'.")}
  if (!is.null(global.data)) {
    if(!"ts" %in% class(global.data)) {stop("Data must be of class 'ts'.")}
  }
  
  # Check if data is named
  if(is.null(names(country.data))){
    stop("Country data must be named. Please, provide a name for each country.")
  }
  
  if (is.null(dimnames(weight.data)[[1]]) | is.null(dimnames(weight.data)[[2]])) {
    stop("Weight matrix rows and columns must both be named.")
  }
  
  if(is.null(dimnames(global.data)[[2]])){
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
    if (dim(weight.data)[3] > t) {
      warning("The weight array does not contain as much periods as the country sample. Trying to correct this.")
      dim_w <- as.numeric(dimnames(weight.data)[[3]])
      w.t <- stats::ts(dim_w, start = dim_w[1], frequency = frequency(country.data[[1]]))
      weight.data <- weight.data[,,-which(!w.t %in% time(country.data[[1]]))]
    }
  }
  
  # Trim global variables if neccessary
  used.global.vars <- unique(unlist(lapply(country.specs, function(x) {return(x$global.variables)})))
  if (!all(is.na(used.global.vars))){
    global.data <- as.matrix(global.data[, used.global.vars])
    dimnames(global.data)[[2]] <- used.global.vars
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
    names.X <- dimnames(X)[[2]]
    names.g <- dimnames(global.data)[[2]]
    temp <- stats::na.omit(cbind(X, global.data))
    X <- temp[, 1:k]
    X.global <- stats::as.ts(as.matrix(temp[, -(1:k)]))
    stats::tsp(X.global) <- stats::tsp(X)
    dimnames(X)[[2]] <- names.X
    dimnames(X.global)[[2]] <- names.g
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
  
  X.index <- stats::tsp(X)
  data <- c()
  for (i in 1:length(Z)){
    # Split z into domestic and foreign
    x <- stats::ts(t(matrix(Z[[i]][country.specs[[i]]$domestic.variables,], length(country.specs[[i]]$domestic.variables))),
                   start = X.index[1], frequency = X.index[3], class = c("mts", "ts", "matrix"))
    dimnames(x)[[2]] <- country.specs[[i]]$domestic.variables
    
    x.star <- stats::ts(t(matrix(Z[[i]][-(1:length(country.specs[[i]]$domestic.variables)),], length(country.specs[[i]]$foreign.variables))),
                        start = X.index[1], frequency = X.index[3], class = c("mts", "ts", "matrix"))
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
    if (!any(is.na(data[[i]]$specs$rank$rank))) {
      if (data[[i]]$specs$type == "VEC" & any(r_temp == "rr")) {
        data[[i]]$specs$rank$rank <- 0:length(data[[i]]$specs$domestic.variables)
      } 
    }
    for (j in 1:length(data[[i]]$specs$lags$lags$domestic)) {
      for (k in 1:length(data[[i]]$specs$lags$lags$foreign)) {
        for (l in 1:length(data[[i]]$specs$lags$lags$global)) {
          for (m in 1:length(data[[i]]$specs$rank$rank)) { # For multiple rank specifications of country VEC models.
            temp <- data[[i]]
            temp$specs$lags$lags$domestic <- temp$specs$lags$lags$domestic[j]
            temp$specs$lags$lags$foreign <- temp$specs$lags$lags$foreign[k]
            temp$specs$lags$lags$global <- temp$specs$lags$lags$global[l]
            temp$specs$rank$rank <- temp$specs$rank$rank[m]
            data.temp <- c(data.temp, list(temp))
            names.temp <- c(names.temp, names(data)[i])  
          } 
        }
      }
    }
  }
  names(data.temp) <- names.temp
  data <- data.temp
  
  #### Add priors to specifications ####
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
      if (det == "IV" || det == "V" || det == "VI"){
        n.det <- 2
      }
      p <- temp$specs$lags$lags$domestic
      p.star <- temp$specs$lags$lags$foreign
      p.global <- temp$specs$lags$lags$global
    }
    if (temp.type == "VEC"){
      if (det == "III" || det == "IV"){
        n.det <- 1
      }
      if (det == "V"){
        n.det <- 2
      } 
      p <- temp$specs$lags$lags$domestic - 1
      p.star <- temp$specs$lags$lags$foreign - 1
      p.global <- temp$specs$lags$lags$global - 1
    }
    if (is.na(p.global)){
      p.global <- 0
    }
    
    r <- temp$specs$rank$rank
    n.ect <- n * (n + n.s + n.g)
    if (det == "II" || det == "IV"){
      n.ect <- n.ect + n
    }
    if (det == "VI"){
      n.ect <- n.ect + 2 * n
    }
    
    if (!is.na(r)) {
      n.alpha <- n * r 
    }
    n.A <- n * (n * p + n.s * (1 + p.star) + n.g * (1 + p.global) + n.det)
    n.A0 <- n * (n - 1) / 2
    
    if (is.null(prior)){
      prior <- standard_priors()
    }
    
    A_mu_prior <- matrix(prior$A$constant[1], n.A)
    A_V_i_prior <- diag(prior$A$constant[2], n.A)
    A_Q_df_prior <- matrix(prior$A$tvp[1], n.A)
    A_Q_V_prior <- diag(prior$A$tvp[2], n.A)
    
    Pi_tvp <- temp$specs$tvp[2]
    if (is.na(r)) {
      Pi_mu_prior <- matrix(prior$Pi$non_rr$constant[1], n.ect)
      Pi_V_i_prior <- diag(prior$Pi$non_rr$constant[2], n.ect)
      if (det == "II" || det == "IV"){
        Pi_mu_prior[(n.ect - n + 1):n.ect] <- prior$Deterministic$constant[1]
        diag(Pi_V_i_prior)[(n.ect - n + 1):n.ect]  <- prior$Deterministic$constant[2]
      }
      if (det == "VI"){
        Pi_mu_prior[(n.ect - 2 * n + 1):n.ect] <- prior$Deterministic$constant[1]
        diag(Pi_V_i_prior)[(n.ect - 2 * n + 1):n.ect]  <- prior$Deterministic$constant[2]
      }
      Pi_Q_df_prior <- matrix(prior$Pi$non_rr$tvp[1], n.ect)
      Pi_Q_V_prior <- diag(prior$Pi$non_rr$tvp[2], n.ect)
      rho <- NA
    } else {
      Pi_mu_prior <- matrix(prior$Pi$reduced_rank$constant[[1]], n.alpha)
      Pi_V_i_prior <- list("v_i" = prior$Pi$reduced_rank$constant[[2]],
                           "P_i" = diag(prior$Pi$reduced_rank$constant[[3]], n.ect / n),
                           "G_i" = NA)
      if (prior$Pi$reduced_rank$constant[[4]] == "Omega_i") {
        Pi_V_i_prior$G_i <- prior$Pi$reduced_rank$constant[[4]]
      } else {
        Pi_V_i_prior$G_i <- diag(prior$Pi$reduced_rank$constant[[4]], n) 
      }
      Pi_Q_df_prior <- matrix(prior$Pi$reduced_rank$tvp$alpha[1], n.alpha)
      Pi_Q_V_prior <- diag(prior$Pi$reduced_rank$tvp$alpha[2], n.alpha)
      rho <- prior$Pi$reduced_rank$tvp$rho
    }
    
    if (n.det > 0){
      pos.det <- n.A - (n.det * n) + 1:(n.det * n)
      A_mu_prior[pos.det, 1] <- prior$Deterministic$constant[1]
      diag(A_V_i_prior)[pos.det] <- prior$Deterministic$constant[2]
      A_Q_df_prior[pos.det, 1] <- prior$Deterministic$tvp[1]
      diag(A_Q_V_prior)[pos.det] <- prior$Deterministic$tvp[2]
    }
    
    if (n.A0 > 0) {
      A0_mu_prior <- matrix(prior$A0$constant[1], n.A0)
      A0_V_i_prior <- diag(prior$A0$constant[2], n.A0)
      A0_Q_df_prior <- matrix(prior$A0$tvp[1], n.A0)
      A0_Q_V_prior <- diag(prior$A0$tvp[2], n.A0) 
    } else {
      A0_mu_prior <- NULL
      A0_V_i_prior <- NULL
      A0_Q_df_prior <- NULL
      A0_Q_V_prior <- NULL
    }
    
    Omega_df_prior <- matrix(prior$Omega$constant[1], n)
    Omega_V_i_prior <- diag(prior$Omega$constant[2], n)
    
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
          A.restricted.variables <- matrix(1:(n.A - n.det * n), n.A - n.det * n)
        } else {
          A.restricted.variables <- matrix(1:n.A, n.A) 
        }
        A.lpr.include.prior <- log(matrix(prior$Shrinkage$spec, n.A))
        A.lpr.exclude.prior <- log(matrix(1 - prior$Shrinkage$spec, n.A))
        
        if (n.A0 > 0) {
          A0.restricted.variables <- matrix(1:n.A0, n.A0)
          A0.lpr.include.prior <- log(matrix(prior$Shrinkage$spec, n.A0))
          A0.lpr.exclude.prior <- log(matrix(1 - prior$Shrinkage$spec, n.A0)) 
        } else {
          A0.restricted.variables <- NULL
          A0.lpr.include.prior <- NULL
          A0.lpr.exclude.prior <- NULL
        }
        
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
    
    data[[i]]$priors <- list(A = list("constant" = list(A_mu_prior, A_V_i_prior),
                                      "tvp" = list(A_Q_df_prior, A_Q_V_prior)),
                             A0 = list("constant" = list(A0_mu_prior, A0_V_i_prior),
                                       "tvp" = list(A0_Q_df_prior, A0_Q_V_prior)),
                             Omega = list("constant" = list(Omega_df_prior, Omega_V_i_prior)),
                             Sigma = list("sv" = log.V),
                             Shrinkage = shrinkage)
    
    if (temp.type == "VEC"){
      names(data[[i]]$priors)[1] <- "B"
      if (is.na(r)) {
        data[[i]]$priors$Pi <- list("constant" =  list("mu" = Pi_mu_prior, "V_i" = Pi_V_i_prior),
                                    "tvp" = list(Pi_Q_df_prior, Pi_Q_V_prior, rho))
      }
    }
  }
  data <- list(country.data = data, global.data = list(X = X, X.global = X.global, W = W, index = index))
  return(data)
}