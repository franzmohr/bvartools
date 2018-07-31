#' Priors for Bayesian Variable Selection
#' 
#' Adds prior hyperparameters for Bayesian variable selection to every country model.
#' 
#' @param data a list of country models as produced by \code{\link{country_models}}.
#' @param kappa_domestic number between 0 and 1.
#' @param kappa_foreign number between 0 and 1.
#' @param kappa_global number between 0 and 1.
#' @param kappa_A0 number between 0 and 1.
#' @param kappa_deterministic number between 0 and 1.
#' @param exclude_deterministic logical; exclude deterministic terms from the BVS algorithm?
#' @param exclude_own_lags logical; exclude the first own lags from the BVS algorithm?
#' @param structural logical; should BVS should also be applied to structual coefficients?
#' 
#' @export
add_BVS_prior <- function(data, kappa_domestic = .9, kappa_foreign = .5, kappa_global = .5, kappa_A0 = .5, kappa_deterministic = .8,
                          exclude_deterministic = TRUE, exclude_own_lags = FALSE, structural = FALSE){
  for (i in 1:length(data$country.data)){
    temp <- data$country.data[[i]]
    
    if (any(diag(temp$priors$A$constant[[2]])==0) | any(diag(temp$priors$A0$constant[[2]])==0)) {
      stop("BVS requires at least a mildly informative prior.")
    }
    
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
    n.A <- n * (n * p + n.s * (1 + p.star) + n.g * (1 + p.global) + n.det)
    
    if (exclude_deterministic){
      A.restricted.variables <- matrix(1:(n.A - n.det * n), n.A - n.det * n)
    } else {
      A.restricted.variables <- matrix(1:n.A, n.A)
    }
    
    if (exclude_own_lags) {
      if (p > 0) {
        A.restricted.variables <- matrix(A.restricted.variables[-diag(matrix(A.restricted.variables, n)[,1:n]),]) 
      }
    }
    
    if (p > 0) {
      A.d.prior <- matrix(NA, n, n * p)
      for (j in 1:p) {
        A.d.prior[, (j - 1) * n + 1:n] <- kappa_domestic / (1 + j)
        diag(A.d.prior[, (j - 1) * n + 1:n]) <- kappa_domestic / j
      } 
    }
    
    x_names <- dimnames(temp$x)[[2]]
    x_star_names <- gsub("s.", "", dimnames(temp$x.star)[[2]])
    A.s.prior <- matrix(NA, n, n.s * (p.star + 1))
    dimnames(A.s.prior) <- list(x_names, rep(x_star_names, p.star + 1))
    for (j in 0:p.star) {
      A.s.prior[,j * n.s + 1:n.s] <- kappa_foreign / (1 + j)
    }
    
    if (p > 0) {
      A.prior <- cbind(A.d.prior, A.s.prior) 
    } else {
      A.prior <- A.s.prior
    }
    
    if (n.g > 0) {
      A.g.prior <- matrix(NA, n, n.g * (p.global + 1))
      for (j in 0:p.global) {
          A.g.prior[, j * n.g + 1:n.g] <- kappa_global / (1 + j)
      }
      A.prior <- cbind(A.prior, A.g.prior)
    }
    
    if (det %in% c("II", "III", "IV", "V")) {
      A.det.prior <- matrix(kappa_deterministic, n, n.det)
      A.prior <- cbind(A.prior, A.det.prior)
    }
    
    A.lpr.include.prior <- log(matrix(A.prior))
    A.lpr.exclude.prior <- log(matrix(1 - A.prior))
    
    shrinkage <- list("type" = "BVS",
                      "spec" = list("A" = list(A.restricted.variables, A.lpr.include.prior, A.lpr.exclude.prior)))
    
    if (structural & n > 1) {
      n_A0 <- n * (n - 1) / 2
      A0_restricted_variables <- matrix(1:n.A0, n_A0)
      A0_lpr_include_prior <- log(matrix(kappa_A0, n_A0))
      A0_lpr_exclude_prior <- log(matrix(1 - kappa_A0, n_A0))
      shrinkage$spec$A0 <- list(A0_restricted_variables, A0_lpr_include_prior, A0_lpr_exclude_prior)
    } else {
      shrinkage$spec$A0 <- list(NULL)
    }

    if (temp.type == "VEC"){
      names(shrinkage$spec)[1] <- "B"
    }
    
    temp$priors$Shrinkage <- shrinkage
    data$country.data[[i]] <- temp
  }
  return(data)
}