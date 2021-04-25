#' Dynamic Factor Model Input
#' 
#' \code{gen_dfm} produces the input for the estimation of a dynamic factor model (DFM).
#' 
#' @param x a time-series object of stationary endogenous variables.
#' @param p an integer vector of the lag order of the measurement equation. See 'Details'.
#' @param n an integer vector of the number of factors. See 'Details'.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 5000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' 
#' @details The function produces the variable matrices of dynamic factor
#' models (DFM) with measurement equation
#' \deqn{x_t = \lambda f_t + u_t,}
#' where
#' \eqn{x_t} is an \eqn{M \times 1} vector of observed variables,
#' \eqn{f_t} is an \eqn{N \times 1} vector of unobserved factors and
#' \eqn{\lambda} is the corresponding \eqn{M \times N} matrix of factor loadings.
#' \eqn{u_t} is an \eqn{M \times 1} error term.
#' 
#' The transition equation is
#' \deqn{f_t = \sum_{i=1}^{p} A_i f_{t - i} + v_t,}
#' where
#' \eqn{A_i} is an \eqn{N \times N} coefficient matrix and
#' \eqn{v_t} is an \eqn{N \times 1} error term.
#' 
#' If integer vectors are provided as arguments \code{p} or \code{n}, the function will
#' produce a distinct model for all possible combinations of those specifications.
#' 
#' @return An object of class \code{'dfmodel'}, which contains the following elements:
#' \item{data}{A list of data objects, which can be used for posterior simulation. Element
#' \code{X} is a time-series object of normalised observable variables, i.e. each column has
#' zero mean and unity variance.}
#' \item{model}{A list of model specifications.}
#' 
#' @examples 
#' 
#' # Load data
#' data("bem_dfmdata")
#' 
#' # Generate model data
#' model <- gen_dfm(x = bem_dfmdata, p = 1, n = 1,
#'                  iterations = 5000, burnin = 1000)
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). \emph{Bayesian Econometric Methods}
#' (2nd ed.). Cambridge: University Press.
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
gen_dfm <- function(x, p = 2, n = 1, iterations = 50000, burnin = 5000) {
  
  # Check data ----
  if (!"ts" %in% class(x)) {
    stop("Argument 'data' must be an object of class 'ts'.")
  }
  
  if (any(p < 0)) {
    stop("Argument 'p' must be at least 0.")
  }
  
  if (any(n < 1)) {
    stop("Argument 'n' must be at least 1.")
  }
  
  if (is.null(dimnames(x))) {
    tsp_temp <- stats::tsp(x)
    data <- stats::ts(as.matrix(x), class = c("mts", "ts", "matrix"))
    stats::tsp(x) <- tsp_temp
    dimnames(x)[[2]] <- "y"
  }
  
  # Normalise every column of x
  x <- scale(x)
  
  data_name <- dimnames(x)[[2]]
  m <- NCOL(x)
  p_max <- max(p)
  
  model <- NULL
  model$type <- "DFM"
  model$factors <- list("variables" = dimnames(x)[[2]],
                        "number" = 1,
                        "lags" = 1)
  
  tt <- nrow(x)
  
  model$iterations <- iterations
  model$burnin <- burnin
  
  result <- NULL
  for (j in n) {
    for (i in p) {
    
      model_i <- model
      model_i$factors$lags <- i
      model_i$factors$number <- j
      
      result_i <- list("data" = list("X" = x),
                       "model" = model_i)
      
      result <- c(result, list(result_i)) 
    }
  }
  
  if (length(result) == 1) {
    result <- result[[1]]
  }
  
  class(result) <- append("dfmodel", class(result))
  return(result)
}