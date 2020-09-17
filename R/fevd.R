#' Forecast Error Variance Decomposition
#' 
#' Produces the forecast error variance decomposition of a Bayesian VAR model.
#' 
#' @param object an object of class \code{"bvar"}, usually, a result of a call to \code{\link{bvar}}
#' or \code{\link{bvec_to_bvar}}.
#' @param response name of the response variable.
#' @param n.ahead number of steps ahead.
#' @param type type of the impulse responses used to calculate forecast error variable decompositions.
#' Possible choices are orthogonalised \code{oir} (default) and generalised \code{gir} impulse responses.
#' @param normalise_gir logical. Should the GIR-based FEVD be normalised?
#' 
#' @details The function produces forecast error variance decompositions (FEVD) for the VAR model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}. For non-structural models matrix \eqn{A_0} is set to the identiy matrix
#' and can therefore be omitted, where not relevant.
#' 
#' If the FEVD is based on the orthogonalised impulse resonse (OIR), the FEVD will be calculated as
#' \deqn{\omega^{OIR}_{jk, h} = \frac{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i P e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma \Phi_i^{\prime} e_j )},}
#' where \eqn{\Phi_i} is the forecast error impulse response for the \eqn{i}th period,
#' \eqn{P} is the lower triangular Choleski decomposition of the variance-covariance
#' matrix \eqn{\Sigma}, \eqn{e_j} is a selection vector for the response variable and
#' \eqn{e_k} a selection vector for the impulse variable.
#'
#' If \code{type = "sir"}, the structural FEVD will be
#' calculated as \deqn{\omega^{SIR}_{jk, h} = \frac{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} A_0^{-1\prime} \Phi_i^{\prime} e_j )},}
#' where \eqn{\sigma_{jj}} is the diagonal element of the \eqn{j}th variable of the variance covariance matrix.
#'
#' If \code{type = "gir"}, the generalised FEVD will be
#' calculated as \deqn{\omega^{GIR}_{jk, h} = \frac{\sigma^{-1}_{jj} \sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma \Phi_i^{\prime} e_j )},}
#' where \eqn{\sigma_{jj}} is the diagonal element of the \eqn{j}th variable of the variance covariance matrix.
#' 
#' If \code{type = "sgir"}, the structural generalised FEVD will be
#' calculated as \deqn{\omega^{SGIR}_{jk, h} = \frac{\sigma^{-1}_{jj} \sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} \Sigma e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} \Sigma A_0^{-1\prime} \Phi_i^{\prime} e_j )}}.
#' 
#' Since GIR-based FEVDs do not add up to unity, they can be normalised by setting \code{normalise_gir = TRUE}.
#' 
#' @return A time-series object of class \code{"bvarfevd"}.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate models
#' model <- gen_var(e1, p = 2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Obtain FEVD
#' vd <- fevd(object, response = "cons")
#' 
#' # Plot FEVD
#' plot(vd)
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' @export
fevd <- function(object, response = NULL, n.ahead = 5, type = "oir", normalise_gir = FALSE) {
  
  # Dev specs
  # rm(list = ls()[-which(ls() == "object")]); response = "invest"; n.ahead = 20; type = "oir"; normalise_gir = FALSE
  
  if (!"bvar" %in% class(object)) {
    stop("Object must be of class 'bvar'.")
  }
  if (is.null(object$y) | is.null(dimnames(object$y)[[2]])) {
    stop("Argument 'object' must include a named matrix of endogenous variables.")
  }
  if (is.null(object$Sigma)) {
    stop("Argument 'object' must include draws of the variance-covariance matrix Sigma.")
  }
  if (!type %in% c("oir", "sir", "gir", "sgir")) {
    stop("The specified type of the used impulse response is not known.")
  }
  if(is.null(response)) {
    stop("Please provide a valid response variable.")
  }
  if (is.null(object[["A"]]) & !type %in% c("sir", "sgir")) {
    stop("Variance decompositions only supported for models with p > 0, i.e. argument 'object' must contain element 'A', or structural models.")
  }
  
  need_A0 <- FALSE
  if (type %in% c("sgir", "sir")) {
    if (is.null(object$A0)) {
      stop("Structural FEVD requires that draws of 'A0' are contained in the 'bvar' object.")
    }
    need_A0 <- TRUE
  }
  
  response <- which(dimnames(object$y)[[2]] == response)
  if (length(response) == 0){stop("Response variable not available.")}
  
  k <- ncol(object$y)
  if (!is.null(object[["A"]])) {
    store <- NROW(object[["A"]]) # Number of draws 
  } else {
    store <- NROW(object[["A0"]])
  }
  
  A <- NULL
  for (i in 1:store) {
    temp <- NULL
    if (!is.null(object[["A"]])) {
      temp[["A"]] <- matrix(object[["A"]][i, ], k)
    } else {
      temp[["A"]] <- matrix(0, k, k)
    }
    temp[["Sigma"]] <- matrix(object[["Sigma"]][i, ], k)
    if (need_A0) {
      temp[["A0"]] <- matrix(object[["A0"]][i, ], k)
    }
    A[[i]] <- temp
  }
  
  phi <- lapply(A, .vardecomp, h = n.ahead, type = type, response = response)
  
  result <- matrix(rowMeans(matrix(unlist(phi), (n.ahead + 1) * k)), n.ahead + 1)
  
  if (type %in% c("gir", "sgir")) {
    if (normalise_gir) {
      result <- t(apply(result, 1, function(x) {x / sum(x)}))
    }
  }
  
  result <- stats::ts(result, start = 0, frequency = 1)
  dimnames(result) <- list(NULL, dimnames(object$y)[[2]]) # Name columns
  
  class(result) <- append("bvarfevd", class(result))
  return(result)
}