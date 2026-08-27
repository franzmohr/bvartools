#' Forecast Error Variance Decomposition
#' 
#' Produces the forecast error variance decomposition for an object of class 'bvarmodel'.
#' 
#' @param x an object of class 'bvarmodel'.
#' @param response name of the response variable.
#' @param n_ahead number of steps ahead.
#' @param type type of the impulse responses used to calculate forecast error variable decompositions.
#' Possible choices are orthogonalised \code{oir} (default) and generalised \code{gir} impulse responses.
#' @param normalise_gir logical. Should the GIR-based FEVD be normalised?
#' @param period integer. Index of the period, for which the variance decomposition should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param ... further arguments passed to or from other methods.
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
#' @return A time-series object of class 'bvarfevd'.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' e1 <- window(e1, end = c(1978, 4))
#' 
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' object <- add_posterior_coefficients(model)
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
fevd.bvarmodel <- function(x, response = NULL, n_ahead = 5, type = "oir", normalise_gir = FALSE, period = NULL, ...) {
  
  
  if (is.null(x[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Argument 'object' must include draws of the variance-covariance matrix Sigma.")
  }
  
  if (!type %in% c("oir", "sir", "gir", "sgir")) {
    stop("The specified type of the used impulse response is not known.")
  }
  
  if(is.null(response)) {
    stop("Please provide a valid response variable.")
  }
  
  if (x[["model"]][["p"]] == 0 & !x[["model"]][["structural"]]) {
    stop("Variance decompositions only supported for models with p > 0 or structural models.")
  }

  need_A0 <- FALSE
  if (type %in% c("sgir", "sir")) {
    if (!x[["model"]][["structural"]]) {
      stop("Structural FEVD requires a structural model as input.")
    }
    need_A0 <- TRUE
  }
  
  varnames <- x[["model"]][["endogen"]]
  response <- which(varnames == response)
  if (length(response) == 0){stop("Response variable not available.")}
  
  k <- x[["model"]][["k"]]
  kk <- k * k
  p <- x[["model"]][["p"]]
  tt <- nrow(x[["data"]][["train"]][["y"]]) / k
  tvp <- x[["model"]][["tvp"]]
  tvp_and_covar <- tvp & x[["model"]][["error"]] %in% c("gamma", "gamma+covar")
  if (tvp) {
    nparams <- ncol(x[["data"]][["train"]][["z"]])
  }
  sv <- x[["model"]][["error"]] %in% c("sv", "sv+covar")
  if (tvp || sv) {
    if (is.null(period)) {
      period <- tt
    } else {
      if (period > tt | period < 1) {
        stop("Implausible specification of argument 'period'.")
      }
    }
  }
  
  if (need_A0) {
    n_struct <- k * (k - 1) / 2
    
    if (tvp) {
      pos_a <- nparams * period - n_struct + 1:n_struct
    } else {
      n_a <- ncol(x[["posterior"]][["a"]][["coeffs"]])
      pos_a <- n_a - n_struct + 1:n_struct 
    }
    
    pos_a0 <- t(matrix(1:kk, k , k))
    pos_a0 <- pos_a0[upper.tri(pos_a0)]
  }
  
  store <- nrow(x[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  
  A <- NULL
  for (i in 1:store) {
    temp <- NULL
    if (p > 0) {
      if (tvp) {
        temp[["A"]] <- matrix(x[["posterior"]][["a"]][["coeffs"]][i, (period - 1) * nparams + 1:(kk * p)], k)
      } else {
        temp[["A"]] <- matrix(x[["posterior"]][["a"]][["coeffs"]][i, 1:(kk * p)], k) 
      }
    } else {
      temp[["A"]] <- matrix(0, k, k)
    }
    
    if (need_A0) {
      a0_temp <- diag(1, k)
      a0_temp[pos_a0] <- x[["posterior"]][["a"]][["coeffs"]][i, pos_a]
      temp[["A0"]] <- a0_temp
      temp[["A"]] <- solve(a0_temp) %*% temp[["A"]]
    }
    
    if (sv | tvp_and_covar) {
      temp[["Sigma"]] <- solve(matrix(x[["posterior"]][["u_sigma_inv"]][["coeffs"]][i, (period - 1) * kk + 1:kk], k))
    } else {
      temp[["Sigma"]] <- solve(matrix(x[["posterior"]][["u_sigma_inv"]][["coeffs"]][i, ], k))
    }
    
    A[[i]] <- temp
  }
  
  phi <- lapply(A, .vardecomp, h = n_ahead, type = type, response = response)
  
  result <- matrix(rowMeans(matrix(unlist(phi), (n_ahead + 1) * k)), n_ahead + 1)
  
  if (type %in% c("gir", "sgir")) {
    if (normalise_gir) {
      result <- t(apply(result, 1, function(x) {x / sum(x)}))
    }
  }
  
  result <- stats::ts(result, start = 0, frequency = 1)
  dimnames(result) <- list(NULL, varnames) # Name columns
  
  class(result) <- append("bvarfevd", class(result))
  return(result)
}