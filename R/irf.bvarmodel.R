#' Impulse Response Function
#' 
#' Computes the impulse response coefficients for an object of class 'bvarmodel'.
#' 
#' @param x an object of class 'bvarmodel'.
#' @param impulse name of the impulse variable.
#' @param response name of the response variable.
#' @param n_ahead number of steps ahead.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param shock size of the shock.
#' @param type type of the impulse response. Possible choices are forecast error \code{"feir"}
#' (default), orthogonalised \code{"oir"}, structural \code{"sir"}, generalised \code{"gir"},
#' and structural generalised \code{"sgir"} impulse responses.
#' @param cumulative logical specifying whether a cumulative IRF should be calculated.
#' @param keep_draws logical specifying whether the function should return all draws of
#' the posterior impulse response function. Defaults to \code{FALSE} so that
#' the median and the credible intervals of the posterior draws are returned.
#' @param period integer. Index of the period, for which the IR should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The function produces different types of impulse responses for the VAR model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}.
#' 
#' Forecast error impulse responses \eqn{\Phi_i} are obtained by recursions
#' \deqn{\Phi_i = \sum_{j = 1}^{i} \Phi_{i-j} A_j,   i = 1, 2,...,h}
#' with \eqn{\Phi_0 = I_K}.
#' 
#' Orthogonalised impulse responses \eqn{\Theta^o_i} are calculated as \eqn{\Theta^o_i = \Phi_i P},
#' where P is the lower triangular Choleski decomposition of \eqn{\Sigma}.
#' 
#' Structural impulse responses \eqn{\Theta^s_i} are calculated as \eqn{\Theta^s_i = \Phi_i A_0^{-1}}.
#' 
#' (Structural) Generalised impulse responses for variable \eqn{j}, i.e. \eqn{\Theta^g_ji} are calculated as
#' \eqn{\Theta^g_{ji} = \sigma_{jj}^{-1/2} \Phi_i A_0^{-1} \Sigma e_j}, where \eqn{\sigma_{jj}} is the variance
#' of the \eqn{j^{th}} diagonal element of \eqn{\Sigma} and \eqn{e_i} is a selection vector containing
#' one in its \eqn{j^{th}} element and zero otherwise. If the \code{"bvarmodel"} object does not contain draws
#' of \eqn{A_0}, it is assumed to be an identity matrix.
#' 
#' @return A time-series object of class 'bvarirf' or, if \code{keep_draws = TRUE}, a simple matrix.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' model <- add_posterior_coefficients(model)
#' 
#' # Obtain IR
#' ir <- irf(model, impulse = "invest", response = "cons")
#' 
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Pesaran, H. H., Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' @export
#' @method irf bvarmodel
irf.bvarmodel <- function(x, impulse = NULL, response = NULL, n_ahead = 5, ci = .95, shock = 1,
                          type = "feir", cumulative = FALSE, keep_draws = FALSE, period = NULL, ...) {
  
  if (!type %in% c("feir", "oir", "gir", "sir", "sgir")) {
    stop("Argument 'type' not known.")
  }
  
  if (x[["model"]][["p"]] == 0 & !x[["model"]][["structural"]]) {
    stop("Impulse responses only supported for models with p > 0 or structural models.")
  }
  
  need_A0 <- FALSE
  if (type %in% c("sgir", "sir")) {
    if (!x[["model"]][["structural"]]) {
      stop("Structural IR requires a structural model as input.")
    }
    need_A0 <- TRUE
  }
  
  if (!(is.numeric(shock) | shock %in% c("sd", "nsd"))) {
    stop("Invalid specification of argument 'shock'.")
  }
  
  if (type %in% c("oir", "gir", "sgir") | shock %in% c("sd", "nsd")) {
    need_Sigma <- TRUE
  } else {
    need_Sigma <- FALSE
  }
  
  varnames <- x[["model"]][["endogen"]]
  impulse <- which(varnames == impulse)
  if (length(impulse) == 0){stop("Impulse variable not available.")}
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
  if (tvp | sv) {
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
      #temp[["A"]] <- solve(a0_temp) %*% temp[["A"]]
    }
    
    
    if (need_Sigma) {
      if (sv | tvp_and_covar) {
        temp[["Sigma"]] <- solve(matrix(x[["posterior"]][["u_sigma_inv"]][["coeffs"]][i, (period - 1) * kk + 1:kk], k))
      } else {
        temp[["Sigma"]] <- solve(matrix(x[["posterior"]][["u_sigma_inv"]][["coeffs"]][i, ], k))
      }
    }
    
    # Shock
    if (is.numeric(shock)) {
      temp[["shock"]] <- shock
    } else {
      if (type == "oir") {
        temp[["shock"]] <- diag(chol(temp[["Sigma"]]))[impulse]
      } else {
        temp[["shock"]] <- sqrt(diag(temp[["Sigma"]])[impulse]) 
      }
      
      if (shock == "nsd") {
        temp[["shock"]] <- -temp[["shock"]]
      } 
    }
    
    A[[i]] <- temp
  }
  
  result <- lapply(A, .ir, h = n_ahead, type = type,
                   impulse = impulse, response = response)
  
  result <- t(matrix(unlist(result), n_ahead + 1))
  
  if (cumulative) {
    result <- t(apply(result, 1, cumsum))
  }
  
  if (!keep_draws) {
    ci_low <- (1 - ci) / 2
    ci_high <- 1 - ci_low
    pr <- c(ci_low, .5, ci_high)
    result <- stats::ts(t(apply(result, 2, stats::quantile, probs = pr)), start = 0, frequency = 1) 
  }
  
  class(result) <- append("bvarirf", class(result))
  return(result)
}