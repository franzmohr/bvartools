#' Impulse Response Function
#' 
#' Computes the impulse response coefficients of an object of class \code{"bvar"} for
#' \code{n.ahead} steps.
#' 
#' @param x an object of class \code{"bvar"}, usually, a result of a call to
#' \code{\link{bvar}} or \code{\link{bvec_to_bvar}}.
#' @param impulse name of the impulse variable.
#' @param response name of the response variable.
#' @param n.ahead number of steps ahead.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param type type of the impulse resoponse. Possible choices are forecast error \code{"feir"}
#' (default), orthogonalised \code{"oir"}, structural \code{"sir"}, generalised \code{"gir"},
#' and structural generalised \code{"sgir"} impulse responses.
#' @param cumulative logical specifying whether a cumulative IRF should be calculated.
#' @param keep_draws logical specifying whether the function should return all draws of
#' the posterior impulse response function. Defaults to \code{FALSE} so that
#' the median and the credible intervals of the posterior draws are returned.
#' @param period integer. Index of the period of a TVP VAR, for which IRFs should be generated.
#' Only used for TVP models. Default is \code{NULL} so that only the most recent time period
#' is used.
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
#' one in its \eqn{j^{th}} element and zero otherwise. If the \code{"bvar"} object does not contain draws
#' of \eqn{A_0}, it is assumed to be an identity matrix.
#' 
#' @return A time-series object of class \code{"bvarirf"} and if \code{keep_draws = TRUE} a simple matrix.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model data
#' model <- gen_var(e1, p = 2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Obtain IR
#' ir <- irf(object, impulse = "invest", response = "cons")
#' 
#' # Plot IR
#' plot(ir)
#'
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Pesaran, H. H., Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' @export
irf <- function(x, impulse = NULL, response = NULL, n.ahead = 5, ci = .95, shock = 1,
                type = "feir", cumulative = FALSE, keep_draws = FALSE, period = NULL) {
  
  # Dev specs
  # rm(list = ls()[-which(ls() == "x")]); impulse = "r"; response = "Dp"; n.ahead = 20; ci = .95; type = "oir"; cumulative = FALSE; keep_draws = FALSE; period <- NULL
  
  if (!type %in% c("feir", "oir", "gir", "sir", "sgir")) {
    stop("Argument 'type' not known.")
  }
  
  if (!"bvar" %in% class(x)) {
    stop("Argument 'x' must be of class 'bvar'.")
  }
  
  if (is.null(x$y) | is.null(dimnames(x$y)[[2]])) {
    stop("Argument 'x' must include a named matrix of endogenous variables.")
  }
  
  if (is.null(x[["A"]]) & !type %in% c("sir", "sgir")) {
    stop("Impulse responses only supported for models with p > 0, i.e. argument 'x' must contain element 'A', or structural models.")
  }
  
  need_A0 <- FALSE
  if (type %in% c("sgir", "sir")) {
    if (is.null(x[["A0"]])) {
      stop("Structural IR requires that draws of 'A0' are contained in the 'bvar' x.")
    }
    need_A0 <- TRUE
  }
  
  if (!(is.numeric(shock) | shock %in% c("sd", "nsd"))) {
    stop("Invalid specification of argument 'shock'.")
  }
  
  if (type %in% c("oir", "gir", "sgir") | shock %in% c("sd", "nsd")) {
    if (is.null(x[["Sigma"]])) {
      stop("OIR, GIR, SGIR require that the 'bvar' x contains draws of 'Sigma'.")
    }
    need_Sigma <- TRUE
  } else {
    need_Sigma <- FALSE
  }
  
  impulse <- which(dimnames(x$y)[[2]] == impulse)
  if (length(impulse) == 0){stop("Impulse variable not available.")}
  response <- which(dimnames(x$y)[[2]] == response)
  if (length(response) == 0){stop("Response variable not available.")}
  
  k <- NCOL(x$y)
  tt <- NROW(x$y)
  tvp <- x[["specifications"]][["tvp"]]
  if (any(unlist(tvp))) {
    if (is.null(period)) {
      period <- tt
    } else {
      if (period > tt | period < 1) {
        stop("Implausible specification of argument 'period'.")
      }
    }
  }
  
  store <- NA
  vars <- c("A0", "A", "B", "C", "Sigma")
  for (i in vars) {
    if (is.na(store)) {
      if (!is.null(x[[i]])) {
        if (x[["specifications"]][["tvp"]][[i]]) {
          store <- nrow(x[[i]][[1]])
        } else {
          store <- nrow(x[[i]]) 
        }
      }   
    }
  }
  
  A <- NULL
  for (i in 1:store) {
    temp <- NULL
    if (!is.null(x[["A"]])) {
      if (x[["specifications"]][["tvp"]][["A"]]) {
        temp[["A"]] <- matrix(x[["A"]][[period]][i, ], k)
      } else {
        temp[["A"]] <- matrix(x[["A"]][i, ], k) 
      }
    } else {
      temp[["A"]] <- matrix(0, k, k)
    }
    if (need_Sigma) {
      if (x[["specifications"]][["tvp"]][["Sigma"]]) {
        temp[["Sigma"]] <- matrix(x[["Sigma"]][[period]][i, ], k)
      } else {
        temp[["Sigma"]] <- matrix(x[["Sigma"]][i, ], k) 
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
      
    if (need_A0) {
      if (x[["specifications"]][["tvp"]][["A0"]]) {
        temp[["A0"]] <- matrix(x[["A0"]][[period]][i, ], k)
      } else {
        temp[["A0"]] <- matrix(x[["A0"]][i, ], k) 
      }
    }
    
    A[[i]] <- temp
  }
  
  result <- lapply(A, .ir, h = n.ahead, type = type,
                   impulse = impulse, response = response)
  
  result <- t(matrix(unlist(result), n.ahead + 1))
  
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