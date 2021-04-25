#' Posterior Simulation for BVEC Models
#' 
#' Produces draws from the posterior distributions of Bayesian VEC models.
#' 
#' @param object an object of class \code{"bvarmodel"}, usually, a result of a call to \code{\link{gen_var}}
#' in combination with \code{\link{add_priors}}.
#' 
#' @details The function implements a posterior simulation algorithm, which
#' places identifying restrictions on the cointegration space. Th algorithm is also able to employ
#' stochastic search variable selection (SSVS) as proposed by Geroge et al. (2008) or Bayesian
#' variable selection (BVS) as proposed in Korobilis (2013). Both SSVS and BVS can also be applied to the
#' covariances of the error term. However, the algorithms cannot be applied to cointegration related
#' coefficients, i.e. to the loading matrix \eqn{\alpha} or the cointegration matrix \eqn{beta}.
#' 
#' The implementation primarily follows the description in Koop et al. (2010). However, Chan et al. (2019),
#' George et al. (2008) and Korobilis (2013) were used to implement variable selection algorithms.
#' For all approaches the SUR form of a VEC model is used to obtain posterior draws. The algorithm is implemented
#' in C++ to reduce calculation time.
#' 
#' The function also supports structural BVEC models, where the structural coefficients are estimated from
#' contemporary endogenous variables, which corresponds to the so-called (A-model). Currently, only
#' specifications are supported, where the structural matrix contains ones on its diagonal and all lower
#' triangular elements are freely estimated. Since posterior draws are obtained based on the SUR form of
#' the VAR model, the structural coefficients are drawn jointly with the other coefficients. No identifying
#' restrictionare are made regarding the cointegration matrix.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' @examples
#' 
#' # Get data
#' data("e6")
#' 
#' # Create model
#' model <- gen_vec(e6, p = 4, r = 1,
#'                  const = "unrestricted", seasonal = "unrestricted",
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#'
#' # Obtain posterior draws 
#' object <- bvecpost(model)
#' 
#' @export
bvecpost <- function(object) {
  
  k <- NCOL(object$data[["Y"]])
  p <- object$model$endogen[["lags"]]
  s <- 0
  r <- object$model[["rank"]]
  
  object <- .bvecalg(object) # Use C++ code to draw posteriors

  A0 <- NULL
  if (object$model$structural) {
    pos <- which(lower.tri(diag(1, k)))
    if ("list" %in% class(object$posteriors[["a0"]])) {
      if ("lambda" %in% names(object$posteriors[["a0"]])) {
        A0 <- list("coeffs" = matrix(diag(1, k), k * k, object$model$iterations),
                   "lambda" = matrix(diag(1, k), k * k, object$model$iterations))
        A0[["coeffs"]][pos, ] <- object$posteriors[["a0"]][["coeffs"]]
        A0[["lambda"]][pos, ] <- object$posteriors[["a0"]][["lambda"]]
      }
    } else {
      A0 <- matrix(diag(1, k), k * k, object$model$iterations)
      A0[pos, ] <- object$posteriors[["a0"]]
    }
  }

  tsp_temp <- stats::tsp(object$data$Y)

  w <- stats::ts(as.matrix(object$data$W[, 1:k]), class = c("mts", "ts", "matrix"))
  stats::tsp(w) <- tsp_temp
  dimnames(w)[[2]] <- dimnames(object$data$W)[[2]][1:k]

  m <- 0
  if (!is.null(object$model$exogen)) {
    m <- length(object$model$exogen$variables)
    w_x <- stats::ts(as.matrix(object$data$W[, k + 1:m]), class = c("mts", "ts", "matrix"))
    stats::tsp(w_x) <- tsp_temp
    dimnames(w_x)[[2]] <- dimnames(object$data$W)[[2]][k + 1:m]
  } else {
    w_x <- NULL
  }

  if (!is.null(object$model$deterministic$restricted)) {
    n_r <- length(object$model$deterministic$restricted)
    w_d <- stats::ts(as.matrix(object$data$W[, k + m + 1:n_r]), class = c("mts", "ts", "matrix"))
    stats::tsp(w_d) <- tsp_temp
    dimnames(w_d)[[2]] <- dimnames(object$data$W)[[2]][k + m + n_r]
  } else {
    w_d <- NULL
  }

  x <- NULL
  x_x <- NULL
  x_d <- NULL
  if (!is.null(object$data$X)) {

    if (!is.null(object$model[["endogen"]])) {
      if (object$model$endogen$lags > 1) {
        x <- stats::ts(as.matrix(object$data[["X"]][, 1:(k * (p - 1))]), class = c("mts", "ts", "matrix"))
        stats::tsp(x) <- tsp_temp
        dimnames(x)[[2]] <- dimnames(object$data[["X"]])[[2]][1:(k * (p - 1))]
      }
    }

    if (!is.null(object$model[["exogen"]])) {
      m <- length(object$model$exogen$variables)
      if (object$model$exogen$lags > 1) {
        s <- object$model$exogen$lags
        x_x <- stats::ts(as.matrix(object$data[["X"]][, k * (p - 1) + 1:(m * s)]), class = c("mts", "ts", "matrix"))
        stats::tsp(x_x) <- tsp_temp
        dimnames(x_x)[[2]] <- dimnames(object$data[["X"]])[[2]][(k * (p - 1)) + 1:(m * s)]
      }
    }

    if (!is.null(object$model$deterministic[["unrestricted"]])) {
      n_ur <- length(object$model$deterministic[["unrestricted"]])
      x_d <- stats::ts(as.matrix(object$data[["X"]][, k * (p - 1) + m * s + 1:n_ur]), class = c("mts", "ts", "matrix"))
      stats::tsp(x_d) <- tsp_temp
      dimnames(x_d)[[2]] <- dimnames(object$data[["X"]])[[2]][k * (p - 1) + m * s + 1:n_ur]
    }
  }

  # Create bvar object
  object <- bvec(y = object$data$Y,
                 w = w,
                 w_x = w_x,
                 w_d = w_d,
                 alpha = object$posteriors[["alpha"]],
                 beta = object$posteriors[["beta"]],
                 Pi = object$posteriors[["pi"]],
                 Pi_x = object$posteriors[["pi_exo"]],
                 Pi_d = object$posteriors[["pi_d"]],
                 x = x,
                 x_x = x_x,
                 x_d = x_d,
                 r = r,
                 A0 = A0,
                 Gamma = object$posteriors[["gamma"]],
                 Upsilon = object$posteriors[["upsilon"]],
                 C = object$posteriors[["c"]],
                 Sigma = object$posteriors[["sigma"]],
                 data = NULL, exogen = NULL)

  return(object)
}