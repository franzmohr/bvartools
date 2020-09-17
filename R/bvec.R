#' Bayesian Vector Error Correction Objects
#' 
#' `bvec` is used to create objects of class \code{"bvec"}.
#' 
#' @param data the original time-series object of endogenous variables.
#' @param exogen the original time-series object of unmodelled variables.
#' @param y a time-series object of differenced endogenous variables,
#' usually, a result of a call to \code{\link{gen_vec}}.
#' @param w a time-series object of lagged endogenous variables in levels, which enter the
#' cointegration term, usually, a result of a call to \code{\link{gen_vec}}.
#' @param w_x a time-series object of lagged unmodelled, non-deterministic variables in levels, which enter the
#' cointegration term, usually, a result of a call to \code{\link{gen_vec}}.
#' @param w_d a time-series object of deterministic terms, which enter the
#' cointegration term, usually, a result of a call to \code{\link{gen_vec}}.
#' @param x a time-series object of \eqn{K(p - 1)} differenced endogenous variables.
#' @param x_x a time-series object of \eqn{Ms} differenced unmodelled regressors.
#' @param x_d a time-series object of \eqn{N^{UR}} deterministic terms that do not enter the
#' cointegration term.
#' @param r an integer of the rank of the cointegration matrix.
#' @param A0 either a \eqn{K^2 \times S} matrix of MCMC coefficient draws of structural parameters or
#' a named list, where element \code{coeffs} contains a \eqn{K^2 \times S} matrix of MCMC coefficient
#' draws of structural parameters and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed.
#' @param alpha a \eqn{Kr \times S} matrix of MCMC coefficient draws of the loading matrix \eqn{\alpha}.
#' @param beta a \eqn{((K + M + N^{R})r) \times S} matrix of MCMC coefficient draws of cointegration matrix \eqn{\beta}.
#' @param Pi a \eqn{K^2 \times S} matrix of MCMC coefficient draws of endogenous varaibles in the cointegration matrix.
#' @param Pi_x a \eqn{KM \times S} matrix of MCMC coefficient draws of unmodelled, non-deterministic variables in the cointegration matrix.
#' @param Pi_d a \eqn{KN^{R} \times S} matrix of MCMC coefficient draws of restricted deterministic terms.
#' @param Gamma a \eqn{(p-1)K^2 \times S} matrix of MCMC coefficient draws of differenced lagged endogenous variables or
#' a named list, where element \code{coeffs} contains a \eqn{(p - 1)K^2 \times S} matrix of MCMC coefficient draws
#' of lagged differenced endogenous variables and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed.
#' @param Upsilon an \eqn{sMK \times S} matrix of MCMC coefficient draws of differenced unmodelled, non-deterministic variables or
#' a named list, where element \code{coeffs} contains a \eqn{sMK \times S} matrix of MCMC coefficient draws of
#' unmodelled, non-deterministic variables and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed.
#' @param C an \eqn{KN^{UR} \times S} matrix of MCMC coefficient draws of unrestricted deterministic terms or
#' a named list, where element \code{coeffs} contains a \eqn{KN^{UR} \times S} matrix of MCMC coefficient draws of
#' deterministic terms and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed.
#' @param Sigma a \eqn{K^2 \times S} matrix of MCMC draws for the error variance-covariance matrix or
#' a named list, where element \code{coeffs} contains a \eqn{K^2 \times S} matrix of MCMC draws for the
#' error variance-covariance matrix and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed to the covariances.
#' 
#' @details For the vector error correction model with unmodelled exogenous variables (VECX) 
#' \deqn{A_0 \Delta y_t = \Pi^{+} \begin{pmatrix} y_{t-1} \\ x_{t-1} \\ d^{R}_{t-1} \end{pmatrix} +
#' \sum_{i = 1}^{p-1} \Gamma_i \Delta y_{t-i} +
#' \sum_{i = 0}^{s-1} \Upsilon_i \Delta x_{t-i} +
#' C^{UR} d^{UR}_t + u_t}
#' the function collects the \eqn{S} draws of a Gibbs sampler in a standardised object,
#' where \eqn{\Delta y_t} is a K-dimensional vector of differenced endogenous variables
#' and \eqn{A_0} is a \eqn{K \times K} matrix of structural coefficients.
#' \eqn{\Pi^{+} = \left[ \Pi, \Pi^{x}, \Pi^{d} \right]} is
#' the coefficient matrix of the error correction term, where
#' \eqn{y_{t-1}}, \eqn{x_{t-1}} and \eqn{d^{R}_{t-1}} are the first lags of endogenous,
#' exogenous variables in levels and restricted deterministic terms, respectively.
#' \eqn{\Pi}, \eqn{\Pi^{x}}, and \eqn{\Pi^{d}} are the corresponding coefficient matrices, respectively.
#' \eqn{\Gamma_i} is a coefficient matrix of lagged differenced endogenous variabels.
#' \eqn{\Delta x_t} is an M-dimensional vector of unmodelled, non-deterministic variables
#' and \eqn{\Upsilon_i} its corresponding coefficient matrix. \eqn{d_t} is an
#' \eqn{N^{UR}}-dimensional vector of unrestricted deterministics and \eqn{C^{UR}}
#' the corresponding coefficient matrix.
#' \eqn{u_t} is an error term with \eqn{u_t \sim N(0, \Sigma_u)}.
#' 
#' The draws of the different coefficient matrices provided in \code{alpha}, \code{beta},
#' \code{Pi}, \code{Pi_x}, \code{Pi_d}, \code{A0}, \code{Gamma}, \code{Ypsilon},
#' \code{C} and \code{Sigma} have to correspond to the same MCMC iteration.
#' 
#' @return An object of class \code{"gvec"} containing the following components, if specified:
#' \item{data}{the original time-series object of endogenous variables.}
#' \item{exogen}{the original time-series object of unmodelled variables.}
#' \item{y}{a time-series object of differenced endogenous variables.}
#' \item{w}{a time-series object of lagged endogenous variables in levels, which enter the
#' cointegration term.}
#' \item{w_x}{a time-series object of lagged unmodelled, non-deterministic variables in levels, which enter the
#' cointegration term.}
#' \item{w_d}{a time-series object of deterministic terms, which enter the
#' cointegration term.}
#' \item{x}{a time-series object of \eqn{K(p - 1)} differenced endogenous variables}
#' \item{x_x}{a time-series object of \eqn{Ms} differenced unmodelled regressors.}
#' \item{x_d}{a time-series object of \eqn{N^{UR}} deterministic terms that do not enter the
#' cointegration term.}
#' \item{A0}{an \eqn{S \times K^2} "mcmc" object of coefficient draws of structural parameters.}
#' \item{A0_lambda}{an \eqn{S \times K^2} "mcmc" object of inclusion parameters for coefficients
#' corresponding to structural parameters.}
#' \item{alpha}{an \eqn{S \times Kr} "mcmc" object of coefficient draws of loading parameters.}
#' \item{beta}{an \eqn{S \times ((K + M + N^{R})r)} "mcmc" object of coefficient draws of cointegration parameters.}
#' \item{Pi}{an \eqn{S \times K^2} "mcmc" object of coefficient draws of endogenous variables in the cointegration matrix.}
#' \item{Pi_x}{an \eqn{S \times KM} "mcmc" object of coefficient draws of unmodelled, non-deterministic variables in the cointegration matrix.}
#' \item{Pi_d}{an \eqn{S \times KN^{R}} "mcmc" object of coefficient draws of unrestricted deterministic variables in the cointegration matrix.}
#' \item{Gamma}{an \eqn{S \times (p-1)K^2} "mcmc" object of coefficient draws of differenced lagged endogenous variables.}
#' \item{Gamma_lamba}{an \eqn{S \times (p-1)K^2} "mcmc" object of inclusion parameters for coefficients
#' corresponding to differenced lagged endogenous variables.}
#' \item{Upsilon}{an \eqn{S \times sMK} "mcmc" object of coefficient draws of differenced unmodelled variables.}
#' \item{Upsilon_lambda}{an \eqn{S \times sMK} "mcmc" object of inclusion parameters for coefficients
#' corresponding to differenced unmodelled, non-deterministic variables.}
#' \item{C}{an \eqn{S \times KN^{UR}} "mcmc" object of coefficient draws of deterministic terms that
#' do not enter the cointegration term.}
#' \item{C_lambda}{an \eqn{S \times KN^{UR}} "mcmc" object of inclusion parameters for coefficients
#' corresponding to deterministic terms, that do not enter the conitegration term.}
#' \item{Sigma}{an \eqn{S \times K^2} "mcmc" object of variance-covariance draws.}
#' \item{Sigma_lambda}{an \eqn{S \times K^2} "mcmc" object inclusion parameters for the variance-covariance
#' matrix.}
#' \item{specifications}{a list containing information on the model specification.}
#' 
#' @examples 
#' 
#' # Load data
#' data("e6")
#' # Generate model
#' data <- gen_vec(e6, p = 4, r = 1, const = "unrestricted", season = "unrestricted")
#' # Obtain data matrices
#' y <- t(data$data$Y)
#' w <- t(data$data$W)
#' x <- t(data$data$X)
#' 
#' # Reset random number generator for reproducibility
#' set.seed(1234567)
#' 
#' iterations <- 400 # Number of iterations of the Gibbs sampler
#' # Chosen number of iterations should be much higher, e.g. 30000.
#' 
#' burnin <- 100 # Number of burn-in draws
#' draws <- iterations + burnin
#' 
#' r <- 1 # Set rank
#' 
#' tt <- ncol(y) # Number of observations
#' k <- nrow(y) # Number of endogenous variables
#' k_w <- nrow(w) # Number of regressors in error correction term
#' k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms
#' 
#' k_alpha <- k * r # Number of elements in alpha
#' k_beta <- k_w * r # Number of elements in beta
#' k_gamma <- k * k_x
#' 
#' # Set uninformative priors
#' a_mu_prior <- matrix(0, k_x * k) # Vector of prior parameter means
#' a_v_i_prior <- diag(0, k_x * k) # Inverse of the prior covariance matrix
#' 
#' v_i <- 0
#' p_tau_i <- diag(1, k_w)
#' 
#' u_sigma_df_prior <- r # Prior degrees of freedom
#' u_sigma_scale_prior <- diag(0, k) # Prior covariance matrix
#' u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom
#' 
#' # Initial values
#' beta <- matrix(c(1, -4), k_w, r)
#' u_sigma_i <- diag(1 / .0001, k)
#' g_i <- u_sigma_i
#' 
#' # Data containers
#' draws_alpha <- matrix(NA, k_alpha, iterations)
#' draws_beta <- matrix(NA, k_beta, iterations)
#' draws_pi <- matrix(NA, k * k_w, iterations)
#' draws_gamma <- matrix(NA, k_gamma, iterations)
#' draws_sigma <- matrix(NA, k^2, iterations)
#' 
#' # Start Gibbs sampler
#' for (draw in 1:draws) {
#'   # Draw conditional mean parameters
#'   temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = u_sigma_i,
#'                          v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
#'                          gamma_mu_prior = a_mu_prior,
#'                          gamma_v_i_prior = a_v_i_prior)
#'   alpha <- temp$alpha
#'   beta <- temp$beta
#'   Pi <- temp$Pi
#'   gamma <- temp$Gamma
#'   
#'   # Draw variance-covariance matrix
#'   u <- y - Pi %*% w - matrix(gamma, k) %*% x
#'   u_sigma_scale_post <- solve(tcrossprod(u) +
#'      v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
#'   u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
#'   u_sigma <- solve(u_sigma_i)
#'   
#'   # Update g_i
#'   g_i <- u_sigma_i
#'   
#'   # Store draws
#'   if (draw > burnin) {
#'     draws_alpha[, draw - burnin] <- alpha
#'     draws_beta[, draw - burnin] <- beta
#'     draws_pi[, draw - burnin] <- Pi
#'     draws_gamma[, draw - burnin] <- gamma
#'     draws_sigma[, draw - burnin] <- u_sigma
#'   }
#' }
#' 
#' # Number of non-deterministic coefficients
#' k_nondet <- (k_x - 4) * k
#' 
#' # Generate bvec object
#' bvec_est <- bvec(y = data$data$Y, w = data$data$W,
#'                  x = data$data$X[, 1:6],
#'                  x_d = data$data$X[, 7:10],
#'                  Pi = draws_pi,
#'                  Gamma = draws_gamma[1:k_nondet,],
#'                  C = draws_gamma[(k_nondet + 1):nrow(draws_gamma),],
#'                  Sigma = draws_sigma)
#' 
#' @export
bvec <- function(y, alpha = NULL, beta = NULL, r = NULL,
                 Pi = NULL,
                 Pi_x = NULL, Pi_d = NULL,
                 w = NULL, w_x = NULL, w_d = NULL,
                 Gamma = NULL, Upsilon = NULL, C = NULL,
                 x = NULL, x_x = NULL, x_d = NULL, 
                 A0 = NULL, Sigma = NULL,
                 data = NULL, exogen = NULL) {

  result <- NULL
  if (is.null(Pi) & is.null(Pi_x) & is.null(Pi_d) & is.null(r)) {
    stop("If arguments 'Pi', 'Pi_x' and 'Pi_d' are not specified, at least argument 'r' must be set to zero.")
  }
  if (is.null(Pi) & is.null(Pi_x) & is.null(Pi_d) & !is.null(r)) {
    if (r > 0) {
      stop("Non of the arguments 'Pi', 'Pi_x' and 'Pi_d' is specified although 'r' is larger than zero.") 
    }
  }
  if (is.null(r)) {
    r <- NA
  }
  
  if (!"ts" %in% class(y)) {
    stop("'Argument 'y' must be of class 'ts'.")
  }
  result$y <- y
  k <- NCOL(y)
  
  if(!is.null(w)) {
    if (!"ts" %in% class(w)) {
      stop("'Argument 'w' must be of class 'ts'.")
    }
    result$w <- w
  }
  
  m <- NULL
  s <- NULL
  if(!is.null(w_x)) {
    if (!"ts" %in% class(w_x)) {
      stop("'Argument 'w_x' must be of class 'ts'.")
    }
    result$w_x <- w_x
    m <- ncol(w_x)
  }
  
  if(!is.null(w_d)) {
    if (!"ts" %in% class(w_d)) {
      stop("'Argument 'w_d' must be of class 'ts'.")
    }
    result$w_d <- w_d
  }
  
  if(!is.null(x)) {
    if (!"ts" %in% class(x)) {
      stop("'Argument 'x' must be of class 'ts'.")
    }
    result$x <- x
  }
  
  if(!is.null(x_x)) {
    if (!"ts" %in% class(x_x)) {
      stop("'Argument 'x_x' must be of class 'ts'.")
    }
    result$x_x <- x_x
  }
  
  if(!is.null(x_d)) {
    if (!"ts" %in% class(x_d)) {
      stop("'Argument 'x_d' must be of class 'ts'.")
    }
    result$x_d <- x_d
  }
  
  if(!is.null(A0)) {
    if ("list" %in% class(A0)) {
      result[["A0"]] <- coda::mcmc(t(A0$coeffs))
      result[["A0_lambda"]] <- coda::mcmc(t(A0$lambda))
    } else {
      result[["A0"]] <- coda::mcmc(t(A0)) 
    }
    if (NCOL(result[["A0"]]) != k^2) {
      stop("Row number of argument 'A0' is not the square of the number of endogenous variables.")
    }
  }
  
  if(!is.null(alpha)) {
    if ("list" %in% class(alpha)) {
      result[["alpha"]] <- coda::mcmc(t(alpha$coeffs))
      result[["alpha_lambda"]] <- coda::mcmc(t(alpha$lambda))
    } else {
      result[["alpha"]] <- coda::mcmc(t(alpha)) 
    }
  }
  
  if(!is.null(beta)) {
    if ("list" %in% class(beta)) {
      result[["beta"]] <- coda::mcmc(t(beta$coeffs))
      result[["beta_lambda"]] <- coda::mcmc(t(beta$lambda))
    } else {
      result[["beta"]] <- coda::mcmc(t(beta)) 
    }
  }
  
  if(!is.null(Pi)) {
    if ("list" %in% class(Pi)) {
      result[["Pi"]] <- coda::mcmc(t(Pi$coeffs))
      result[["Pi_lambda"]] <- coda::mcmc(t(Pi$lambda))
    } else {
      result[["Pi"]] <- coda::mcmc(t(Pi)) 
    }
    if (NCOL(result[["Pi"]]) != k^2) {
      stop("Row number of argument 'Pi' is not the square of the number of endogenous variables.")
    }
  }
  
  if(!is.null(Pi_x)) {
    if ("list" %in% class(Pi_x)) {
      result[["Pi_x"]] <- coda::mcmc(t(Pi_x$coeffs))
      result[["Pi_x_lambda"]] <- coda::mcmc(t(Pi_x$lambda))
    } else {
      result[["Pi_x"]] <- coda::mcmc(t(Pi_x)) 
    }
    n_x <- NCOL(result[["Pi_x"]])
    if (n_x %% k == 0) {
      if (is.null(m)) {
        m <- n_x / k 
      }
    } else {
      stop("Row number of argument 'Pi_x' is not a multiple of the number of endogenous variables.")
    }
    s <- 0
  }
  
  if(!is.null(Pi_d)) {
    if ("list" %in% class(Pi_d)) {
      result[["Pi_d"]] <- coda::mcmc(t(Pi_d$coeffs))
      result[["Pi_d_lambda"]] <- coda::mcmc(t(Pi_d$lambda))
    } else {
      result[["Pi_d"]] <- coda::mcmc(t(Pi_d)) 
    }
    n_c_r <- NCOL(result[["Pi_d"]])
    if (n_c_r %% k == 0) {
      n_r <- n_c_r / k
    } else {
      stop("Row number of argument 'Pi_d' is not a multiple of the number of endogenous variables.")
    }
  }
  
  if(!is.null(Gamma)) {
    if ("list" %in% class(Gamma)) {
      result[["Gamma"]] <- coda::mcmc(t(Gamma$coeffs))
      result[["Gamma_lambda"]] <- coda::mcmc(t(Gamma$lambda))
    } else {
      result[["Gamma"]] <- coda::mcmc(t(Gamma)) 
    }
    n_gamma <- NCOL(result[["Gamma"]])
    if (n_gamma %% k == 0) {
      p <- n_gamma / k^2
    } else {
      stop("Row number of argument 'Gamma' is not a multiple of the number of endogenous variables.")
    }
  } else {
    p <- 0
  }
  
  if(!is.null(Upsilon)) {
    if ("list" %in% class(Upsilon)) {
      result[["Upsilon"]] <- coda::mcmc(t(Upsilon$coeffs))
      result[["Upsilon_lambda"]] <- coda::mcmc(t(Upsilon$lambda))
    } else {
      result[["Upsilon"]] <- coda::mcmc(t(Upsilon))
    }
    n_x <- NCOL(result[["Upsilon"]])
    if (n_x %% k == 0) {
      if (!is.null(m)) {
        s <- n_x / k / m - 1  
      }
    } else {
      stop("Row number of argument 'Upsilon' is not a multiple of the number of endogenous variables.")
    }
  }
  if (!is.null(s)) {
    s <- s + 1
  }
  
  if(!is.null(C)) {
    if ("list" %in% class(C)) {
      result[["C"]] <- coda::mcmc(t(C$coeffs))
      result[["C_lambda"]] <- coda::mcmc(t(C$lambda))
    } else {
      result[["C"]] <- coda::mcmc(t(C))
    }
    n_x <- NCOL(result[["C"]])
    if (n_x %% k == 0) {
      n_ur <- n_x / k
    } else {
      stop("Row number of argument 'C' is not a multiple of the number of endogenous variables.")
    }
  }
  
  if(!is.null(Sigma)) {
    if ("list" %in% class(Sigma)) {
      result[["Sigma"]] <- coda::mcmc(t(Sigma$coeffs))
      result[["Sigma_lambda"]] <- coda::mcmc(t(Sigma$lambda))
    } else {
      result[["Sigma"]] <- coda::mcmc(t(Sigma))
    }
    n_x <- NCOL(result[["Sigma"]])
    if (n_x != k^2) {
      stop("Row number of argument 'Sigma' is not a multiple of the number of endogenous variables.")
    }
  }
  
  result$specifications <- list("dims" = c("K" = k, "M" = m),
                                "lags" = c("p" = p + 1, "s" = s),
                                "rank" = r)
  
  class(result) <- "bvec"
  return(result)
}