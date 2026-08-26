#' Bayesian Vector Error Correction Objects
#'
#' \code{bvec} collects the posterior draws of a vector error correction model
#' in an object of class 'bvecmodel'.
#'
#' @param data the original time-series object of endogenous variables in levels.
#' If \code{NULL} (default), it is reconstructed from arguments \code{y} and \code{w}.
#' @param exogen the original time-series object of unmodelled variables in levels.
#' If \code{NULL} (default), it is reconstructed from arguments \code{w_x} and \code{x_x}.
#' @param z a \eqn{TK \times K (r + (p-1)K + sM + N)} data matrix, usually, a result of a
#' call to \code{\link{create_bvecmodel}}. If \code{NULL} (default), it is generated
#' from the provided data.
#' @param y a time-series object of differenced endogenous variables,
#' usually, a result of a call to \code{\link{create_bvecmodel}}.
#' @param w a time-series object of lagged endogenous variables in levels, which enter the
#' cointegration term, usually, a result of a call to \code{\link{create_bvecmodel}}.
#' @param w_x a time-series object of lagged unmodelled, non-deterministic variables in levels, which enter the
#' cointegration term, usually, a result of a call to \code{\link{create_bvecmodel}}.
#' @param w_d a time-series object of deterministic terms, which enter the
#' cointegration term, usually, a result of a call to \code{\link{create_bvecmodel}}.
#' @param x a time-series object of \eqn{K(p - 1)} differenced endogenous variables.
#' @param x_x a time-series object of \eqn{Ms} differenced unmodelled regressors.
#' @param x_d a time-series object of \eqn{N^{UR}} deterministic terms that do not enter the
#' cointegration term.
#' @param r an integer of the rank of the cointegration matrix. If \code{NULL}
#' (default), it is obtained from argument \code{alpha}.
#' @param A0 either a \eqn{K^2 \times S} matrix of MCMC coefficient draws of structural parameters or
#' a named list, where element \code{coeffs} contains a \eqn{K^2 \times S} matrix of MCMC coefficient
#' draws of structural parameters and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed.
#' @param alpha a \eqn{Kr \times S} matrix of MCMC coefficient draws of the loading matrix \eqn{\alpha}.
#' @param beta a \eqn{Kr \times S} matrix of MCMC coefficient draws of cointegration matrix \eqn{\beta}
#' corresponding to the endogenous variables of the model.
#' @param beta_x a \eqn{Mr \times S} matrix of MCMC coefficient draws of cointegration matrix \eqn{\beta}
#' corresponding to unmodelled, non-deterministic variables.
#' @param beta_d a \eqn{N^{R}r \times S} matrix of MCMC coefficient draws of cointegration matrix \eqn{\beta}
#' corresponding to restricted deterministic terms.
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
#' parameters in case variable selection algorithms were employed to the covariances. For models with
#' stochastic volatility the matrix must be \eqn{TK^2 \times S}.
#' @param error a character specifying the model that was used for the estimation
#' of the covariance matrix of the error term. If \code{NULL} (default), it is
#' inferred from the draws in argument \code{Sigma}, which can only distinguish
#' a constant covariance matrix from a time varying one. See
#' \code{\link{create_bvecmodel}}.
#' @param varsel a character specifying the variable selection algorithm that was
#' employed. If \code{NULL} (default), it is inferred from the presence of draws
#' of inclusion parameters. See \code{\link{create_bvecmodel}}.
#' @param iterations an integer of the number of MCMC draws in the provided
#' coefficient matrices. If \code{NULL} (default), it is obtained from those
#' matrices.
#' @param burnin an integer of the number of MCMC draws that were used to
#' initialise the sampler. Defaults to zero, since the draws provided to the
#' function are assumed to exclude them already.
#'
#' @details For the vector error correction model with unmodelled exogenous variables (VECX)
#' \deqn{A_0 \Delta y_t = \alpha \beta^\prime \begin{pmatrix} y_{t-1} \\ x_{t-1} \\ d^{R}_{t-1} \end{pmatrix} +
#' \sum_{i = 1}^{p-1} \Gamma_i \Delta y_{t-i} +
#' \sum_{i = 0}^{s-1} \Upsilon_i \Delta x_{t-i} +
#' C^{UR} d^{UR}_t + u_t}
#' the function collects the \eqn{S} draws of a Gibbs sampler in a standardised object,
#' where \eqn{\Delta y_t} is a K-dimensional vector of differenced endogenous variables
#' and \eqn{A_0} is a \eqn{K \times K} matrix of structural coefficients.
#' \eqn{\alpha} is the \eqn{K \times r} loading matrix and \eqn{\beta} the
#' \eqn{(K + M + N^{R}) \times r} cointegration matrix of the error correction term,
#' where \eqn{y_{t-1}}, \eqn{x_{t-1}} and \eqn{d^{R}_{t-1}} are the first lags of endogenous,
#' exogenous variables in levels and restricted deterministic terms, respectively.
#' \eqn{\Gamma_i} is a coefficient matrix of lagged differenced endogenous variabels.
#' \eqn{\Delta x_t} is an M-dimensional vector of unmodelled, non-deterministic variables
#' and \eqn{\Upsilon_i} its corresponding coefficient matrix. \eqn{d_t} is an
#' \eqn{N^{UR}}-dimensional vector of unrestricted deterministics and \eqn{C^{UR}}
#' the corresponding coefficient matrix.
#' \eqn{u_t} is an error term with \eqn{u_t \sim N(0, \Sigma_u)}.
#'
#' For time varying parameter and stochastic volatility models the respective coefficients and
#' error covariance matrix of the above model are assumed to be time varying, respectively.
#'
#' The draws of the different coefficient matrices provided in \code{alpha}, \code{beta},
#' \code{beta_x}, \code{beta_d}, \code{A0}, \code{Gamma}, \code{Upsilon},
#' \code{C} and \code{Sigma} have to correspond to the same MCMC iteration.
#'
#' The result is the same kind of object as the output of
#' \code{\link{create_bvecmodel}} in combination with
#' \code{\link{draw_posterior}}, so it can be used with the methods of class
#' 'bvecmodel' such as \code{\link{vec_to_var.bvecmodel}} or
#' \code{\link{summary.bvecmodel}}. Accordingly, the draws are stored in the
#' parameterisation those methods expect. The error correction term is described
#' by \eqn{\alpha} and \eqn{\beta} and not by \eqn{\Pi = \alpha \beta^\prime},
#' since the latter cannot be decomposed into the former. The draws of
#' \eqn{\alpha} are collected in the same matrix as the coefficients of the
#' remaining regressors and the draws of \eqn{\beta} of all cointegration
#' relations are collected in a matrix of their own. Draws of the error term are
#' stored as draws of the inverse of \eqn{\Sigma_u} and draws of \eqn{A_0} are
#' reduced to the free elements below its diagonal, which are the elements that
#' are estimated.
#'
#' Since a model that is constructed from posterior draws does not contain the
#' specification of its priors, the elements \code{priors} and \code{initial} of
#' the resulting object are empty. Use \code{\link{add_priors}} and
#' \code{\link{add_initial_values}} if the model should be re-estimated.
#'
#' @return An object of class \code{"bvecmodel"} containing the following components:
#' \item{model}{a list containing information on the model specification.}
#' \item{data}{a list of data objects. Element \code{original} contains the
#' original time-series objects of endogenous and unmodelled variables in levels.
#' Element \code{train} contains the time-series objects \code{y} of differenced
#' endogenous variables, \code{w} of the regressors of the error correction term
#' and \code{x} of the remaining regressors as well as the matrix \code{z} of
#' regressors in SUR form.}
#' \item{posterior}{a list of posterior draws. Element \code{a} contains an
#' \eqn{S \times (Kr + K((p - 1)K + sM + N^{UR}) + K(K - 1) / 2)} "mcmc" object of
#' the draws of the loading matrix and the coefficients of the remaining
#' regressors and, if provided, the corresponding draws of inclusion parameters
#' in element \code{lambda} and of the error covariance matrix of the state
#' equation in element \code{sigma}. Element \code{beta} contains an
#' \eqn{S \times ((K + M + N^{R})r)} "mcmc" object of the draws of the
#' cointegration matrix. Element \code{u_sigma_inv} contains an
#' \eqn{S \times K^2} "mcmc" object of the draws of the inverse of the error
#' variance-covariance matrix.}
#'
#' @examples
#'
#' # Load data
#' data("e6")
#' # Generate model
#' model <- create_bvecmodel(e6, p = 4, r = 1, const = "unrestricted", seasonal = "unrestricted")
#' # Obtain data matrices
#' y <- t(model$data$train$y)
#' w <- t(model$data$train$w)
#' x <- t(model$data$train$x)
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
#' a_mu_prior <- matrix(0, k_alpha + k_gamma) # Vector of prior parameter means
#' a_v_i_prior <- diag(0, k_alpha + k_gamma) # Inverse of the prior covariance matrix
#'
#' v_i <- 0
#' p_tau_i <- diag(1, k_w)
#'
#' u_sigma_df_prior <- r # Prior degrees of freedom
#' u_sigma_scale_prior <- diag(0.01, k) # Prior covariance matrix
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
#' draws_gamma <- matrix(NA, k_gamma, iterations)
#' draws_sigma <- matrix(NA, k^2, iterations)
#'
#' # Start Gibbs sampler
#' for (draw in 1:draws) {
#'   # Draw conditional mean parameters
#'   temp <- post_coint_kls(y = y, beta = beta, w = w, sigma_i = u_sigma_i,
#'                          v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
#'                          x = x,
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
#'     draws_gamma[, draw - burnin] <- gamma
#'     draws_sigma[, draw - burnin] <- u_sigma
#'   }
#' }
#'
#' # Number of non-deterministic coefficients
#' k_nondet <- (k_x - 4) * k
#'
#' # Generate bvecmodel object
#' bvec_est <- bvec(y = model$data$train$y, w = model$data$train$w,
#'                  x = model$data$train$x[, 1:6],
#'                  x_d = model$data$train$x[, 7:10],
#'                  r = r,
#'                  alpha = draws_alpha,
#'                  beta = draws_beta,
#'                  Gamma = draws_gamma[1:k_nondet,],
#'                  C = draws_gamma[(k_nondet + 1):nrow(draws_gamma),],
#'                  Sigma = draws_sigma)
#'
#' @export
bvec <- function(y, alpha = NULL, beta = NULL, beta_x = NULL, beta_d = NULL, r = NULL,
                 w = NULL, w_x = NULL, w_d = NULL,
                 Gamma = NULL, Upsilon = NULL, C = NULL,
                 x = NULL, x_x = NULL, x_d = NULL,
                 A0 = NULL, Sigma = NULL,
                 data = NULL, exogen = NULL, z = NULL,
                 error = NULL, varsel = NULL, iterations = NULL, burnin = 0) {

  for (i in c("y", "w", "w_x", "w_d", "x", "x_x", "x_d")) {
    temp <- get(i)
    if (!is.null(temp)) {
      if (!"ts" %in% class(temp)) {
        stop("Argument '", i, "' must be of class 'ts'.")
      }
    }
  }

  y <- .name_ts(y, "d.y")
  k <- NCOL(y)
  tt <- NROW(y)

  # Regressors of the error correction term ----

  w <- .name_ts(w, "l.y")
  w_x <- .name_ts(w_x, "l.x")
  w_d <- .name_ts(w_d, "l.d")

  m <- ifelse(is.null(w_x), 0L, as.integer(NCOL(w_x)))
  n_restricted <- ifelse(is.null(w_d), 0L, as.integer(NCOL(w_d)))

  if (!is.null(w)) {
    if (NCOL(w) != k) {
      stop("Argument 'w' must contain one column per endogenous variable.")
    }
  }
  k_beta <- as.integer(k + m + n_restricted)

  # Cointegration rank ----

  block_alpha <- .split_draws(alpha, "alpha")
  if (is.null(block_alpha)) {
    if (is.null(r)) {
      stop("Either argument 'alpha' or argument 'r' must be specified.")
    }
    if (r > 0) {
      stop("Argument 'alpha' must be specified if argument 'r' is larger than zero.")
    }
    rank <- 0L
  } else {
    block_alpha <- .block_dims(block_alpha, tt, k, "alpha")
    rank <- as.integer(block_alpha[["n"]] / k)
    if (!is.null(r)) {
      if (r != rank) {
        stop("The number of rows of argument 'alpha' does not correspond to the rank in argument 'r'.")
      }
    }
    if (is.null(beta)) {
      stop("Argument 'beta' must be specified if argument 'alpha' is specified.")
    }
    if (is.null(w)) {
      stop("Argument 'w' must be specified for models with a cointegration rank larger than zero.")
    }
  }

  # Remaining coefficient draws ----

  block_gamma <- .block_dims(.split_draws(Gamma, "Gamma"), tt, k * k, "Gamma")
  block_upsilon <- .block_dims(.split_draws(Upsilon, "Upsilon"), tt, k, "Upsilon")
  block_c <- .block_dims(.split_draws(C, "C"), tt, k, "C")
  block_sigma <- .block_dims(.split_draws(Sigma, "Sigma"), tt, k * k, "Sigma")

  block_a0 <- .split_draws(A0, "A0")
  structural <- !is.null(block_a0)
  if (structural) {
    if (k == 1) {
      stop("Argument 'A0' cannot be used for models with a single endogenous variable.")
    }
    block_a0 <- .block_dims(block_a0, tt, k * k, "A0")
    block_a0 <- .a0_free_elements(block_a0, k)
    block_a0[["n"]] <- k * (k - 1) / 2
  }

  draws <- .count_draws(list(block_alpha, block_gamma, block_upsilon,
                             block_c, block_a0, block_sigma))

  # Dimensions of the model ----

  # The lag orders of a VEC model are those of the corresponding VAR model, so
  # they exceed the number of blocks of differenced regressors by one
  p <- if (is.null(block_gamma)) 1L else as.integer(block_gamma[["n"]] / (k * k) + 1)
  n <- if (is.null(block_c)) 0L else as.integer(block_c[["n"]] / k)

  s <- 0L
  if (!is.null(block_upsilon)) {
    if (m == 0) {
      stop("Argument 'w_x' must be specified if argument 'Upsilon' is specified.")
    }
    s <- as.integer(block_upsilon[["n"]] / (k * m))
  }

  # Data ----

  ect <- .cbind_ts(list(w, w_x, w_d))
  regressors <- .cbind_ts(list(x, x_x, x_d))

  if (!is.null(regressors)) {
    if (NCOL(regressors) != k * (p - 1) + m * s + n) {
      stop("The number of regressors does not correspond to the number of provided coefficient draws.")
    }
  }

  # The levels of the endogenous variables are the sum of their differences and
  # their first lag in the error correction term
  if (is.null(data) & !is.null(w)) {
    data <- .subset_ts(y + w, 1:k, sub("^l\\.", "", dimnames(w)[[2]]))
  }
  data <- .name_ts(data, "y")

  # The same holds for the unmodelled variables, whose current differences are
  # the first block of the differenced unmodelled regressors
  if (is.null(exogen) & m > 0 & s > 0 & !is.null(x_x)) {
    exogen <- .subset_ts(w_x + x_x[, 1:m, drop = FALSE], 1:m, sub("^l\\.", "", dimnames(w_x)[[2]]))
  }
  exogen <- .name_ts(exogen, "x")

  if (is.null(z)) {
    z <- NULL
    if (!is.null(regressors)) {
      z <- kronecker(regressors, diag(1, k))
    }
    # The columns of the loading matrix are estimated conditional on the
    # cointegration matrix, so their regressors are not known in advance
    if (rank > 0) {
      z <- cbind(matrix(NA_real_, tt * k, rank * k), z)
    }
    if (structural) {
      y_A0 <- kronecker(-y, diag(1, k))
      pos <- NULL
      for (i in 1:k) {
        pos <- c(pos, (i - 1) * k + 1:i)
      }
      z <- cbind(z, y_A0[, -pos])
    }
    if (!is.null(z)) {
      dimnames(z) <- NULL
    }
  }

  # Model specification ----

  sv <- isTRUE(block_sigma[["tvp"]])
  tvp <- any(vapply(list(block_alpha, block_gamma, block_upsilon, block_c, block_a0),
                    function(x) {isTRUE(x[["tvp"]])}, logical(1)))

  error <- .check_error_spec(error, sv, k)

  model <- NULL
  model[["type"]] <- ifelse(structural, "SVECX", ifelse(k == 1, "EC", "VEC"))
  model[["algorithm"]] <- .vec_algorithm(error, tvp)
  model[["k"]] <- k
  model[["p"]] <- p
  model[["m"]] <- m
  model[["s"]] <- s
  model[["n"]] <- n
  model[["n_restricted"]] <- n_restricted
  model[["rank"]] <- rank
  model[["k_beta"]] <- k_beta
  # Filled in below, once it is known whether inclusion parameters were provided
  model[["varsel"]] <- "none"
  model[["endogen"]] <- dimnames(data)[[2]]
  if (m > 0) {
    model[["exogen"]] <- dimnames(exogen)[[2]]
  }
  model[["structural"]] <- structural
  model[["error"]] <- error
  model[["tvp"]] <- tvp
  model[["iterations"]] <- if (is.null(iterations)) draws else as.integer(iterations)
  model[["burnin"]] <- as.integer(burnin)

  # Posterior draws ----

  blocks <- list(block_alpha, block_gamma, block_upsilon, block_c, block_a0)

  posterior <- NULL

  coeffs <- .bind_blocks(blocks, tt, tvp)
  if (!is.null(coeffs)) {
    posterior[["a"]][["coeffs"]] <- coda::mcmc(coeffs)

    # Coefficients without inclusion parameters are always part of the model
    lambda <- .bind_blocks(blocks, tt, tvp, "lambda", fill = 1)
    if (!is.null(lambda)) {
      posterior[["a"]][["lambda"]] <- coda::mcmc(lambda)
    }

    # Coefficients without a state equation do not vary over time
    sigma <- .bind_blocks(blocks, tt, tvp, "sigma", fill = 0)
    if (!is.null(sigma)) {
      posterior[["a"]][["sigma"]] <- coda::mcmc(sigma)
    }
  }

  if (rank > 0) {
    beta_blocks <- list(.split_draws(beta, "beta"),
                        .split_draws(beta_x, "beta_x"),
                        .split_draws(beta_d, "beta_d"))
    posterior[["beta"]][["coeffs"]] <- coda::mcmc(t(.bind_beta(beta_blocks,
                                                               c(k, m, n_restricted),
                                                               rank, k_beta)))
  }

  if (!is.null(block_sigma)) {
    posterior[["u_sigma_inv"]][["coeffs"]] <- coda::mcmc(.invert_sigma_draws(block_sigma[["coeffs"]], k))
    if (!is.null(block_sigma[["lambda"]])) {
      posterior[["psi"]][["lambda"]] <- coda::mcmc(t(block_sigma[["lambda"]]))
    }
  }

  if (is.null(varsel)) {
    varsel <- ifelse(is.null(posterior[["a"]][["lambda"]]) & is.null(posterior[["psi"]][["lambda"]]),
                     "none", "bvs")
  }
  model[["varsel"]] <- varsel

  # Result ----

  result <- list("model" = model,
                 "data" = list("original" = list("endogen" = data,
                                                 "exogen" = exogen),
                               "train" = list("y" = y,
                                              "w" = ect,
                                              "x" = regressors,
                                              "z" = z)))

  if (!is.null(posterior)) {
    result[["posterior"]] <- posterior
  }

  class(result) <- c("bvecmodel", "list")

  return(result)
}
