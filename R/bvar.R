#' Bayesian Vector Autoregression Objects
#'
#' \code{bvar} collects the posterior draws of a vector autoregressive model in
#' an object of class 'bvarmodel'.
#'
#' @param data the original time-series object of endogenous variables. If
#' \code{NULL} (default), the object provided in argument \code{y} is used.
#' @param exogen the original time-series object of unmodelled variables. If
#' \code{NULL} (default), it is reconstructed from the current values of those
#' variables in argument \code{x}.
#' @param y a time-series object of endogenous variables with \eqn{T} observations,
#' usually, a result of a call to \code{\link{create_bvarmodel}}.
#' @param x a time-series object of \eqn{(pK + (1+s)M + N)} regressor variables, usually, a result of a
#' call to \code{\link{create_bvarmodel}}.
#' @param z a \eqn{TK \times K (pK + (1+s)M + N + K)} data matrix, usually, a result of a
#' call to \code{\link{create_bvarmodel}}. If \code{NULL} (default), it is
#' generated from argument \code{x}.
#' @param A0 either a \eqn{K^2 \times S} matrix of MCMC coefficient draws of structural parameters or
#' a named list, where element \code{coeffs} contains a \eqn{K^2 \times S} matrix of MCMC coefficient
#' draws of structural parameters and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed. For time varying parameter models
#' the coefficient matrix must be \eqn{TK^2 \times S}. Draws of the error covariance matrix of the state
#' equation can be provided as a \eqn{K^2 \times S} matrix in an additional list element.
#' @param A either a \eqn{pK^2 \times S} matrix of MCMC coefficient draws of lagged endogenous variables or
#' a named list, where element \code{coeffs} contains a \eqn{pK^2 \times S} matrix of MCMC coefficient draws
#' of lagged endogenous variables and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed. For time varying parameter models
#' the coefficient matrix must be \eqn{pTK^2 \times S}. Draws of the error covariance matrix of the state
#' equation can be provided as a \eqn{pK^2 \times S} matrix in an additional list element.
#' @param B either a \eqn{((1 + s)MK) \times S} matrix of MCMC coefficient draws of unmodelled, non-deterministic variables
#' or a named list, where element \code{coeffs} contains a \eqn{((1 + s)MK) \times S} matrix of MCMC coefficient draws of
#' unmodelled, non-deterministic variables and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed. For time varying parameter models
#' the coefficient matrix must be \eqn{(1 + s)TMK \times S}. Draws of the error covariance matrix of the state
#' equation can be provided as a \eqn{(1 + s)MK \times S} matrix in an additional list element.
#' @param C either a \eqn{KN \times S} matrix of MCMC coefficient draws of deterministic terms or
#' a named list, where element \code{coeffs} contains a \eqn{KN \times S} matrix of MCMC coefficient draws of
#' deterministic terms and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed. For time varying parameter models
#' the coefficient matrix must be \eqn{TKN \times S}. Draws of the error covariance matrix of the state
#' equation can be provided as a \eqn{KN \times S} matrix in an additional list element.
#' @param Sigma a \eqn{K^2 \times S} matrix of MCMC draws for the error variance-covariance matrix or
#' a named list, where element \code{coeffs} contains a \eqn{K^2 \times S} matrix of MCMC draws for the
#' error variance-covariance matrix and element \code{lambda} contains the corresponding draws of inclusion
#' parameters in case variable selection algorithms were employed to the covariances. For models with
#' stochastic volatility the matrix must be \eqn{TK^2 \times S}.
#' @param error a character specifying the model that was used for the estimation
#' of the covariance matrix of the error term. If \code{NULL} (default), it is
#' inferred from the draws in argument \code{Sigma}, which can only distinguish
#' a constant covariance matrix from a time varying one. See
#' \code{\link{create_bvarmodel}}.
#' @param varsel a character specifying the variable selection algorithm that was
#' employed. If \code{NULL} (default), it is inferred from the presence of draws
#' of inclusion parameters. See \code{\link{create_bvarmodel}}.
#' @param iterations an integer of the number of MCMC draws in the provided
#' coefficient matrices. If \code{NULL} (default), it is obtained from those
#' matrices.
#' @param burnin an integer of the number of MCMC draws that were used to
#' initialise the sampler. Defaults to zero, since the draws provided to the
#' function are assumed to exclude them already.
#'
#' @details For the VARX model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_i y_{t-i} + \sum_{i = 0}^{s} B_i x_{t - i} + C d_t + u_t}
#' the function collects the S draws of a Gibbs sampler (after the burn-in phase) in a standardised object,
#' where \eqn{y_t} is a K-dimensional vector of endogenous variables,
#' \eqn{A_0} is a \eqn{K \times K} matrix of structural coefficients.
#' \eqn{A_i} is a \eqn{K \times K} coefficient matrix of lagged endogenous variabels.
#' \eqn{x_t} is an M-dimensional vector of unmodelled, non-deterministic variables
#' and \eqn{B_i} its corresponding coefficient matrix.
#' \eqn{d_t} is an N-dimensional vector of deterministic terms
#' and \eqn{C} its corresponding coefficient matrix.
#' \eqn{u_t} is an error term with \eqn{u_t \sim N(0, \Sigma_u)}.
#'
#' For time varying parameter and stochastic volatility models the respective coefficients and
#' error covariance matrix of the above model are assumed to be time varying, respectively.
#'
#' The draws of the different coefficient matrices provided in \code{A0}, \code{A},
#' \code{B}, \code{C} and \code{Sigma} have to correspond to the same MCMC iterations.
#'
#' The result is the same kind of object as the output of
#' \code{\link{create_bvarmodel}} in combination with
#' \code{\link{draw_posterior}}, so it can be used with the methods of class
#' 'bvarmodel' such as \code{\link{irf.bvarmodel}},
#' \code{\link{predict.bvarmodel}} or \code{\link{summary.bvarmodel}}.
#' Accordingly, the draws are stored in the parameterisation those methods
#' expect: the coefficients of the conditional mean are collected in a single
#' matrix \code{a} in the order in which their regressors enter the data matrix
#' and the draws of the error term are stored as draws of the inverse of
#' \eqn{\Sigma_u}. Draws of \eqn{A_0} are reduced to the free elements below its
#' diagonal, which are the elements that are estimated.
#'
#' Since a model that is constructed from posterior draws does not contain the
#' specification of its priors, the elements \code{priors} and \code{initial} of
#' the resulting object are empty. Use \code{\link{add_priors}} and
#' \code{\link{add_initial_values}} if the model should be re-estimated.
#'
#' @return An object of class \code{"bvarmodel"} containing the following components:
#' \item{model}{a list containing information on the model specification.}
#' \item{data}{a list of data objects. Element \code{original} contains the
#' original time-series objects of endogenous, unmodelled and deterministic
#' variables. Element \code{train} contains the time-series objects \code{y} and
#' \code{x} of dependent variables and regressors as well as the matrix \code{z}
#' of regressors in SUR form.}
#' \item{posterior}{a list of posterior draws. Element \code{a} contains an
#' \eqn{S \times (K(pK + (1 + s)M + N) + K(K - 1) / 2)} "mcmc" object of the
#' draws of the coefficients of the conditional mean and, if provided, the
#' corresponding draws of inclusion parameters in element \code{lambda} and of
#' the error covariance matrix of the state equation in element \code{sigma}.
#' Element \code{u_sigma_inv} contains an \eqn{S \times K^2} "mcmc" object of
#' the draws of the inverse of the error variance-covariance matrix. For time
#' varying parameter and stochastic volatility models the draws of all periods
#' are appended to each other.}

#' @examples
#'
#' # Get data
#' data("e1")
#' e1 <- diff(log(e1))
#' e1 <- window(e1, end = c(1978, 4))
#'
#' # Generate model data
#' data <- create_bvarmodel(e1, p = 2, deterministic = "const")
#'
#' # Add priors
#' model <- add_priors(data,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     sigma = list(df = 0, scale = .00001))
#'
#' # Set RNG seed for reproducibility
#' set.seed(1234567)
#'
#' iterations <- 400 # Number of iterations of the Gibbs sampler
#' # Chosen number of iterations and burnin should be much higher.
#' burnin <- 100 # Number of burn-in draws
#' draws <- iterations + burnin # Total number of MCMC draws
#'
#' y <- t(model$data$train$y)
#' x <- t(model$data$train$x)
#' tt <- ncol(y) # Number of observations
#' k <- nrow(y) # Number of endogenous variables
#' m <- k * nrow(x) # Number of estimated coefficients
#'
#' # Priors
#' a_mu_prior <- model$priors$a$mu # Vector of prior parameter means
#' a_v_i_prior <- model$priors$a$v_inv # Inverse of the prior covariance matrix
#'
#' u_sigma_df_prior <- model$priors$u_sigma$df # Prior degrees of freedom
#' u_sigma_scale_prior <- model$priors$u_sigma$scale # Prior covariance matrix
#' u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom
#'
#' # Initial values
#' u_sigma_i <- diag(1 / .00001, k)
#'
#' # Data containers for posterior draws
#' draws_a <- matrix(NA, m, iterations)
#' draws_sigma <- matrix(NA, k^2, iterations)
#'
#' # Start Gibbs sampler
#' for (draw in 1:draws) {
#'  # Draw conditional mean parameters
#'  a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior)
#'
#'  # Draw variance-covariance matrix
#'  u <- y - matrix(a, k) %*% x # Obtain residuals
#'  u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
#'  u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
#'
#'  # Store draws
#'  if (draw > burnin) {
#'   draws_a[, draw - burnin] <- a
#'   draws_sigma[, draw - burnin] <- solve(u_sigma_i)
#'  }
#' }
#'
#' # Generate bvarmodel object
#' bvar_est <- bvar(data = e1,
#'                  y = model$data$train$y, x = model$data$train$x,
#'                  A = draws_a[1:18,], C = draws_a[19:21, ],
#'                  Sigma = draws_sigma)
#'
#' @export
bvar <- function(data = NULL, exogen = NULL, y, x = NULL, z = NULL,
                 A0 = NULL, A = NULL, B = NULL,
                 C = NULL, Sigma = NULL,
                 error = NULL, varsel = NULL, iterations = NULL, burnin = 0) {

  if (!"ts" %in% class(y)) {
    stop("Argument 'y' must be an object of class time-series")
  }
  if (!is.null(x)) {
    if (!"ts" %in% class(x)) {
      stop("Argument 'x' must be an object of class time-series.")
    }
  }

  y <- .name_ts(y, "y")
  k <- NCOL(y)
  tt <- NROW(y)

  # Coefficient draws ----

  block_a0 <- .split_draws(A0, "A0")
  block_a <- .split_draws(A, "A")
  block_b <- .split_draws(B, "B")
  block_c <- .split_draws(C, "C")
  block_sigma <- .split_draws(Sigma, "Sigma")

  structural <- !is.null(block_a0)
  if (structural & k == 1) {
    stop("Argument 'A0' cannot be used for models with a single endogenous variable.")
  }
  if (!is.null(block_b) & is.null(x) & is.null(exogen)) {
    stop("Please specify either argument 'x' or 'exogen' when using exogenous regressors.")
  }

  block_a0 <- .block_dims(block_a0, tt, k * k, "A0")
  block_a <- .block_dims(block_a, tt, k * k, "A")
  block_b <- .block_dims(block_b, tt, k, "B")
  block_c <- .block_dims(block_c, tt, k, "C")
  block_sigma <- .block_dims(block_sigma, tt, k * k, "Sigma")

  # Only the elements below the diagonal of A0 are estimated
  if (structural) {
    block_a0 <- .a0_free_elements(block_a0, k)
    block_a0[["n"]] <- k * (k - 1) / 2
  }

  draws <- .count_draws(list(block_a0, block_a, block_b, block_c, block_sigma))

  # Dimensions of the model ----

  p <- if (is.null(block_a)) 0L else as.integer(block_a[["n"]] / (k * k))
  n <- if (is.null(block_c)) 0L else as.integer(block_c[["n"]] / k)

  m <- 0L
  s <- 0L
  exogen_names <- NULL
  exogen_pos <- NULL
  if (!is.null(block_b)) {
    if (is.null(x) | is.null(dimnames(x)[[2]])) {
      if (is.null(exogen)) {
        stop("Argument 'exogen' is required if the columns of argument 'x' are not named.")
      }
      exogen <- .name_ts(exogen, "exogen")
      m <- NCOL(exogen)
      exogen_names <- dimnames(exogen)[[2]]
      s <- block_b[["n"]] / (k * m) - 1
    } else {
      # The lag of an unmodelled variable is the trailing number of its column
      # name, as generated by create_bvarmodel
      pos <- k * p + 1:(block_b[["n"]] / k)
      x_names <- dimnames(x)[[2]][pos]
      x_lags <- suppressWarnings(as.integer(sub("^.*\\.l", "", x_names)))
      if (any(is.na(x_lags))) {
        stop("The lag order of the unmodelled variables cannot be obtained from the column names of argument 'x'. Please specify argument 'exogen' instead.")
      }
      exogen_names <- unique(sub("\\.l[0-9]+$", "", x_names))
      exogen_pos <- pos[x_lags == 0]
      m <- length(exogen_pos)
      s <- max(x_lags)
    }
  }
  m <- as.integer(m)
  s <- as.integer(s)

  # Data ----

  if (is.null(data)) {
    data <- y
  }
  data <- .name_ts(data, "y")

  # The current values of the unmodelled variables are part of the regressors,
  # so the original series can be recovered from them
  if (is.null(exogen) & !is.null(exogen_pos)) {
    exogen <- .subset_ts(x, exogen_pos, exogen_names)
  }
  exogen <- .name_ts(exogen, "exogen")

  det_data <- NULL
  if (n > 0 & !is.null(x)) {
    det_data <- .subset_ts(x, NCOL(x) - n + 1:n)
  }

  if (is.null(z) & !is.null(x)) {
    z <- kronecker(x, diag(1, k))
    if (structural) {
      y_A0 <- kronecker(-y, diag(1, k))
      pos <- NULL
      for (i in 1:k) {
        pos <- c(pos, (i - 1) * k + 1:i)
      }
      z <- cbind(z, y_A0[, -pos])
    }
    dimnames(z) <- NULL
  }

  # Model specification ----

  sv <- isTRUE(block_sigma[["tvp"]])
  tvp <- any(vapply(list(block_a0, block_a, block_b, block_c),
                    function(x) {isTRUE(x[["tvp"]])}, logical(1)))

  error <- .check_error_spec(error, sv, k)

  model <- NULL
  model[["type"]] <- ifelse(structural, "SVAR", ifelse(k == 1, "AR", "VAR"))
  model[["algorithm"]] <- .var_algorithm(error, tvp)
  model[["k"]] <- k
  model[["p"]] <- p
  model[["m"]] <- m
  model[["s"]] <- s
  model[["n"]] <- n
  # Filled in below, once it is known whether inclusion parameters were provided
  model[["varsel"]] <- "none"
  model[["endogen"]] <- dimnames(data)[[2]]
  if (m > 0) {
    model[["exogen"]] <- exogen_names
  }
  model[["structural"]] <- structural
  model[["error"]] <- error
  model[["tvp"]] <- tvp
  model[["iterations"]] <- if (is.null(iterations)) draws else as.integer(iterations)
  model[["burnin"]] <- as.integer(burnin)

  # Posterior draws ----

  blocks <- list(block_a, block_b, block_c, block_a0)

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
                                                 "exogen" = exogen,
                                                 "deterministic" = det_data),
                               "train" = list("y" = y,
                                              "x" = x,
                                              "z" = z)))

  if (!is.null(posterior)) {
    result[["posterior"]] <- posterior
  }

  class(result) <- c("bvarmodel", "list")

  return(result)
}
