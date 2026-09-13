# Helpers of generate_artificial_var and generate_artificial_vec

# Checks that 'x' is a single whole number of at least 'min'
.artificial_count <- function(x, name, min) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < min || x != round(x)) {
    stop(paste0("Argument '", name, "' must be a whole number of at least ", min, "."))
  }
}

# Checks that the named arguments in list 'x' are TRUE or FALSE
.artificial_flags <- function(x) {
  for (i in names(x)) {
    if (!isTRUE(x[[i]]) & !isFALSE(x[[i]])) {
      stop(paste0("Argument '", i, "' must be TRUE or FALSE."))
    }
  }
}

# Checks that 'x' is a share between 0 and 1
.artificial_share <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0 || x > 1) {
    stop(paste0("Argument '", name, "' must be between 0 and 1."))
  }
}

# Checks that 'x' contains two numbers and returns them in increasing order
.artificial_range <- function(x, name, positive = FALSE, non_negative = FALSE) {
  if (!is.numeric(x) || length(x) != 2 || any(!is.finite(x))) {
    stop(paste0("Argument '", name, "' must be a numeric vector with two finite elements."))
  }
  if (positive && any(x <= 0)) {
    stop(paste0("The elements of argument '", name, "' must be positive."))
  }
  if (non_negative && any(x < 0)) {
    stop(paste0("The elements of argument '", name, "' must not be negative."))
  }
  return(c(min(x), max(x)))
}

# Draws 'n' uniform values from 'range' and sets a share 'zeros' of them to zero
.artificial_draw <- function(n, range, zeros = 0) {
  x <- stats::runif(n, range[1], range[2])
  x[stats::runif(n) < zeros] <- 0
  return(x)
}

# Moduli of the eigenvalues of the companion matrix of the K x Kp lag coefficients 'a'
.artificial_roots <- function(a, k, p) {
  companion <- matrix(0, k * p, k * p)
  companion[1:k, ] <- a
  if (p > 1) {
    companion[k + 1:(k * (p - 1)), 1:(k * (p - 1))] <- diag(1, k * (p - 1))
  }
  return(Mod(eigen(companion, only.values = TRUE)$values))
}

# Checks whether the VAR process with the K x Kp lag coefficients 'a' is stable
.artificial_var_stable <- function(a, k, p) {
  return(all(.artificial_roots(a, k, p) < 1))
}

# Checks whether the VEC process with loadings 'alpha', cointegration vectors 'beta',
# whose first K rows belong to the endogenous variables, and the K x K(p - 1) matrix
# 'gamma' is integrated of order one with cointegration rank 'r', i.e. whether its
# VAR representation has exactly K - r unit roots and all other roots inside the
# unit circle
.artificial_vec_stable <- function(alpha, beta, gamma, k, p, r) {
  a <- matrix(0, k, k * p)
  a[, 1:k] <- diag(1, k)
  if (r > 0) {
    a[, 1:k] <- a[, 1:k] + alpha %*% t(beta[1:k, , drop = FALSE])
  }
  if (p > 1) {
    for (i in 1:(p - 1)) {
      gamma_i <- gamma[, (i - 1) * k + 1:k, drop = FALSE]
      a[, (i - 1) * k + 1:k] <- a[, (i - 1) * k + 1:k] + gamma_i
      a[, i * k + 1:k] <- a[, i * k + 1:k] - gamma_i
    }
  }
  roots <- .artificial_roots(a, k, p)
  tol <- 1e-6
  return(sum(roots > 1 - tol) == k - r && all(roots < 1 + tol))
}

# State variances drawn from 'range' for the elements of matrix 'x'. Elements that
# are zero get a state variance of zero, so that they stay zero.
.artificial_state_variance <- function(x, range) {
  result <- matrix(stats::runif(length(x), range[1], range[2]), nrow(x), ncol(x),
                   dimnames = dimnames(x))
  return(result * (x != 0))
}

# One step of the random walks of the coefficient matrices in list 'current' with
# the state variances in list 'variance'. Innovations, for which 'is_stable' is
# FALSE, are drawn again. After 100 failed attempts the coefficients are kept.
.artificial_tvp_step <- function(current, variance, is_stable) {
  for (i in 1:100) {
    candidate <- Map(function(x, v) {x + stats::rnorm(length(x)) * sqrt(v)},
                     current, variance[names(current)])
    if (is_stable(candidate)) {
      return(candidate)
    }
  }
  return(current)
}

# Array of the matrices named 'name' in the list of periods 'hist', with the
# dimensions and names of matrix 'template'
.artificial_path <- function(hist, name, template) {
  dim_names <- NULL
  if (!is.null(dimnames(template))) {
    dim_names <- c(dimnames(template), list(NULL))
  }
  x <- unlist(lapply(hist, function(x) {x[[name]]}), use.names = FALSE)
  return(array(x, c(dim(template), length(hist)), dimnames = dim_names))
}

# Parameters of the error term in the first period and the state variances
.artificial_error_init <- function(k, names_series, tvp, sv, range_variance, range_psi,
                                   range_variance_state, range_variance_sv) {

  omega <- stats::runif(k, range_variance[1], range_variance[2])
  psi <- diag(1, k)
  if (k > 1) {
    psi[lower.tri(psi)] <- stats::runif(k * (k - 1) / 2, range_psi[1], range_psi[2])
  }
  dimnames(psi) <- list(names_series, names_series)

  err <- list(k = k, tvp = tvp, sv = sv, omega = omega, psi = psi,
              psi_t = psi, log_omega_t = log(omega))

  if (tvp) {
    psi_state <- matrix(0, k, k, dimnames = dimnames(psi))
    if (k > 1) {
      pos_psi <- lower.tri(psi) & psi != 0
      psi_state[pos_psi] <- stats::runif(sum(pos_psi), range_variance_state[1], range_variance_state[2])
    }
    err[["psi_state"]] <- psi_state
  }
  if (sv) {
    err[["sv_state"]] <- stats::runif(k, range_variance_sv[1], range_variance_sv[2])
    names(err[["sv_state"]]) <- names_series
  }

  return(err)
}

# One step of the random walks of psi and the log-volatilities
.artificial_error_update <- function(err) {
  k <- err[["k"]]
  if (err[["tvp"]]) {
    err[["psi_t"]] <- err[["psi_t"]] + stats::rnorm(k * k) * sqrt(err[["psi_state"]])
  }
  if (err[["sv"]]) {
    err[["log_omega_t"]] <- err[["log_omega_t"]] + stats::rnorm(k) * sqrt(err[["sv_state"]])
  }
  return(err)
}

# Draws an error of the current period from N(0, Psi^-1 Omega Psi^-1')
.artificial_error_draw <- function(err) {
  u <- forwardsolve(err[["psi_t"]], sqrt(exp(err[["log_omega_t"]])) * stats::rnorm(err[["k"]]))
  return(as.numeric(u))
}

# True parameters of the error term from the first period in 'err' and, for time
# varying models, the list of periods 'hist'
.artificial_error_params <- function(err, hist) {

  k <- err[["k"]]
  tvp <- err[["tvp"]]
  sv <- err[["sv"]]

  psi <- err[["psi"]]
  psi_inv <- forwardsolve(psi, diag(1, k))
  omega <- diag(err[["omega"]], k)
  dimnames(omega) <- dimnames(psi)
  sigma <- psi_inv %*% omega %*% t(psi_inv)
  dimnames(sigma) <- dimnames(psi)

  if (tvp | sv) {
    nobs <- length(hist)
    dim_names <- c(dimnames(psi), list(NULL))
    omega_path <- array(0, c(k, k, nobs), dimnames = dim_names)
    sigma_path <- array(NA_real_, c(k, k, nobs), dimnames = dim_names)
    for (j in 1:nobs) {
      omega_t <- diag(exp(hist[[j]][["log_omega_t"]]), k)
      psi_inv <- forwardsolve(hist[[j]][["psi_t"]], diag(1, k))
      omega_path[, , j] <- omega_t
      sigma_path[, , j] <- psi_inv %*% omega_t %*% t(psi_inv)
    }
  }

  params <- list(psi_coef = if (tvp) .artificial_path(hist, "psi_t", psi) else psi,
                 u_omega = if (sv) omega_path else omega,
                 u_sigma = if (tvp | sv) sigma_path else sigma)
  if (k == 1) {
    params["psi_coef"] <- list(NULL)
  }

  states <- list()
  if (tvp & k > 1) {
    states[["psi_state_variance"]] <- err[["psi_state"]]
  }
  if (sv) {
    states[["u_state_variance"]] <- err[["sv_state"]]
  }

  return(list(params = params, states = states))
}
