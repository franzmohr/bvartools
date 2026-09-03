# Posterior draws of the coefficients and the error covariance, one list entry
# per draw, in the shape the C++ workers .vardecomp and .spillover_table expect.
#
# Lifted out of fevd.bvarmodel so that the variance decomposition and the
# spillover index cannot drift apart in how they index a draw. The awkward part
# is not the loop but the indexing that precedes it: a time varying model stores
# one coefficient vector per period end to end, a stochastic volatility model
# one covariance per period, and a structural model keeps its contemporaneous
# block at the end of `a` -- so which slice of a row belongs to `period` depends
# on three flags at once.
#
# `x` is a bvarmodel with posterior draws, `period` the period to read for a TVP
# or SV model, and `need_A0` whether the contemporaneous block has to be split
# off and inverted. Returns a list of lists with elements `A`, `Sigma` and, when
# asked for, `A0`.

# Number of periods in the estimation sample.
#
# Taken from the SUR regressor matrix where there is one, because its row count
# is unambiguously `k` times the sample length. `y` is only fallen back on when
# there are no regressors, and then its layout has to be decided: a single
# column holds the stacked series, anything wider is one row per period. The
# stacked reading used to be the only one -- it is what gen_var produced -- and
# create_bvarmodel returns the wide one, so both have to be handled.
.train_periods <- function(x, k) {
  z <- x[["data"]][["train"]][["z"]]
  if (!is.null(z)) {
    return(as.integer(round(NROW(z) / k)))
  }
  y <- x[["data"]][["train"]][["y"]]
  if (NCOL(y) == 1) {
    return(as.integer(round(NROW(y) / k)))
  }
  as.integer(NROW(y))
}

.collect_draws <- function(x, period = NULL, need_A0 = FALSE) {

  k <- x[["model"]][["k"]]
  kk <- k * k
  p <- x[["model"]][["p"]]
  tt <- .train_periods(x, k)
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

    pos_a0 <- t(matrix(1:kk, k, k))
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

  return(A)
}
