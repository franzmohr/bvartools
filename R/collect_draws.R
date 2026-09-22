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
# off and inverted. `need_Sigma` is there for irf(), which does not need the
# error covariance for a forecast error or structural response and would
# otherwise pay for one inversion per draw to obtain it. `impact` carries a
# caller supplied identification through to the workers, which read it under
# type "custom". Returns a list of lists with elements `A`, `Sigma` and, when
# asked for, `A0` and `P`.

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

# Per-draw impact matrices from whatever the caller supplied.
#
# One k x k matrix identifies every draw the same way, which is what a fixed
# identification looks like. A list of them is one matrix per draw, which is
# what an identification drawn alongside the coefficients produces. Both are
# normalised to a list of length `store` here, so that .collect_draws has a
# single case to attach and the workers never learn which of the two the caller
# had.
#
# A function is the third form and is left alone: it is called with the draw
# once the draw exists, because an identification can depend on quantities --
# the error covariance above all -- that are only assembled there. Returning
# NULL from it drops that draw, which is how a sign restricted identification
# reports that it found no admissible rotation.
.impact_draws <- function(impact, k, store) {

  if (is.function(impact)) {
    return(impact)
  }

  if (is.matrix(impact)) {
    impact <- rep(list(impact), store)
  }

  if (!is.list(impact)) {
    stop("Argument 'impact' must be a matrix, a list of matrices or a function.")
  }

  if (length(impact) != store) {
    stop("Argument 'impact' must contain either a single matrix or one per posterior draw (",
         store, "), but it contains ", length(impact), ".")
  }

  ok <- vapply(impact, .is_impact_matrix, logical(1), k = k)

  if (!all(ok)) {
    stop("Argument 'impact' must contain finite numeric ", k, " x ", k,
         " matrices. Element ", which(!ok)[1], " is not one.")
  }

  impact
}

.is_impact_matrix <- function(p_i, k) {
  is.matrix(p_i) && is.numeric(p_i) && nrow(p_i) == k && ncol(p_i) == k &&
    all(is.finite(p_i))
}

# The impact matrices of a sign restricted identification, in the form
# .collect_draws calls per draw.
#
# A closure rather than a list, because the impact matrix of a draw is the
# Choleski factor of that draw's own covariance times its own rotation, and the
# covariance is only assembled inside .collect_draws -- which is also the only
# place that knows which slice of a row belongs to `period`. Returning NULL for
# a draw that add_sign_restrictions() could not identify is what drops it from
# the sample.
.sign_impact <- function(x, caller) {

  rotations <- x[["posterior"]][["q"]][["coeffs"]]

  if (is.null(rotations)) {
    stop(caller, " of type \"sign\" need an identified model: run add_sign_restrictions() ",
         "on it first.")
  }

  k <- x[["model"]][["k"]]

  function(draw, i) {
    rotation <- rotations[i, ]
    if (anyNA(rotation)) {
      return(NULL)
    }
    t(chol(draw[["Sigma"]])) %*% matrix(rotation, k)
  }
}

# A quantile VAR has no error covariance to identify shocks from. What its
# posterior keeps in u_sigma_inv are the per-period precisions of the latent
# mixture that makes the asymmetric Laplace likelihood conditionally normal:
# diagonal, redrawn every period, and not a covariance of anything. Everything
# that factorises or decomposes Sigma is therefore refused, as forecasting and
# add_sign_restrictions() already refuse the model.
.refuse_quantile_covariance <- function(x, caller) {
  if (identical(x[["model"]][["error"]], "ald")) {
    stop(caller, " need the error covariance of the model, which a quantile VAR ",
         "(error = \"ald\") does not estimate: its u_sigma_inv holds the latent precisions ",
         "of the asymmetric Laplace errors, period by period.", call. = FALSE)
  }
}

# The period a sign restricted identification is used at.
#
# Each rotation was found against the covariance of its draw in the period
# add_sign_restrictions() was given, and satisfies the restrictions there. Used
# with the covariance of another period of a model whose covariance moves, it
# no longer does -- a share of the accepted draws then violate the very signs
# they were accepted for, and nothing says so. So the stored period is the
# default, and another one is refused where it would differ.
.sign_period <- function(x, period, caller) {

  stored <- x[["model"]][["sign_restrictions"]][["period"]]
  if (is.null(period)) {
    return(stored)
  }

  varies <- isTRUE(x[["model"]][["tvp"]]) ||
    isTRUE(x[["model"]][["error"]] %in% c("sv", "sv+covar"))
  if (!varies) {
    return(period)
  }

  tt <- .train_periods(x, x[["model"]][["k"]])
  if (.check_period(period, tt) != (if (is.null(stored)) tt else stored)) {
    stop(caller, " of type \"sign\" use the rotations add_sign_restrictions() found in period ",
         if (is.null(stored)) tt else stored, ", and the covariance of this model differs ",
         "from period to period, so they do not identify period ", period, ". Leave ",
         "'period' out, or run add_sign_restrictions() with that period.", call. = FALSE)
  }
  period
}

.collect_draws <- function(x, period = NULL, need_A0 = FALSE, need_Sigma = TRUE,
                           impact = NULL) {

  k <- x[["model"]][["k"]]
  kk <- k * k
  p <- x[["model"]][["p"]]
  tt <- .train_periods(x, k)
  tvp <- x[["model"]][["tvp"]]
  if (tvp) {
    nparams <- ncol(x[["data"]][["train"]][["z"]])
  }
  sv <- .error_varies_by_period(x[["model"]][["error"]])
  sigma_path <- .u_sigma_is_path(x, k, tt)
  if (tvp || sv || sigma_path) {
    if (is.null(period)) {
      period <- tt
    } else {
      period <- .check_period(period, tt)
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

    # The free elements are stored column by column -- (2,1), (3,1), ..., (k,1),
    # (3,2), ... -- which is the order of the contemporaneous regressors in z and
    # of which(lower.tri()). They used to be read row by row, which gives the
    # same matrix for up to three variables and swaps elements from four on.
    pos_a0 <- which(lower.tri(diag(k)))
  }

  store <- nrow(x[["posterior"]][["u_sigma_inv"]][["coeffs"]])

  if (!is.null(impact)) {
    impact <- .impact_draws(impact, k, store)
  }

  # Pre-allocated, because a draw that the identification declines is dropped
  # at the end rather than while the list is being built: assigning NULL into a
  # list removes the element and shifts every later one onto the wrong index.
  A <- vector("list", store)
  keep <- rep(TRUE, store)

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

    if (need_Sigma) {
      if (sigma_path) {
        temp[["Sigma"]] <- solve(matrix(x[["posterior"]][["u_sigma_inv"]][["coeffs"]][i, (period - 1) * kk + 1:kk], k))
      } else {
        temp[["Sigma"]] <- solve(matrix(x[["posterior"]][["u_sigma_inv"]][["coeffs"]][i, ], k))
      }
    }

    if (!is.null(impact)) {
      if (is.function(impact)) {
        p_i <- impact(temp, i)
        if (is.null(p_i)) {
          keep[i] <- FALSE
        } else {
          if (!.is_impact_matrix(p_i, k)) {
            stop("Argument 'impact' returned something other than a finite numeric ",
                 k, " x ", k, " matrix for draw ", i, ".")
          }
          temp[["P"]] <- p_i
        }
      } else {
        temp[["P"]] <- impact[[i]]
      }
    }

    A[[i]] <- temp
  }

  return(A[keep])
}
