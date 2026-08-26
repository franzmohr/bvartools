#' Plotting Draws of a Bayesian VEC Model
#'
#' A plot function for objects of class 'bvecmodel'.
#'
#' @param x an object of class 'bvecmodel'.
#' @param ci interval used to calculate credible bands for time-varying parameters.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param show_zero_y if \code{TRUE} (default), a horizontal line with y = 0 is
#' added to the plot. Only used for time varying parameters.
#' @param ... further graphical parameters.
#'
#' @details The coefficients of the error correction term are displayed as draws
#' of the cointegration matrix \eqn{\Pi = \alpha \beta^\prime} and not as draws
#' of the loading matrix \eqn{\alpha} and the cointegration matrix \eqn{\beta}
#' separately. The latter two are only identified up to a rotation, so their
#' individual draws are not informative, while their product is.
#'
#' @examples
#'
#' # Load data
#' data("e6")
#' e6 <- e6 * 100
#'
#' # Create model
#' model <- create_bvecmodel(e6, p = 2, r = 1, const = "restricted",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#'
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 1))
#'
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws
#' model <- add_posterior_coefficients(model)
#'
#' # Plot
#' plot(model, type = "hist")
#' plot(model, type = "trace")
#' plot(model, type = "boxplot")
#'
#' @export
plot.bvecmodel <- function(x, ci = 0.95, type = "hist", show_zero_y = TRUE, ...) {

  # 'layout' is called below, so all parameters have to be restored on exit
  orig_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(orig_par))

  if (!type %in% c("hist", "trace", "boxplot")) {
    stop("Argument 'type' must be 'hist', 'trace' or 'boxplot'.")
  }

  if (is.null(x[["posterior"]])) {
    stop("Argument 'x' does not contain posterior draws.")
  }

  k <- x[["model"]][["k"]]
  kk <- k * k
  p <- x[["model"]][["p"]]
  m <- x[["model"]][["m"]]
  s <- x[["model"]][["s"]]
  rank <- x[["model"]][["rank"]]
  k_beta <- k + m + x[["model"]][["n_restricted"]]
  tvp <- x[["model"]][["tvp"]]
  sv <- x[["model"]][["error"]] %in% c("sv", "sv+covar")
  tvp_and_covar <- tvp & x[["model"]][["error"]] == "gamma+covar"
  structural <- x[["model"]][["structural"]]
  if (structural) {
    n_struct <- k * (k - 1) / 2
  } else {
    n_struct <- 0
  }

  tt <- nrow(x[["data"]][["train"]][["y"]])
  draws <- nrow(x[["posterior"]][["a"]][["coeffs"]])

  # Blocks of coefficients. The error correction term contributes the columns of
  # Pi, which are those of the regressors in 'w', and not those of alpha.
  n_alpha <- k * rank
  n_x <- ifelse(is.null(x[["data"]][["train"]][["x"]]), 0, NCOL(x[["data"]][["train"]][["x"]]))
  n_pi <- ifelse(rank > 0, k * k_beta, 0)
  ncoeffs <- n_alpha + k * n_x + n_struct

  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  y_names <- x[["model"]][["endogen"]]
  x_names <- .get_regressor_names_bvecmodel(x, add_block = TRUE)
  lab_size <- .05

  nparams <- ifelse(rank > 0, k_beta, 0) + n_x + ifelse(structural, k, 0) + k

  mat <- matrix(NA_integer_, k + 2 , nparams + 1)
  mat[1, ] <- 1
  mat[-1, 1] <- c(0, 2:(k + 1))
  mat[2, -1] <- (k + 1) + 1:nparams
  mat[-(1:2), -1] <- matrix(1:(k * nparams) + k + nparams + 1, k, nparams)
  graphics::layout(mat,
                   widths = c(lab_size, rep((1 - lab_size) / nparams, nparams)),
                   heights = c(.07, lab_size, rep((1 - lab_size) / k, k)))

  # Title
  title_text <- "Bayesian "
  if (tvp) {
    title_text <- paste0(title_text, "TVP-")
  }
  if (sv) {
    title_text <- paste0(title_text, "SV-")
  }
  if (structural) {
    title_text <- paste0(title_text, "S")
  }
  title_text <- paste0(title_text, "VEC model")
  spec_text <- paste0("p = ", p)
  if (m > 0) {
    spec_text <- c(spec_text, paste0("s = ", s))
  }
  spec_text <- c(spec_text, paste0("r = ", rank))
  title_text <- paste0(title_text, " with ", paste0(spec_text, collapse = ", "))

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new(); graphics::text(0.5, 0.5, labels = title_text, cex = 1.5)
  # Fill rows
  graphics::par(mar = c(3, 0, 0, 0))
  for (j in y_names) {
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = j, adj = 0.5)
  }
  # Fill columns
  graphics::par(mar = c(0, 0, 0, 0))
  for (j in x_names) {
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = j, adj = 0.5)
  }
  for (j in y_names) {
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = paste0("Sigma\n", j), adj = 0.5)
  }

  graphics::par(mar = c(3, 2.1, .5, 1))

  periods <- ifelse(tvp, tt, 1)

  # Draws of a coefficient of the measurement equation, either as a single
  # vector or, for time varying parameters, as one column per period
  coeff_draws <- function(pos) {
    x[["posterior"]][["a"]][["coeffs"]][, ncoeffs * 0:(periods - 1) + pos, drop = !tvp]
  }

  # Cointegration matrix ----

  if (rank > 0) {

    n_beta <- k_beta * rank
    beta_draws <- x[["posterior"]][["beta"]][["coeffs"]]
    # Draws of beta are time varying only if they cover every period
    beta_tvp <- tvp & NCOL(beta_draws) == tt * n_beta

    Pi <- matrix(NA_real_, draws, periods * n_pi)
    for (i in 1:periods) {
      alpha_i <- x[["posterior"]][["a"]][["coeffs"]][, (i - 1) * ncoeffs + 1:n_alpha, drop = FALSE]
      if (beta_tvp) {
        beta_i <- beta_draws[, (i - 1) * n_beta + 1:n_beta, drop = FALSE]
      } else {
        beta_i <- beta_draws[, 1:n_beta, drop = FALSE]
      }
      for (j in 1:draws) {
        Pi[j, (i - 1) * n_pi + 1:n_pi] <- tcrossprod(matrix(alpha_i[j, ], k),
                                                     matrix(beta_i[j, ], ncol = rank))
      }
    }

    for (i in 1:n_pi) {
      if (tvp) {
        temp <- Pi[, n_pi * 0:(tt - 1) + i, drop = FALSE]
        stats::ts.plot(t(apply(temp, 2, stats::quantile, probs = c(ci_low, .5, ci_high))), xlab = "")
        if (show_zero_y) {
          graphics::abline(h = 0)
        }
      } else {
        if (type == "hist") {
          graphics::hist(Pi[, i], plot = TRUE, main = NA)
        }
        if (type == "trace") {
          stats::ts.plot(Pi[, i], xlab = "")
        }
        if (type == "boxplot") {
          graphics::boxplot(Pi[, i])
        }
      }
    }
  }

  # Coefficients of the regressors outside the error correction term ----

  n_nonect <- k * n_x
  if (n_nonect > 0) {
    for (i in 1:n_nonect) {
      if (tvp) {
        temp <- coeff_draws(n_alpha + i)
        stats::ts.plot(t(apply(temp, 2, stats::quantile, probs = c(ci_low, .5, ci_high))), xlab = "")
        if (show_zero_y) {
          graphics::abline(h = 0)
        }
      } else {
        if (type == "hist") {
          graphics::hist(coeff_draws(n_alpha + i), plot = TRUE, main = NA)
        }
        if (type == "trace") {
          stats::ts.plot(coeff_draws(n_alpha + i), xlab = "")
        }
        if (type == "boxplot") {
          graphics::boxplot(coeff_draws(n_alpha + i))
        }
      }
    }
  }

  # Structural coefficients ----

  if (structural) {

    struct_matrix <- matrix(1:kk, k)
    pos_values <- which(lower.tri(struct_matrix))
    pos_zero <- which(upper.tri(struct_matrix))
    pos_one <- struct_matrix[-c(pos_values, pos_zero)]

    # Position of the free elements of A0 in the vector of coefficients
    temp <- matrix(NA, k , k)
    temp[upper.tri(temp)] <- 1:n_struct
    temp <- t(temp)
    pos_a <- n_alpha + n_nonect + temp[lower.tri(temp)]

    pos_i <- 0
    for (i in 1:kk) {
      if (i %in% pos_values) {
        pos_i <- pos_i + 1

        if (tvp) {
          temp <- coeff_draws(pos_a[pos_i])
          stats::ts.plot(t(apply(temp, 2, stats::quantile, probs = c(ci_low, .5, ci_high))), xlab = "")
          if (show_zero_y) {
            graphics::abline(h = 0)
          }
        } else {
          if (type == "hist") {
            graphics::hist(coeff_draws(pos_a[pos_i]), plot = TRUE, main = NA)
          }
          if (type == "trace") {
            stats::ts.plot(coeff_draws(pos_a[pos_i]), xlab = "")
          }
          if (type == "boxplot") {
            graphics::boxplot(coeff_draws(pos_a[pos_i]))
          }
        }
      } else {
        if (i %in% pos_zero) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = 0, adj = 0.5)
        }
        if (i %in% pos_one) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = 1, adj = 0.5)
        }
      }
    }
  }

  # Covariance matrix of the error term ----

  # Obtain inverse and calculate bands
  if (sv | tvp_and_covar) {

    if (k == 1) {
      temp <- matrix(1 / x[["posterior"]][["u_sigma_inv"]][["coeffs"]], ncol = tt)
    } else {
      temp <- x[["posterior"]][["u_sigma_inv"]][["coeffs"]]
      for (i in 1:tt) {
        temp[, (i - 1) * kk + 1:kk] <- t(apply(x[["posterior"]][["u_sigma_inv"]][["coeffs"]][, (i - 1) * kk + 1:kk], 1, function(x, k) {solve(matrix(x, k))}, k = k))
      }
    }
    u_sigma <- t(apply(temp, 2, stats::quantile, probs = c(ci_low, .5, ci_high)))
  } else {
    if (k == 1) {
      u_sigma <- matrix(1 / x[["posterior"]][["u_sigma_inv"]][["coeffs"]])
    } else {
      u_sigma <- t(apply(x[["posterior"]][["u_sigma_inv"]][["coeffs"]], 1, function(x, k) {solve(matrix(x, k))}, k = k))
    }
  }

  for (i in 1:kk) {
    if (sv | tvp_and_covar) {
      pos <- kk * 0:(tt - 1) + i
      stats::plot.ts(u_sigma[pos, ], plot.type = "single")
    } else {
      if (all(u_sigma[, i] == u_sigma[1, i])) {
        graphics::plot.new(); graphics::text(0.5, 0.5, labels = u_sigma[1, i], adj = 0.5)
      } else {
        if (type == "hist") {
          graphics::hist(u_sigma[, i], plot = TRUE, main = NA)
        }
        if (type == "trace") {
          stats::ts.plot(u_sigma[, i], xlab = "")
        }
        if (type == "boxplot") {
          graphics::boxplot(u_sigma[, i])
        }
      }
    }
  }
}
