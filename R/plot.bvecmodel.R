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
#' @param max_cols an integer of the maximum number of regressors per figure. A block
#' with more regressors than this is drawn as several figures of nearly equal width.
#' Defaults to 6.
#' @param ... further graphical parameters.
#'
#' @details The coefficients of the error correction term are displayed as draws
#' of the cointegration matrix \eqn{\Pi = \alpha \beta^\prime} and not as draws
#' of the loading matrix \eqn{\alpha} and the cointegration matrix \eqn{\beta}
#' separately. The latter two are only identified up to a rotation, so their
#' individual draws are not informative, while their product is.
#'
#' The function draws one figure per block of coefficients -- the cointegration
#' matrix, the lagged differenced endogenous variables, the differenced exogenous
#' variables, the deterministic terms, the contemporaneous endogenous variables of
#' a structural model, and the covariance matrix of the error term -- instead of
#' one figure for the whole model. A model with many regressors would otherwise
#' produce panels too small to read.
#'
#' @return A plot per block of coefficients.
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
plot.bvecmodel <- function(x, ci = 0.95, type = "hist", show_zero_y = TRUE,
                           max_cols = 6, ...) {

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

  periods <- ifelse(tvp, tt, 1)

  # Draws of a coefficient of the measurement equation, either as a single
  # vector or, for time varying parameters, as one column per period
  coeff_draws <- function(pos) {
    x[["posterior"]][["a"]][["coeffs"]][, ncoeffs * 0:(periods - 1) + pos, drop = !tvp]
  }

  blocks <- list()
  regressors <- .get_regressor_blocks_bvecmodel(x)

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

    blocks[["Pi"]] <- list(title = regressors[["Pi"]][["title"]],
                           labels = regressors[["Pi"]][["labels"]],
                           panel = function(i) {
                             .plot_coefficient_panel(Pi[, n_pi * 0:(periods - 1) + i, drop = !tvp],
                                                     tvp, type, ci_low, ci_high, show_zero_y)
                           })
  }

  # Coefficients of the regressors outside the error correction term ----

  offset <- 0
  for (i in setdiff(names(regressors), c("Pi", "A0"))) {

    spec <- regressors[[i]]

    blocks[[i]] <- list(title = spec[["title"]],
                        labels = spec[["labels"]],
                        # The loadings come first among the coefficients, and
                        # 'offset' counts the regressors of the preceding blocks.
                        panel = local({
                          pos_0 <- n_alpha + offset * k
                          function(i) {
                            .plot_coefficient_panel(coeff_draws(pos_0 + i), tvp, type,
                                                    ci_low, ci_high, show_zero_y)
                          }
                        }))

    offset <- offset + length(spec[["labels"]])
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
    pos_a <- n_alpha + k * n_x + temp[lower.tri(temp)]

    blocks[["A0"]] <- list(title = regressors[["A0"]][["title"]],
                           labels = regressors[["A0"]][["labels"]],
                           panel = function(i) {
                             if (i %in% pos_values) {
                               .plot_coefficient_panel(coeff_draws(pos_a[sum(pos_values <= i)]),
                                                       tvp, type, ci_low, ci_high, show_zero_y)
                             } else {
                               graphics::plot.new()
                               graphics::text(0.5, 0.5, labels = ifelse(i %in% pos_one, 1, 0),
                                              adj = 0.5)
                             }
                           })
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

  blocks[["Sigma"]] <- list(title = "Covariance matrix of the error term",
                            labels = y_names,
                            panel = function(i) {
                              if (sv | tvp_and_covar) {
                                stats::plot.ts(u_sigma[kk * 0:(tt - 1) + i, ], plot.type = "single")
                              } else {
                                if (all(u_sigma[, i] == u_sigma[1, i])) {
                                  graphics::plot.new()
                                  graphics::text(0.5, 0.5, labels = u_sigma[1, i], adj = 0.5)
                                } else {
                                  .plot_coefficient_panel(u_sigma[, i], FALSE, type,
                                                          ci_low, ci_high, show_zero_y)
                                }
                              }
                            })

  .plot_blocks(blocks, row_names = y_names, title = title_text, max_cols = max_cols)
}
