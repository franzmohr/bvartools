#' Plotting Draws of a Bayesian VAR Model
#'
#' A plot function for objects of class 'bvarmodel'.
#'
#' @param x an object of class 'bvarmodel'.
#' @param ci interval used to calculate credible bands for time-varying parameters.
# @param style the 'layout' of the plot. If \code{style = 1} (default), all parameter draws are displayed in one large plot.
# If \code{style = 2}, multiple panels are generated.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param show_zero_y if \code{TRUE} (default), a horizontal line with y = 0 is
#' added to the plot. Only used for time varying parameters.
#' @param max_cols an integer of the maximum number of regressors per figure. A block
#' with more regressors than this is drawn as several figures of nearly equal width.
#' Defaults to 6.
#' @param ... further graphical parameters.
#'
#' @details The function draws one figure per block of coefficients -- the lags of
#' the endogenous variables, the exogenous variables, the deterministic terms, the
#' contemporaneous endogenous variables of a structural model, and the covariance
#' matrix of the error term -- instead of one figure for the whole model. A model
#' with many regressors would otherwise produce panels too small to read.
#'
#' @return A plot per block of coefficients.
#'
#' @examples
#'
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#'
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
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
#'
#' @export
plot.bvarmodel <- function(x, ci = 0.95, type = "hist", show_zero_y = TRUE,
                           max_cols = 6, ...) {

  # 'layout' is called below, so all parameters have to be restored on exit
  orig_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(orig_par))

  if (!type %in% c("hist", "trace", "boxplot")) {
    stop("Argument 'type' must be 'hist', 'trace' or 'boxplot'.")
  }

  k <- x[["model"]][["k"]]
  kk <- k * k
  p <- x[["model"]][["p"]]
  m <- x[["model"]][["m"]]
  s <- x[["model"]][["s"]]
  n <- x[["model"]][["n"]]
  tvp <- x[["model"]][["tvp"]]
  tvp_and_covar <- tvp & x[["model"]][["error"]] == "gamma+covar"
  sv <- x[["model"]][["error"]] %in% c("sv", "sv+covar")
  structural <- x[["model"]][["structural"]]
  if (structural) {
    n_struct <- k * (k - 1) / 2
  } else {
    n_struct <- 0
  }

  tt <- nrow(x[["data"]][["train"]][["y"]])

  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  y_names <- dimnames(x[["data"]][["original"]][["endogen"]])[[2]]

  n_nonstruct <- k * (k * p + m * (s + 1) + n)
  ncoeffs <- n_nonstruct + n_struct

  # Title
  title_text <- "Bayesian "
  if (tvp) {
    title_text <- paste0(title_text, "TVP-")
  }
  if (sv) {
    title_text <- paste0(title_text, "SV-")
  }
  if (x[["model"]][["error"]] == "ald") {
    title_text <- paste0(title_text, "Quantile-")
  }
  if (structural) {
    title_text <- paste0(title_text, "S")
  }
  title_text <- paste0(title_text, "VAR model")
  p_text <- paste0("p = ", p)
  s_text <- NULL
  if (m > 0) {
    s_text <- paste0("s = ", s)
  }
  # The quantile is part of what the model is, not of how it was fitted, so it
  # is reported beside the lag orders rather than left to the specification.
  q_text <- NULL
  if (!is.null(x[["model"]][["quantile"]])) {
    q_text <- paste0("q = ", x[["model"]][["quantile"]])
  }
  if (any(!is.null(c(p_text, s_text, q_text)))) {
    lag_text <- paste0(c(p_text, s_text, q_text), collapse = " and ")
  } else {
    lag_text <- NULL
  }
  title_text <- paste0(c(title_text, lag_text), collapse = " with ")

  periods <- ifelse(tvp, tt, 1)

  # Draws of one coefficient, either as a vector or, for a time varying model,
  # as one column per period
  coeff_draws <- function(pos) {
    x[["posterior"]][["a"]][["coeffs"]][, ncoeffs * 0:(periods - 1) + pos, drop = !tvp]
  }

  blocks <- list()

  # Coefficients of the regressors ----
  regressors <- .get_regressor_blocks_bvarmodel(x)
  offset <- 0
  for (i in names(regressors)) {

    spec <- regressors[[i]]

    if (i == "A0") {
      next
    }

    blocks[[i]] <- list(title = spec[["title"]],
                        labels = spec[["labels"]],
                        # 'offset' is the number of regressors of the preceding
                        # blocks, and a block holds k coefficients per regressor.
                        panel = local({
                          pos_0 <- offset * k
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

    # Position of the free elements of A0 among the coefficients
    temp <- matrix(NA, k , k)
    temp[upper.tri(temp)] <- 1:n_struct
    temp <- t(temp)
    pos_a <- n_nonstruct + temp[lower.tri(temp)]

    blocks[["A0"]] <- list(title = "Contemporaneous endogenous variables",
                           labels = y_names,
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
