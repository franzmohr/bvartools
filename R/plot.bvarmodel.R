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
#' @param ... further graphical parameters.
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
plot.bvarmodel <- function(x, ci = 0.95, type = "hist", show_zero_y = TRUE, ...) {
  
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
  x_names <- .get_regressor_names_bvarmodel(x, add_block = TRUE)
  lab_size <- .05
  
  nparams <- 0
  if (!is.null(x[["data"]][["train"]][["z"]])) {
    nparams <- ncol(x[["data"]][["train"]][["z"]])
    if (structural) {
      nparams <- nparams - n_struct + k * k
    }
    nparams <- nparams / k
  }
  nparams <- nparams + k
  
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
  sv <- x[["model"]][["error"]] %in% c("sv", "sv+covar") 
  if (sv) {
    title_text <- paste0(title_text, "SV-")
  }
  if (x[["model"]][["structural"]]) {
    title_text <- paste0(title_text, "S")
  }
  title_text <- paste0(title_text, "VAR model")
  p_text <- paste0("p = ", x[["model"]][["p"]])
  s_text <- NULL
  if (x[["model"]][["m"]] > 0) {
    s_text <- paste0("s = ", x[["model"]][["s"]])
  }
  if (any(!is.null(c(p_text, s_text)))) {
    lag_text <- paste0(c(p_text, s_text), collapse = " and ")
  } else {
    lag_text <- NULL
  }
  title_text <- paste0(c(title_text, lag_text), collapse = " with ")
  
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
  
  n_nonstruct <- k * (k * p + m * (s + 1) + n)
  ncoeffs <- n_nonstruct + n_struct
  
  if (n_nonstruct > 0) {
    for (i in 1:n_nonstruct) {
      if (tvp) {
        pos <- ncoeffs * 0:(tt - 1) + i
        temp <- x[["posterior"]][["a"]][["coeffs"]][, pos]
        stats::ts.plot(t(apply(temp, 2, stats::quantile, probs = c(ci_low, .5, ci_high))), xlab = "")
      } else {
        if (type == "hist") {
          graphics::hist(x[["posterior"]][["a"]][["coeffs"]][, i], plot = TRUE, main = NA)  
        }
        if (type == "trace") {
          stats::ts.plot(x[["posterior"]][["a"]][["coeffs"]][, i], xlab = "")
        }
        if (type == "boxplot") {
          graphics::boxplot(x[["posterior"]][["a"]][["coeffs"]][, i])
        } 
      }
    }
  }
  
  if (structural) {
    
    struct_matrix <- matrix(1:(k * k), k)
    pos_values <- which(lower.tri(struct_matrix))
    pos_zero <- which(upper.tri(struct_matrix))
    pos_one <- struct_matrix[-c(pos_values, pos_zero)]
    
    # pos in a
    temp <- matrix(NA, k , k)
    temp[upper.tri(temp)] <- 1:n_struct
    temp <- t(temp)
    pos_a <- n_nonstruct + temp[lower.tri(temp)]
    
    pos_i <- 0
    for (i in 1:(k * k)) {
      if (i %in% pos_values) {
        pos_i <- pos_i + 1
        
        if (tvp) {
          pos <- ncoeffs * 0:(tt - 1) + pos_a[pos_i]
          temp <- x[["posterior"]][["a"]][["coeffs"]][, pos]
          stats::ts.plot(t(apply(temp, 2, stats::quantile, probs = c(ci_low, .5, ci_high))), xlab = "")
        } else {
          if (type == "hist") {
            graphics::hist(x[["posterior"]][["a"]][["coeffs"]][, pos_a[pos_i]], plot = TRUE, main = NA)  
          }
          if (type == "trace") {
            stats::ts.plot(x[["posterior"]][["a"]][["coeffs"]][, pos_a[pos_i]], xlab = "")
          }
          if (type == "boxplot") {
            graphics::boxplot(x[["posterior"]][["a"]][["coeffs"]][, pos_a[pos_i]])
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


