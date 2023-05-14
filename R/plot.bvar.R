#' Plotting Draws of a Bayesian VAR Model
#' 
#' A plot function for objects of class \code{"bvar"}.
#' 
#' @param x an object of class \code{"bvar"}, usually, a result of a call to \code{\link{draw_posterior}}.
#' @param ci interval used to calculate credible bands for time-varying parameters.
# @param style the 'layout' of the plot. If \code{style = 1} (default), all parameter draws are displayed in one large plot.
# If \code{style = 2}, multiple panels are generated.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param ... further graphical parameters.
#' 
#' @examples
#' 
#' # Load data 
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Generate model
#' model <- gen_var(e1, p = 1, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' # Plot draws
#' plot(object)
#' 
#' @export
#' @rdname bvar
plot.bvar <- function(x, ci = 0.95, type = "hist", ...) {
  
  if (!type %in% c("hist", "trace", "boxplot")) {
    stop("Argument 'type' must be 'hist', 'trace' or 'boxplot'.")
  }
  
  k <- x[["specifications"]][["dims"]][["K"]]
  m <- x[["specifications"]][["dims"]][["M"]]
  p <- x[["specifications"]][["lags"]][["p"]]
  s <- x[["specifications"]][["lags"]][["s"]]
  n <- 0
  if (!is.null(x[["C"]])) {
    if (x[["specifications"]][["tvp"]][["C"]]) {
      n <- NCOL(x[["C"]][[1]]) / k
    } else {
      n <- NCOL(x[["C"]]) / k
    } 
  }
  tt <- nrow(x[["y"]])
  tsp_info <- stats::tsp(x[["y"]])
  structural <- !is.null(x[["A0"]])
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  y_names <- dimnames(x[["y"]])[[2]]
  x_names <- .get_regressor_names_var(x, add_block = TRUE)
  lab_size <- .05
  mar_orig <- graphics::par("mar")
  
  n_tot <- 0
  if ("A" %in% names(x)) {
    n_tot <- n_tot + k * p
  }
  if ("B" %in% names(x)) {
    n_tot <- n_tot + m * (s + 1)
  }
  if ("C" %in% names(x)) {
    n_tot <- n_tot + n
  }
  if ("A0" %in% names(x)) {
    n_tot <- n_tot + k
  }
  if ("Sigma" %in% names(x)) {
    n_tot <- n_tot + k
  }
  
  mat <- matrix(NA_integer_, k + 2 , n_tot + 1)
  mat[1, ] <- 1
  mat[-1, 1] <- c(0, 2:(k + 1))
  mat[2, -1] <- (k + 1) + 1:n_tot
  mat[-(1:2), -1] <- matrix(1:(k * n_tot) + k + n_tot + 1, k, n_tot)
  graphics::layout(mat,
                   widths = c(lab_size, rep((1 - lab_size) / n_tot, n_tot)),
                   heights = c(.07, lab_size, rep((1 - lab_size) / k, k)))
  
  # Title
  title_text <- "Bayesian " 
  tvp <- any(unlist(x[["specifications"]][["tvp"]])[c("A0", "A", "B", "C")])
  if (tvp) {
    title_text <- paste0(title_text, "TVP-")
  }
  if (x[["specifications"]][["tvp"]][["Sigma"]]) {
    title_text <- paste0(title_text, "SV-")
  }
  if (x[["specifications"]][["structural"]]) {
    title_text <- paste0(title_text, "S")
  }
  title_text <- paste0(title_text, "VAR model")
  p_text <- paste0("p = ", x[["specifications"]][["lags"]][["p"]])
  s_text <- NULL
  if (x[["specifications"]][["dims"]][["M"]] > 0) {
    s_text <- paste0("s = ", x[["specifications"]][["lags"]][["s"]])
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
  if ("Sigma" %in% names(x)) {
    for (j in y_names) {
      graphics::plot.new(); graphics::text(0.5, 0.5, labels = paste0("Sigma\n", j), adj = 0.5)
    } 
  }
  
  graphics::par(mar = c(3, 2.1, .5, 1))
  
  if ("A" %in% names(x)) {
    for (i in 1:p) {
      var_pos <- ((i - 1) * k * k) + 1:(k * k)
      if (x[["specifications"]][["tvp"]][["A"]]) {
        for (j in var_pos) {
          draws <- .tvpribbon(x[["A"]], j, ci_low, ci_high)
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in var_pos) {
          if (all(x[["A"]][, j] == x[["A"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["A"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["A"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["A"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["A"]][, j])
            }
          } 
        }
      }
    }
  }
  
  if ("B" %in% names(x)) {
    for (i in 1:(s + 1)) {
      var_pos <- ((i - 1) * k * m) + 1:(k * m)
      if (x[["specifications"]][["tvp"]][["B"]]) {
        for (j in var_pos) {
          draws <- .tvpribbon(x[["B"]], j, ci_low, ci_high)
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in var_pos) {
          if (all(x[["B"]][, j] == x[["B"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["B"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["B"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["B"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["B"]][, j])
            }
          }
        }
      }
    }
  }
  
  if ("C" %in% names(x)) {
    if (x[["specifications"]][["tvp"]][["C"]]) {
      for (j in 1:NCOL(x[["C"]][[1]])) {
        draws <- .tvpribbon(x[["C"]], j, ci_low, ci_high)
        stats::tsp(draws) <- tsp_info
        stats::ts.plot(draws, xlab = "")
      }
    } else {
      for (j in 1:NCOL(x[["C"]])) {
        if (all(x[["C"]][, j] == x[["C"]][1, j])) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["C"]][1, j], adj = 0.5)
        } else {
          if (type == "hist") {
            graphics::hist(x[["C"]][, j], plot = TRUE, main = NA)  
          }
          if (type == "trace") {
            stats::ts.plot(x[["C"]][, j], xlab = "")
          }
          if (type == "boxplot") {
            graphics::boxplot(x[["C"]][, j])
          }
        }
      }
    }
  }
  
  if ("A0" %in% names(x)) {
    if (x[["specifications"]][["tvp"]][["A0"]]) {
      for (j in 1:NCOL(x[["A0"]][[1]])) {
        draws <- .tvpribbon(x[["A0"]], j, ci_low, ci_high)
        if (all(draws[, 1] == draws[1, 1])) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = draws[1, 2], adj = 0.5)
        } else {
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "")
        }
      }
    } else {
      for (j in 1:NCOL(x[["A0"]])) {
        if (all(x[["A0"]][, j] == x[["A0"]][1, j])) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["A0"]][1, j], adj = 0.5)
        } else {
          if (type == "hist") {
            graphics::hist(x[["A0"]][, j], plot = TRUE, main = NA)   
          }
          if (type == "trace") {
            stats::ts.plot(x[["A0"]][, j], xlab = "")
          }
          if (type == "boxplot") {
            graphics::boxplot(x[["A0"]][, j])
          }
        }
      }
    }
  }
  
  if ("Sigma" %in% names(x)) {
    var_pos <- 1:(k * k)
    if (x[["specifications"]][["tvp"]][["Sigma"]]) {
      for (j in var_pos) {
        draws <- .tvpribbon(x[["Sigma"]], j, ci_low, ci_high)
        if (all(draws[, 1] == draws[1, 1])) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = draws[1, 2], adj = 0.5)
        } else {
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "") 
        }
      }
    } else {
      for (j in var_pos) {
        if (all(x[["Sigma"]][, j] == x[["Sigma"]][1, j])) {
          graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["Sigma"]][1, j], adj = 0.5)
        } else {
          if (type == "hist") {
            graphics::hist(x[["Sigma"]][, j], plot = TRUE, main = NA)  
          }
          if (type == "trace") {
            stats::ts.plot(x[["Sigma"]][, j], xlab = "")
          }
          if (type == "boxplot") {
            graphics::boxplot(x[["Sigma"]][, j])
          }
        }
      }
    }
  }
  
  graphics::par(mar = mar_orig)
  graphics::layout(matrix(1))
}


