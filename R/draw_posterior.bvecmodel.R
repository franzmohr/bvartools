#' Posterior Simulation for Vector Error Correction Models
#' 
#' Forwards model input to posterior simulation functions for vector error correction models.
#' 
#' @param object a list of model specifications, which should be passed on
#' to function \code{FUN}. Usually, the output of a call to \code{\link{gen_vec}}
#' in combination with \code{\link{add_priors}}.
#' @param posterior_function the function to be applied to each list element in argument \code{object}.
#' If \code{NULL} (default), the internal function \code{\link{bvecpost}} is used.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return n object of the class of the output of the applied posterior
#' simulation function. In case the package's own function is used, this will
#' result in an object of class \code{"bvec"}.
#' 
#' @references
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference in
#' a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
#' \doi{10.1016/j.jeconom.2011.07.007}
#' 
#' @examples
#' 
#' # Load data 
#' data("e6")
#' e6 <- e6 * 100
#' 
#' # Generate model
#' model <- create_vec_model(e6, p = 1, r = 1, const = "restricted",
#'                           iterations = 10, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(model)
#' 
#' @export
draw_posterior.bvecmodel <- function(object, posterior_function = NULL, ...){
  
  model <- object
  if ("posteriors" %in% names(model)) {
    model[["posteriors"]] <- NULL
  }
  
  if (is.null(posterior_function)) {
    object <- try(bvecpost(object))
  } else {
    # Apply own function
    object <- try(posterior_function(object))
  }
  
  # Produce something if estimation fails
  if (inherits(object, "try-error")) {
    object <- c(model, list(error = TRUE))
  } else {
    
    if (is.null(posterior_function)) {
      
      k <- object[["model"]][["k"]]
      p <- object[["model"]][["p"]]
      m <- object[["model"]][["m"]]
      s <- object[["model"]][["s"]]
      n_gamma <- k * k * (p - 1)
      n_upsilon <- k * m * s
      n_unres <- k * object[["model"]][["n_unrestricted"]]
      n_a0 <- ifelse(object[["model"]][["structural"]], k * (k - 1) / 2, 0)
      n_res <- object[["model"]][["n_restricted"]]
      tt <- NROW(object[["data"]][["y"]])
      r <- object[["model"]][["rank"]]
      n_alpha <- r * k
      n_a <- n_alpha + n_gamma + n_upsilon + n_unres + n_a0
      pos <- which(lower.tri(diag(1, k)))
      draws <- object[["model"]][["iterations"]]
      tvp <- object[["model"]][["tvp"]]
      
      # Input data preparation ----
      
      tsp_temp <- stats::tsp(object[["data"]][["y"]])
      
      w <- stats::ts(as.matrix(object[["data"]][["w"]][, 1:k]), class = c("mts", "ts", "matrix"))
      stats::tsp(w) <- tsp_temp
      dimnames(w)[[2]] <- dimnames(object[["data"]][["w"]])[[2]][1:k]
      
      if (m > 0) {
        w_x <- stats::ts(as.matrix(object[["data"]][["w"]][, k + 1:m]), class = c("mts", "ts", "matrix"))
        stats::tsp(w_x) <- tsp_temp
        dimnames(w_x)[[2]] <- dimnames(object[["data"]][["w"]])[[2]][k + 1:m]
      } else {
        w_x <- NULL
      }
      
      if (n_res > 0) {
        w_d <- stats::ts(as.matrix(object[["data"]][["w"]][, k + m + 1:n_res]), class = c("mts", "ts", "matrix"))
        stats::tsp(w_d) <- tsp_temp
        dimnames(w_d)[[2]] <- dimnames(object[["data"]][["w"]])[[2]][k + m + n_res]
      } else {
        w_d <- NULL
      }
      
      x <- NULL
      x_x <- NULL
      x_d <- NULL
      
      if (!is.null(object[["data"]][["x"]])) {
        
        if (p > 1) {
          x <- stats::ts(as.matrix(object[["data"]][["x"]][, 1:(k * (p - 1))]), class = c("mts", "ts", "matrix"))
          stats::tsp(x) <- tsp_temp
          dimnames(x)[[2]] <- dimnames(object[["data"]][["x"]])[[2]][1:(k * (p - 1))]
        }
        
        
        if (m > 0) {
          x_x <- stats::ts(as.matrix(object[["data"]][["x"]][, k * (p - 1) + 1:(m * s)]), class = c("mts", "ts", "matrix"))
          stats::tsp(x_x) <- tsp_temp
          dimnames(x_x)[[2]] <- dimnames(object[["data"]][["x"]])[[2]][(k * (p - 1)) + 1:(m * s)]
        }
        
        if (object[["model"]][["n_unrestricted"]] > 0) {
          n_ur <- object[["model"]][["n_unrestricted"]]
          x_d <- stats::ts(as.matrix(object[["data"]][["x"]][, k * (p - 1) + m * s + 1:n_ur]), class = c("mts", "ts", "matrix"))
          stats::tsp(x_d) <- tsp_temp
          dimnames(x_d)[[2]] <- dimnames(object[["data"]][["x"]])[[2]][k * (p - 1) + m * s + 1:n_ur]
        }
      }
      
      # Coefficient draws ----
      
      ## alpha ----
      alpha <- NULL
      beta <- NULL
      beta_x <- NULL
      beta_d <- NULL
      n_beta <- r * k
      n_beta_x <- r * m
      n_beta_c <- r * n_res
      k_ect <- k + m + n_res
      n_beta_tot <- n_beta + n_beta_x + n_beta_c
      if (r > 0) {
        
        pos_alpha <- 1:n_alpha
        if (tvp) {
          pos_alpha <- rep(pos_alpha, tt) + rep(0:(tt - 1), each = length(pos_alpha)) * n_a 
        }
        alpha <- t(object[["posteriors"]][["a"]][["coeffs"]][, pos_alpha])
        
        pos_beta <- rep((1:r - 1) * k_ect, each = k) + rep(1:k, r)
        if (tvp) {
          pos_beta <- rep(pos_beta, tt) + rep(0:(tt - 1), each = length(pos_beta)) * n_beta_tot
        }
        beta <- t(object[["posteriors"]][["beta"]][["coeffs"]][, pos_beta])
        
        if (m > 0) {
          pos_beta_x <- rep((1:r - 1) * k_ect, each = m) + rep(k + 1:m, r)
          if (tvp) {
            pos_beta_x <- rep(pos_beta_x, tt) + rep(0:(tt - 1), each = length(pos_beta_x)) * n_beta_tot
          }
          beta_x <- t(object[["posteriors"]][["beta"]][["coeffs"]][, pos_beta_x]) 
        }
        
        if (n_res > 0) {
          pos_beta_c <- rep((1:r - 1) * k_ect, each = n_res) + rep(k + m + 1:n_res, r)
          if (tvp) {
            pos_beta_c <- rep(pos_beta_c, tt) + rep(0:(tt - 1), each = length(pos_beta_c)) * n_beta_tot
          }
          beta_d <- t(object[["posteriors"]][["beta"]][["coeffs"]][, pos_beta_c])
        }
      }
      
      # Modelled variables ----
      gamma <- NULL
      if (p > 1) {
        pos_gamma <- n_alpha + 1:n_gamma
        if (tvp) {
          pos_gamma <- rep(pos_gamma, tt) + rep(0:(tt - 1), each = length(pos_gamma)) * n_a
        }
        gamma[["coeffs"]] <- t(object[["posteriors"]][["a"]][["coeffs"]][, pos_gamma])
        
        if ("sigma" %in% names(object[["posteriors"]][["a"]])) {
          gamma[["sigma"]] <- matrix(NA, n_gamma, draws)
          gamma[["sigma"]] <- t(object[["posteriors"]][["a"]][["sigma"]][, pos_gamma[1:n_gamma]])
        }
        
        if ("lambda" %in% names(object[["posteriors"]][["a"]])) {
          gamma[["lambda"]] <- matrix(NA, n_gamma, draws)
          gamma[["lambda"]] <- t(object[["posteriors"]][["a"]][["lambda"]][, pos_gamma[1:n_gamma]])
        }
      }
      
      # Unmodelled variables ----
      upsilon <- NULL
      if (m > 0) {
        pos_upsilon <- n_alpha + n_gamma + 1:n_upsilon
        if (tvp) {
          pos_upsilon <- rep(pos_upsilon, tt) + rep(0:(tt - 1), each = length(pos_upsilon)) * n_a
        }
        upsilon[["coeffs"]] <- t(object[["posteriors"]][["a"]][["coeffs"]][, pos_upsilon])
        
        if ("sigma" %in% names(object[["posteriors"]][["a"]])) {
          upsilon[["sigma"]] <- matrix(NA, n_upsilon, draws)
          upsilon[["sigma"]] <- t(object[["posteriors"]][["a"]][["sigma"]][, pos_upsilon[1:n_upsilon]])
        }
        
        if ("lambda" %in% names(object[["posteriors"]][["a"]])) {
          upsilon[["lambda"]] <- matrix(NA, n_upsilon, draws)
          upsilon[["lambda"]] <- t(object[["posteriors"]][["a"]][["lambda"]][, pos_upsilon[1:n_upsilon]])
        }
      }
      
      # Unrestricted deterministic terms ----
      c_unres <- NULL
      if (n_unres > 0) {
        pos_unres <- n_alpha + n_gamma + n_upsilon + 1:n_unres
        if (tvp) {
          pos_unres <- rep(pos_unres, tt) + rep(0:(tt - 1), each = length(pos_unres)) * n_a
        }
        c_unres[["coeffs"]] <- t(object[["posteriors"]][["a"]][["coeffs"]][, pos_unres])
        
        if ("sigma" %in% names(object[["posteriors"]][["a"]])) {
          c_unres[["sigma"]] <- matrix(NA, n_unres, draws)
          c_unres[["sigma"]] <- t(object[["posteriors"]][["a"]][["sigma"]][, pos_unres[1:n_unres]])
        }
        
        if ("lambda" %in% names(object[["posteriors"]][["a"]])) {
          c_unres[["lambda"]] <- matrix(NA, n_unres, draws)
          c_unres[["lambda"]] <- t(object[["posteriors"]][["a"]][["lambda"]][, pos_unres[1:n_unres]])
        }
      }
      
      # Structural ----
      draws_a0 <- .extract_posterior_a0(object)
      
      # Sigma----
      draws_sigma <- .extract_posterior_sigma(object)
      
      # Create bvar object
      object <- bvec(y = object[["data"]][["y"]],
                     w = w,
                     w_x = w_x,
                     w_d = w_d,
                     alpha = alpha,
                     beta = beta,
                     beta_x = beta_x,
                     beta_d = beta_d,
                     Pi = NULL, Pi_x = NULL, Pi_d = NULL,
                     x = x,
                     x_x = x_x,
                     x_d = x_d,
                     r = r,
                     A0 = draws_a0,
                     Gamma = gamma,
                     Upsilon = upsilon,
                     C = c_unres,
                     Sigma = draws_sigma,
                     data = NULL, exogen = NULL)
    }
  }
  
  return(object)
}