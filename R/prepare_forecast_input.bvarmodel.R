#' Prepare Forecast Input
#' 
#' Generates data matrices serving as input for forecasting simulation for objects
#' of class 'bvarmodel'.
#' 
#' @param object an object of class 'bvarmodel'.
#' @param n_ahead number of steps ahead at which to predict.
#' @param deterministic a time-series object with deterministic data. If not
#' specified, the function will try to identify the deterministic terms
#' automatically. If this is not successful, an error message we be returned.
#' @param exogen a time-series object with unmodelled, non-deterministic data.
#' See 'Details'.
#' @param ... additional arguments.
#' 
#' @return A list.
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
#' # Generate forcast input
#' fcst_input <- prepare_forecast_input(model, n_ahead = 4)
#' 
#' 
#' @export
prepare_forecast_input.bvarmodel <- function(object, n_ahead = 10, deterministic = NULL, exogen = NULL, ...) {
  
  # Input checks
  if (n_ahead < 1) {
    stop("Argument 'n_ahead' must be at least 1.")
  }
  
  if (!is.null(exogen)) {
    if (!"ts" %in% class(exogen)) {
      stop("Argument 'exogen' must be of class 'ts'.")
    }
  }
  
  if (!is.null(deterministic)) {
    if (!"ts" %in% class(deterministic)) {
      stop("Argument 'deterministic' must be of class 'ts'.")
    }
  }
  
  # Model specs
  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n"]]
  y <- object[["data"]][["train"]][["y"]]
  z <- object[["data"]][["train"]][["z"]]
  tt <- nrow(y)
  y_tsp <- stats::tsp(y)
  y_end <- y_tsp[2]
  y_freq <- y_tsp[3]
  
  if (m > 0 & is.null(exogen)) {
    stop("If parameters of unmodeled, non-deterministic variables are estimated, argument 'exogen' must be specified.")
  }
  
  if (n > 0 & is.null(deterministic)) {
    
    # Deducing deterministic terms from model input
    if (p > 0 | m > 0) {
      det_names <- dimnames(object[["data"]][["train"]][["x"]])[[2]][-c(1:(k * p + m * (s + 1)))]
    } else {
      det_names <- dimnames(object[["data"]][["train"]][["x"]])[[2]] 
    }
    
    deterministic <- matrix(NA, n_ahead + 1, n)
    dimnames(deterministic) <- list(NULL, det_names)
    
    ## constant
    if ("const" %in% det_names) {
      if (all(object[["data"]][["train"]][["x"]][, "const"] == 1)) {
        deterministic[, "const"] <- rep(1, n_ahead + 1)
      }
    }
    
    ## trend
    if ("trend" %in% det_names) {
      increment <- object[["data"]][["train"]][["x"]][tt, "trend"] - object[["data"]][["train"]][["x"]][tt - 1, "trend"]
      deterministic[, "trend"] <- seq(from = object[["data"]][["train"]][["x"]][tt, "trend"], by = increment, length.out = n_ahead + 1)
    }
    
    ## seasonal
    if (any(grepl("season", det_names, fixed = TRUE))) {
      # Detect seasonal pattern
      pos_season <- det_names[which(grepl("season", det_names, fixed = TRUE))]
      last_obs <- object[["data"]][["train"]][["x"]][tt, pos_season]
      for (i in 1:tt) {
        if (identical(last_obs, object[["data"]][["train"]][["x"]][i, pos_season])) {
          deterministic[, pos_season] <- object[["data"]][["train"]][["x"]][i + 0:n_ahead, pos_season]
          break
        }
      }
    }
    
    if (any(is.na(deterministic))) {
      stop("Could not identify all deterministic terms. Please specify argument 'deterministic' instead.") 
    }
    
    deterministic <- stats::ts(deterministic, start = y_end, frequency = y_freq)
  }
  
  n_tot <- k * p + m * (s + 1) + n
  
  pred_start <- stats::time(stats::ts(rep(NA, 2), start = y_end, frequency = y_freq))[-1]
  
  z <- NULL
  if (n_tot > 0) {
    z <- matrix(NA, k * n_ahead, n_tot)
    
    x <- stats::ts(matrix(NA, n_ahead, n_tot), start = pred_start, frequency = y_freq)
    x_time <- stats::time(x)
    
    if (p > 0) {
      temp_p <- stats::embed(object[["data"]][["original"]][["endogen"]], p)
      temp_p <- temp_p[nrow(temp_p),]
      for (i in 1:p) {
        if (i <= n_ahead) {
          x[i, ((i - 1) * k + 1):(p * k)] <- temp_p[1:((p - i + 1) * k)] 
        }
      }
    }
    
    if (m > 0) {
      tsp_x <- stats::tsp(exogen)
      temp_x <- stats::ts(stats::embed(exogen, s + 1), end = tsp_x[2], frequency = tsp_x[3])
      temp_x <- stats::window(temp_x, start = x_time[1], end = x_time[length(x_time)])
      x[, k * p + 1:(m * (s + 1))] <- temp_x
    }
    
    if (n > 0) {
      deterministic <- stats::window(deterministic, start = x_time[1], end = x_time[length(x_time)])
      x[, k * p + m * (s + 1) + 1:n] <- deterministic
    }
    
    z <- kronecker(x, diag(1, k))
  }
  
  result <- list("h" = as.integer(n_ahead),
                 "z" = z)
  
  return(result)
}
