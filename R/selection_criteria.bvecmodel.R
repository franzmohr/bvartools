#' Model Selection Criteria
#'
#' Calculates model selection criteria for an object of class 'bvecmodel'.
#'
#' @param object an object of class 'bvecmodel'.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'selcrit'.
#' 
#' @examples
#' 
#' # Load data
#' data("e6")
#' 
#' # Create model
#' model <- create_bvecmodel(e6, p = 4,
#'                           const = "unrestricted",
#'                           seasonal = "unrestricted",
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
#' # Add log-likelihoods
#' model <- add_posterior_loglik(model)
#' 
#' # Calculate selection criteria
#' sel <- selection_criteria(model)
#' sel
#' 
#' 
#' 
#' @export
#' @method selection_criteria bvecmodel
selection_criteria.bvecmodel <- function(object, ci = 0.95, ...){
  
  if (ci < 0 | ci > 1) {
    stop("Argument 'ci' is not within the permitted range of 0 and 1.")
  }
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  
  # Model information
  k <- object[["model"]][["k"]]
  rank <- object[["model"]][["rank"]]
  structural <- object[["model"]][["structrual"]]
  tt <- nrow(object[["data"]][["train"]][["y"]])
  varnames <- dimnames(object[["data"]][["original"]][["endogen"]])[[2]]
  if (!is.null(object[["data"]][["train"]][["x"]])) {
    nparams <- rank + ncol(object[["data"]][["train"]][["x"]])
  } else {
    nparams <- rank
  }
  h <- object[["model"]][["h"]]
  max_n_columns <- k * h
  
  use_ll <- !is.null(object[["posterior"]][["loglik"]])
  
  errors <- get_forecast_errors(object)
  use_fe <- !is.null(errors)
  
  if (!use_ll & !use_fe) {
    stop("Model object must contain at least either posterior draws of the log-likelihood or forecast errors.")
  }
  
  result <- NULL
  result[["model"]] <- object[["model"]]
  
  if (use_ll) {
    
    loglik <- rowSums(object[["posterior"]][["loglik"]])
    
    # Log-likelihood
    ll <- data.frame("mean" = mean(loglik),
                     "median" = stats::median(loglik),
                     "qlower" = stats::quantile(loglik, probs = ci_low),
                     "qupper" = stats::quantile(loglik, probs = ci_high))
    row.names(ll) <- NULL
    result[["LL"]] <- ll
    
    # AIC
    aic <- 2 * nparams - 2 * loglik
    aic <- data.frame("mean" = mean(aic),
                      "median" = stats::median(aic),
                      "qlower" = stats::quantile(aic, probs = ci_low),
                      "qupper" = stats::quantile(aic, probs = ci_high))
    row.names(aic) <- NULL
    result[["AIC"]] <- aic
    
    
    # BIC
    bic <- nparams * log(tt) - 2 * loglik
    bic <- data.frame("mean" = mean(bic),
                      "median" = stats::median(bic),
                      "qlower" = stats::quantile(bic, probs = ci_low),
                      "qupper" = stats::quantile(bic, probs = ci_high))
    row.names(bic) <- NULL
    result[["BIC"]] <- bic
    
    
    # HQ
    hq <- 2 * nparams * log(log(tt)) - 2 * loglik
    hq <- data.frame("mean" = mean(hq),
                     "median" = stats::median(hq),
                     "qlower" = stats::quantile(hq, probs = ci_low),
                     "qupper" = stats::quantile(hq, probs = ci_high))
    row.names(hq) <- NULL
    result[["HQ"]] <- hq
    
  }
  
  
  
  if (use_fe) {
    
    h <- ncol(errors) / k
    
    # Forecast errors
    result[["FE"]] <- as.data.frame(matrix(NA, ncol(errors), 6))
    names(result[["FE"]]) <- c("variable", "h", "mean", "median", "qlower", "qupper")
    result[["FE"]][, "variable"] <- rep(varnames, h)
    result[["FE"]][, "h"] <- rep(1:h, each = k)
    result[["FE"]][, "mean"] <- apply(errors, 2, mean, na.rm = TRUE)
    result[["FE"]][, "median"] <- apply(errors, 2, stats::median, na.rm = TRUE)
    result[["FE"]][, "qlower"] <- apply(errors, 2, stats::quantile, probs = ci_low, na.rm = TRUE)
    result[["FE"]][, "qupper"] <- apply(errors, 2, stats::quantile, probs = ci_high, na.rm = TRUE)
    
    # Absolute errors
    errors <- abs(errors)
    result[["AFE"]] <- as.data.frame(matrix(NA, ncol(errors), 6))
    names(result[["AFE"]]) <- c("variable", "h", "mean", "median", "qlower", "qupper")
    result[["AFE"]][, "variable"] <- rep(varnames, h)
    result[["AFE"]][, "h"] <- rep(1:h, each = k)
    result[["AFE"]][, "mean"] <- apply(errors, 2, mean, na.rm = TRUE)
    result[["AFE"]][, "median"] <- apply(errors, 2, stats::median, na.rm = TRUE)
    result[["AFE"]][, "qlower"] <- apply(errors, 2, stats::quantile, probs = ci_low, na.rm = TRUE)
    result[["AFE"]][, "qupper"] <- apply(errors, 2, stats::quantile, probs = ci_high, na.rm = TRUE)
    
    # Squared errors
    errors <- errors^2
    result[["RSFE"]] <- as.data.frame(matrix(NA, ncol(errors), 6))
    names(result[["RSFE"]]) <- c("variable", "h", "mean", "median", "qlower", "qupper")
    result[["RSFE"]][, "variable"] <- rep(varnames, h)
    result[["RSFE"]][, "h"] <- rep(1:h, each = k)
    result[["RSFE"]][, "mean"] <- sqrt(apply(errors, 2, mean, na.rm = TRUE))
    result[["RSFE"]][, "median"] <- sqrt(apply(errors, 2, stats::median, na.rm = TRUE))
    result[["RSFE"]][, "qlower"] <- sqrt(apply(errors, 2, stats::quantile, probs = ci_low, na.rm = TRUE))
    result[["RSFE"]][, "qupper"] <- sqrt(apply(errors, 2, stats::quantile, probs = ci_high, na.rm = TRUE)) 
  }
  
  attr(result, "ci") <- c(paste0(ci_low * 100, "%"), paste0(ci_high * 100, "%"))
  class(result) <- c("selcrit", "list")
  
  return(result)
}