#' Model Selection Criteria
#'
#' Calculates model selection criteria for an object of class 'bvecmodel'.
#'
#' @param object an object of class 'bvecmodel'.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'selcrit', which also inherits the class of the model, with the
#' element \code{model} and one data frame per criterion, as described in
#' \code{\link{selection_criteria.bvarmodel}}.
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
#' @family model comparison
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
  tt <- nrow(object[["data"]][["train"]][["y"]])
  varnames <- dimnames(object[["data"]][["original"]][["endogen"]])[[2]]
  # Free parameters of the model.
  #
  # Pi = alpha beta' is a k x k_ect matrix of rank 'rank'. A matrix of that
  # shape and rank has rank * (k + k_ect - rank) free elements rather than
  # k * k_ect of them, because alpha and beta are identified only up to an
  # r x r rotation, and it is this term -- and only this term -- that grows
  # with the rank. The remaining coefficients are one per equation and column
  # of x. The error term contributes k(k + 1)/2, whether it is an unrestricted
  # covariance matrix or, in a structural model, a diagonal one of k elements
  # beside the k(k - 1)/2 free elements of A0, which come to the same number.
  n_x <- 0L
  if (!is.null(object[["data"]][["train"]][["x"]])) {
    n_x <- ncol(object[["data"]][["train"]][["x"]])
  }
  k_ect <- 0L
  if (!is.null(object[["data"]][["train"]][["w"]])) {
    k_ect <- ncol(object[["data"]][["train"]][["w"]])
  }
  if (rank > 0 && k_ect == 0) {
    stop("The free parameters of the cointegration term cannot be counted ",
         "without the error correction term in 'data$train$w'.")
  }
  nparams <- rank * (k + k_ect - rank) + k * n_x + k * (k + 1) / 2
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

    # Information criteria
    #
    # The criteria are evaluated at the point estimate of the model rather than
    # averaged over the posterior, which would charge the complexity of the
    # model a second time. See .plugin_deviance().
    deviance <- .plugin_deviance(object)

    # AIC
    result[["AIC"]] <- .point_criterion(deviance + 2 * nparams)

    # BIC
    result[["BIC"]] <- .point_criterion(deviance + log(tt) * nparams)

    # HQ
    result[["HQ"]] <- .point_criterion(deviance + 2 * log(log(tt)) * nparams)

    # WAIC
    #
    # Penalises by the flexibility the fit actually used rather than by a count
    # of parameters, which is what makes constant, time varying and stochastic
    # volatility specifications comparable with each other. See .waic().
    waic <- .waic_data_frame(object[["posterior"]][["loglik"]], ci_low, ci_high)
    if (!is.null(waic)) {
      row.names(waic) <- NULL
      result[["WAIC"]] <- waic
    }

    # LOOIC
    #
    # Estimates the same out-of-sample deviance as WAIC, by reweighting the
    # posterior towards the one that has not seen a period rather than by a
    # correction term, and reports how far that reweighting can be trusted.
    # See .psis_loo().
    looic <- .loo_data_frame(object[["posterior"]][["loglik"]], ci_low, ci_high)
    if (!is.null(looic)) {
      row.names(looic) <- NULL
      result[["LOOIC"]] <- looic
    }

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
  # The classes of the model are maintained, so that methods, which use the model
  # specifications, can be dispatched on them
  class(result) <- c("selcrit", class(object))
  
  return(result)
}
