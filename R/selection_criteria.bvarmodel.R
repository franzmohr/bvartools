#' Model Selection Criteria
#'
#' Calculates model selection criteria for an object of class 'bvarmodel'.
#'
#' @param object an object of class 'bvarmodel'.
#' @param ci a numeric between 0 and 1 specifying the probability of the credible band.
#' Defaults to 0.95.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'selcrit', which also inherits the class of the model, with the
#' element \code{model} and one data frame per criterion. If the model contains
#' \code{posterior$loglik}, these are \code{LL}, \code{AIC}, \code{BIC}, \code{HQ},
#' \code{WAIC} and \code{LOOIC}, each with the columns \code{mean}, \code{median},
#' \code{qlower} and \code{qupper}, where bands that do not apply are \code{NA}. If it
#' contains \code{posterior$forecast_errors}, these are \code{FE}, \code{AFE} and
#' \code{RSFE}, the forecast errors and their absolute and root squared values, with the
#' columns \code{variable}, \code{h}, \code{mean}, \code{median}, \code{qlower} and
#' \code{qupper}.
#'
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 1:2, deterministic = "const",
#'                           iterations = 10, burnin = 10)
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
#' @method selection_criteria bvarmodel
selection_criteria.bvarmodel <- function(object, ci = 0.95, ...){
  
  if (ci < 0 | ci > 1) {
    stop("Argument 'ci' is not within the permitted range of 0 and 1.")
  }
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  
  
  # Model information
  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n"]]
  structural <- object[["model"]][["structrual"]]
  tt <- nrow(object[["data"]][["train"]][["y"]])
  varnames <- dimnames(object[["data"]][["original"]][["endogen"]])[[2]]
  # Free parameters of the model. 'k * p + m * (s + 1) + n' is the number of
  # regressors of one equation, and there are k equations, so the coefficients
  # come to k times that. The error term contributes k(k + 1)/2, whether it is
  # an unrestricted covariance matrix or, in a structural model, a diagonal one
  # of k elements beside the k(k - 1)/2 free elements of A0, which come to the
  # same number.
  nparams <- k * (k * p + m * (s + 1) + n) + k * (k + 1) / 2
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
    deviance <- .plugin_deviance(object[["posterior"]][["loglik"]])

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
