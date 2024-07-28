
.add_priors_check_coef <- function(object, coef) {
  
  if (!all(names(coef) %in% c("v_i", "v_i_det", "coint_var", "const", "minnesota", "max_var", "shape", "rate", "rate_det"))) {
    stop("At least one element in argument 'coef' is not valid.")
  }
  
  if (!is.null(coef[["v_i"]])) {
    if (coef[["v_i"]] < 0) {
      stop("Argument 'v_i' must be at least 0.")
    } 
    # Define "v_i_det" if not specified (needed for a check later)
    if (is.null(coef[["v_i_det"]])) {
      coef[["v_i_det"]] <- coef[["v_i"]]
    }
  } else {
    if (!any(c("minnesota", "ssvs") %in% names(coef))) {
      stop("If 'coef$v_i' is not specified, at least 'coef$minnesota' or 'coef$ssvs' must be specified.")
    }
  }
  
  # Tests for specifications used in TVP models
  if (object[["model"]][["tvp"]]) {
    if (!"shape" %in% names(coef)) {
      stop("Argument 'coef$shape' must be specified for TVP models.")
    }
    if (!"rate" %in% names(coef)) {
      stop("Argument 'coef$rate' must be specified for TVP models.")
    }
    
    if (!is.null(coef[["const"]])) {
      if ("character" %in% class(coef[["const"]])) {
        if (!coef[["const"]] %in% c("first", "mean")) {
          stop("Invalid specificatin of coef$const.")
        }
      }
    }
  }
  
  if (!is.null(coef[["minnesota"]])) {
    if (!"list" %in% class(coef[["minnesota"]])) {
      stop("Argument coef$minnesota must be a named list.")
    }
    if (is.null(names(coef[["minnesota"]]))) {
      stop("Argument coef$minnesota must be a named list.")
    }
    if (!all(c("kappa0", "kappa1", "kappa3") %in% names(coef[["minnesota"]]))) {
      stop("Argument coeff$minnesota must contain at least the elements 'kappa0', 'kappa1' and 'kappa3'.")
    }
    if (object$model$error %in% c("gamma+covar", "sv+covar") & is.null(coef[["v_i"]])) {
      stop("If error covarances should be estimated, argument coef$v_i must be provided also when the Minnesota prior is used.")
    }
  }
}


.add_priors_check_sigma <- function(object, sigma) {
  
  if (is.null(sigma)) {
    stop("Argument 'sigma' may not be NULL.")
  }
  
  if (length(sigma) < 2) {
    stop("Argument 'sigma' must be at least of length 2.")
  } else {
    error_prior <- NULL
    
    if (object$model$error %in% c("gamma", "gamma+covar")) {
      if (all(c("shape", "rate") %in% names(sigma))) {
        error_prior <- "gamma"
      } else {
        stop("Gamma prior requires specification of elements 'shape' and 'rate' in 'sigma'.")
      }
      if (sigma$shape < 0) {
        stop("Argument 'sigma$shape' must be at least 0.")
      }
      if (sigma$rate <= 0) {
        stop("Argument 'sigma$rate' must be larger than 0.")
      } 
    }
    
    if (object$model$error %in% c("sv", "sv+covar")) {
      if (any(!c("mu", "v_i", "shape", "rate", "state_variance", "offset") %in% names(sigma))) {
        stop("Missing prior specifications for stochastic volatility prior.")
      }
      error_prior <- "sv"
    }
    
    if (object$model$error == "wishart") {
      if (all(c("df", "scale") %in% names(sigma))) {
        error_prior <- "wishart"
      } else {
        stop("Wishart prior requires specification of elements 'df' and 'scale' in 'sigma'.")
      }
      if (any(unlist(lapply(object, function(x) {x$model$structural})))) {
        stop("Structural models may not use a Wishart prior. Consider using a gamma prior instead.")
      }
      if (sigma$df < 0) {
        stop("Argument 'sigma$df' must be at least 0.")
      }
      if (sigma$scale <= 0) {
        stop("Argument 'sigma$scale' must be larger than 0.")
      } 
    }
    
    if (is.null(error_prior)) {
      stop("Invalid specification for argument 'sigma'.")
    }
  }
  
  return(error_prior)
}

.add_priors_check_ssvs <- function(object, ssvs) {
  
  use_ssvs_error <- FALSE
  
  if (is.null(ssvs[["inprior"]])) {
    stop("Argument 'ssvs$inprior' must be specified for SSVS.")
  }
  if (is.null(ssvs[["tau"]]) & is.null(ssvs[["semiautomatic"]])) {
    stop("Either argument 'ssvs$tau' or 'ssvs$semiautomatic' must be specified for SSVS.")
  }
  if (!is.null(ssvs[["covar"]])) {
    if (ssvs[["covar"]]) {
      use_ssvs_error <- TRUE 
    }
  }
  if (object[["model"]][["error"]] == "gamma+covar" & use_ssvs_error & is.null(ssvs[["tau"]])) {
    stop("If SSVS should be applied to error covariances, argument 'ssvs$tau' must be specified.")
  }
  if (object[["model"]][["structural"]] & is.null(ssvs[["tau"]])) {
    stop("If SSVS should be used with structural models, argument 'ssvs$tau' must be specified.")
  }
  
  if (!is.null(ssvs[["semiautomatic"]])) {
    if (!"numeric" %in% class(ssvs[["semiautomatic"]])) {
      stop("Argument 'ssvs$semiautomatic' must be a numeric vector.")
    }
    if (length(ssvs[["semiautomatic"]]) != 2) {
      stop("Argument 'ssvs$semiautomatic' must be a numeric vector with two elements.")
    }
  }
}