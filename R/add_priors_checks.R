
.add_priors_check_varsel_names <- function(varsel) {

  allowed_varsel_arguments <- c("inprior", "covar", "exclude_det", "minnesota",
                                "tau", "semiautomatic")
  for (i in names(varsel)) {
    if (!i %in% allowed_varsel_arguments) {
      stop(paste0("Element '", i, "' in argument 'varsel' is not recognised."))
    }
  }
}

.add_priors_check_bvs <- function(object, bvs) {

  use_bvs_error <- FALSE

  if (object[["model"]][["varsel"]] == "bvs") {

    if (is.null(bvs)) {
      stop("BVS was chosen as variable selection algorithm, but no prior specification was provided.")
    }

    .add_priors_check_varsel_names(bvs)

    if (is.null(bvs[["inprior"]])) {
      stop("Argument 'bvs$inprior' must be specified for BVS")
    }
    
    if (!is.null(bvs[["covar"]])) {
      if (bvs[["covar"]]) {
        use_bvs_error <- TRUE 
      }
    }
  }
}

.add_priors_check_coef <- function(object, coef) {
  
  allowed_coef_arguments <- c("v_i", "v_i_det", "const", "minnesota",
                              "max_var", "shape", "rate", "rate_det")
  # Only the Minnesota prior of a VAR model can centre the first own lag on one.
  if (inherits(object, "bvarmodel")) {
    allowed_coef_arguments <- c(allowed_coef_arguments, "coint_var")
  }
  for (i in names(coef)) {
    if (!i %in% allowed_coef_arguments) {
      stop(paste0("Element '", i, "' in argument 'coef' is not recognised."))
    }
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
    if (!"minnesota" %in% names(coef)) {
      stop("If 'coef$v_i' is not specified, 'coef$minnesota' must be specified.")
    }
  }
  
  # Tests for specifications used in TVP models
  if (!is.null(object[["model"]][["tvp"]])) {
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
  }
  
  if (!is.null(coef[["minnesota"]])) {
    if (!"list" %in% class(coef[["minnesota"]])) {
      stop("Argument coef$minnesota must be a named list.")
    }
    if (is.null(names(coef[["minnesota"]]))) {
      stop("Argument coef$minnesota must be a named list.")
    }
    for (i in names(coef[["minnesota"]])) {
      if (!i %in% c("kappa1", "kappa2", "kappa3", "kappa4")) {
        stop(paste0("Element '", i, "' in argument 'coef$minnesota' is not recognised."))
      }
    }
    if (!all(c("kappa1", "kappa2", "kappa4") %in% names(coef[["minnesota"]]))) {
      stop("Argument coeff$minnesota must contain at least the elements 'kappa1', 'kappa2' and 'kappa4'.")
    }
    if (object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar") & is.null(coef[["v_i"]])) {
      stop("If error covarances should be estimated, argument coef$v_i must be provided also when the Minnesota prior is used.")
    }
  }
}


.add_priors_check_sigma <- function(object, sigma) {
  
  if (is.null(sigma)) {
    stop("Argument 'sigma' may not be NULL.")
  }

  allowed_sigma_arguments <- c("df", "scale", "shape", "rate", "mu", "v_i",
                               "state_variance", "offset")
  for (i in names(sigma)) {
    if (!i %in% allowed_sigma_arguments) {
      stop(paste0("Element '", i, "' in argument 'sigma' is not recognised."))
    }
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
      .add_priors_check_per_equation(sigma, c("shape", "rate"), object$model$k)
      if (is.numeric(sigma$shape) && any(sigma$shape < 0)) {
        stop("Argument 'sigma$shape' must be at least 0.")
      }
      if (any(sigma$rate <= 0)) {
        stop("Argument 'sigma$rate' must be larger than 0.")
      }
    }
    
    if (object$model$error %in% c("sv", "sv+covar")) {
      if (any(!c("mu", "v_i", "shape", "rate", "state_variance", "offset") %in% names(sigma))) {
        stop("Missing prior specifications for stochastic volatility prior.")
      }
      error_prior <- "sv"
    }
    
    if (object$model$error == "ald") {
      if (all(c("shape", "rate") %in% names(sigma))) {
        error_prior <- "ald"
      } else {
        stop("An asymmetric Laplace model requires specification of elements 'shape' and 'rate' in 'sigma'.")
      }
      .add_priors_check_per_equation(sigma, c("shape", "rate"), object$model$k)
      if (any(sigma$shape <= 0)) {
        stop("Argument 'sigma$shape' must be larger than 0.")
      }
      if (any(sigma$rate <= 0)) {
        stop("Argument 'sigma$rate' must be larger than 0.")
      }
    }
    
    if (object$model$error == "wishart") {
      if (all(c("df", "scale") %in% names(sigma))) {
        error_prior <- "wishart"
      } else {
        stop("Wishart prior requires specification of elements 'df' and 'scale' in 'sigma'.")
      }
      if (object[["model"]][["structural"]]) {
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




# The gamma priors of the error term are stored one per equation, so each of
# their elements is either one value for every equation or one per equation.
# A vector used to reach a scalar comparison and stop with "the condition has
# length > 1".
.add_priors_check_per_equation <- function(sigma, elements, k) {
  for (i in elements) {
    if (is.numeric(sigma[[i]]) && !length(sigma[[i]]) %in% c(1, k)) {
      stop("Argument 'sigma$", i, "' must contain either one value or one per endogenous ",
           "variable (", k, "), but it contains ", length(sigma[[i]]), ".", call. = FALSE)
    }
  }
}


.add_priors_check_ssvs <- function(object, ssvs) {
  
  use_ssvs_error <- FALSE

  .add_priors_check_varsel_names(ssvs)

  if (is.null(ssvs[["inprior"]])) {
    stop("Argument 'varsel$inprior' must be specified for SSVS.")
  }
  if (is.null(ssvs[["tau"]]) & is.null(ssvs[["semiautomatic"]])) {
    stop("Either argument 'varsel$tau' or 'varsel$semiautomatic' must be specified for SSVS.")
  }
  if (!is.null(ssvs[["covar"]])) {
    if (ssvs[["covar"]]) {
      use_ssvs_error <- TRUE 
    }
  }
  if (object[["model"]][["error"]] == "gamma+covar" & use_ssvs_error & is.null(ssvs[["tau"]])) {
    stop("If SSVS should be applied to error covariances, argument 'varsel$tau' must be specified.")
  }
  if (!is.null(ssvs[["semiautomatic"]])) {
    if (!"numeric" %in% class(ssvs[["semiautomatic"]])) {
      stop("Argument 'varsel$semiautomatic' must be a numeric vector.")
    }
    if (length(ssvs[["semiautomatic"]]) != 2) {
      stop("Argument 'varsel$semiautomatic' must be a numeric vector with two elements.")
    }
  }
}