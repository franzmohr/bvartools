#' Add Priors to Bayesian Models
#'
#' Adds prior specifications to a BVAR model, which was produced by
#' function \code{\link{create_bvarmodel}}.
#'
#' @param object a list of class 'bvarmodel'.
#' @param coef a named list of prior specifications for the coefficients of the
#' models. For the default specification, all prior means are set to zero and the diagonal elements of
#' the inverse prior variance-covariance matrix are set to 1 for coefficients corresponding to non-deterministic
#' and structural terms. For deterministic coefficients the prior variances are set to 10 via \code{v_i_det = 0.1}.
#' The variances need to be specified as precisions, i.e. as inverses of the variances.
#' For further specification choices for TVP or the Minnesota prior see 'Details'.
#' @param sigma a named list of prior specifications for the error variance-covariance matrix.
#' For the default specification of an inverse Wishart distribution
#' the prior degrees of freedom are set to the number of endogenous variables and
#' the prior variances to 1. See 'Details'.
#' @param varsel optional; a named list of prior specifications for the variable
#' selection algorithm. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The arguments of the function require named lists. Possible
#' specifications are described in the following. Note that it is important to specify the
#' priors in the correct list. Otherwise, the provided specification will be disregarded
#' and default values will be used.
#' 
#' Argument \code{coef} can contain the following elements:
#' \describe{
#'   \item{\code{v_i}}{a numeric specifying the prior precision of the coefficients. Default is 1.
#'   Will be ignored if \code{minnesota} is specified.}
#'   \item{\code{v_i_det}}{a numeric specifying the prior precision of coefficients
#'   corresponding to deterministic terms. Default is 0.1. Will be ignored if \code{minnesota} is specified.}
#'   \item{\code{coint_var}}{a logical specifying whether the prior mean of the first own lag of an
#'   endogenous variable should be set to 1, which is commonly used for cointegrated
#'   VAR models. Default is \code{FALSE}.}
#'   \item{\code{const}}{a numeric or character specifying the prior mean of coefficients, which correspond
#'   to the intercept. If a numeric is provided, all prior means are set to this value.
#'   If \code{coef$const = "mean"}, the mean of the respective endogenous variable is used as prior mean.
#'   If \code{coef$const = "first"}, the first values of the respective endogenous variable is used as prior mean.
#'   This can be useful for trend estimation when using a TVP model.}
#'   \item{\code{minnesota}}{a list containing for parameters for the calculation of
#'   the Minnesota prior, where the element names must be \code{kappa1}, \code{kappa2}, \code{kappa3} and \code{kappa4}.
#'   For the endogenous variable \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#'   \deqn{ \frac{\kappa_{1}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#'   \deqn{ \frac{\kappa_{1} \kappa_{2}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{1} \kappa_{3}}{(l+1)^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for exogenous variables,}}
#'   \deqn{ \kappa_{1} \kappa_{4} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#'   where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#'   LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#'   If the model does not contain exogenous variables, \code{kappa2} will be ignored.}
#'   \item{\code{max_var}}{a numeric specifying the maximum prior variance that is allowed for
#'   non-deterministic coefficients. This argument can be used to overwrite prior variances that are
#'   based on LS estimates, but should not exceed a certain threshold. This can
#'   be useful when using the BVS algorithm.}
#'   \item{\code{shape}}{a numeric specifying the prior shape parameter of the error term of the
#'   state equation. Only used for models with time varying parameters. Default is 3.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error term of the
#'   state equation. Only used for models with time varying parameters. Default is 0.0001.}
#'   \item{\code{rate_det}}{a numeric specifying the prior rate parameter of the error term of the
#'   state equation for coefficients, which correspond to deterministic terms.
#'   Only used for models with time varying parameters. Default is 0.01.}
#' }
#' If \code{minnesota} is specified, \code{v_i} and \code{v_i_det} are ignored.
#' 
#' Argument \code{sigma} can contain the following elements:
#' \describe{
#'   \item{\code{df}}{an integer or character specifying the prior degrees of freedom of the error term. Only
#'   used, if the prior is inverse Wishart.
#'   Default is \code{"k"}, which indicates the amount of endogenous variables in the respective model.
#'   \code{"k + 3"} can be used to set the prior to the amount of endogenous variables plus 3. If an integer
#'   is provided, the degrees of freedom are set to this value in all models.}
#'   \item{\code{scale}}{a numeric specifying the prior error variance of endogenous variables.
#'   Default is 1.}
#'   \item{\code{shape}}{a numeric or character specifying the prior shape parameter of the error term. Only
#'   used, if the prior is inverse gamma or if time varying volatilities are estimated.
#'   For models with constant volatility the default is \code{"k"}, which indicates the amount of endogenous
#'   variables of the model. \code{"k + 3"} can be used to set the prior to the amount of
#'   endogenous variables plus 3. If a numeric is provided, the shape parameters are set to this value in all
#'   models. For models with stochastic volatility this prior refers to the error variance of the state
#'   equation.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error term. Only used, if the
#'   prior is inverse gamma or if time varying volatilities are estimated. For models with stochastic
#'   volatility this prior refers to the error variance of the state equation.}
#'   \item{\code{mu}}{numeric of the prior mean of the initial state of the log-volatilities.
#'   Only used for models with time varying volatility.}
#'   \item{\code{v_i}}{numeric of the prior precision of the initial state of the log-volatilities.
#'   Only used for models with time varying volatility.}
#'   \item{\code{state_variance}}{numeric of the initial draw for the variance of the log-volatilities.
#'   Only used for models with time varying volatility.}
#'   \item{\code{offset}}{numeric of the constant, which is added before taking the log of the squared errors.
#'   Only used for models with time varying volatility.}
#' }
#' For structural models only a gamma prior or stochastic volatility specification is allowed.
#' 
#' Argument \code{varsel} can contain the following elements:
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model.}
#'   \item{\code{covar}}{logical indicating if the variable selection algorithm
#'   should also be applied to the error covariance matrix.}
#'   \item{\code{exclude_det}}{logical indicating if deterministic terms should
#'   be excluded from the variable selection algorithm.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters
#'   for the calculation of the Minnesota-like inclusion priors. See below.}
#'   \item{\code{tau}}{a numeric vector of two elements containing the prior standard errors
#'   of restricted variables (\eqn{\tau_0}) as its first element and unrestricted variables (\eqn{\tau_1})
#'   as its second. Only used for SSVS.}
#'   \item{\code{semiautomatic}}{an numeric vector of two elements containing the
#'   factors by which the standard errors associated with an unconstrained least squares
#'   estimate of the model are multiplied to obtain the prior standard errors
#'   of restricted (\eqn{\tau_0}) and unrestricted (\eqn{\tau_1}) variables, respectively.
#'   This is the semiautomatic approach described in George et al. (2008). Only used for SSVS.}
#' }
#' In the case of SSVS, either \code{tau} or \code{semiautomatic} must be specified.
#' 
#' If \code{varsel$minnesota} is specified, prior
#' inclusion probabilities are calculated in a Minnesota-like fashion as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_{1}}{l}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_{2}}{l}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_{3}}{1 + l}} \tab for exogenous variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables, 
#' }
#' for lag \eqn{l} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3},
#' \eqn{\kappa_4} as the first, second, third and forth element in
#' \code{varsel$minnesota}, respectively.
#' 
#' @return An object of class 'bvarmodel'.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @examples 
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 10, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' @export
add_priors.bvarmodel <- function(object,
                                 coef,
                                 sigma,
                                 varsel = NULL,
                                 ...){
  
  # Input checks
  ## Coefficient priors ----
  if (!is.null(coef)) {
    
    .add_priors_check_coef(object, coef)
    
    if (!is.null(coef[["const"]])) {
      if ("character" %in% class(coef[["const"]])) {
        if (!coef[["const"]] %in% c("first", "mean")) {
          stop("Invalid specificatin of coef$const.")
        }
      }
    }
    
    if (is.null(coef[["v_i_det"]])) {
      coef[["v_i_det"]] <- coef[["v_i"]]
    }
    
  }
  
  ## Sigma priors ----
  error_prior <- .add_priors_check_sigma(object, sigma)
  
  # Check Minnesota ----
  minnesota <- FALSE # Minnesota prior?
  if (!is.null(coef[["minnesota"]])) {
    minnesota <- TRUE
  }
  
  # Check coint VAR ----
  coint_var <- FALSE # Cointegrated VAR?
  if (!is.null(coef[["coint_var"]])) {
    if (coef[["coint_var"]]) {
      coint_var <- TRUE 
    }
  }
  
  if (!is.null(varsel) & object[["model"]][["varsel"]] == "none") {
    stop("Argument 'varsel' was provided, but according to argument 'object$model$varsel' no variable selection algorithm should be used.")
  }
  
  # Check SSVS ----
  use_ssvs <- FALSE
  use_ssvs_error <- FALSE
  use_ssvs_semi <- FALSE
  if (object[["model"]][["varsel"]] == "ssvs") {
    
    .add_priors_check_ssvs(object, varsel)
    
    if (!is.null(varsel[["covar"]])) {
      if (varsel[["covar"]]) {
        use_ssvs_error <- TRUE 
      }
    }
    
    if (is.null(varsel[["exclude_det"]])) {
      varsel[["exclude_det"]] <- FALSE
    }
    # In case ssvs is specified, check if the semi-automatic approach of 
    # George et al. (2008) should be used
    if (!is.null(varsel[["semiautomatic"]])) {
      use_ssvs_semi <- TRUE
    }
    
    use_ssvs <- TRUE
    if (minnesota) {
      minnesota <- FALSE
      warning("Minnesota prior specification overwritten by SSVS.")
    }
  }
  
  # BVS ----
  use_bvs <- FALSE
  use_bvs_error <- FALSE
  if (object[["model"]][["varsel"]] == "bvs") {
    use_bvs <- TRUE
    
    .add_priors_check_bvs(object, varsel)
    
    if (is.null(varsel[["exclude_det"]])) {
      varsel[["exclude_det"]] <- FALSE
    }
    if (!is.null(varsel[["covar"]])) {
      if (varsel[["covar"]]) {
        use_bvs_error <- TRUE 
      }
    }
    if (coef[["v_i"]] == 0 | (coef[["v_i_det"]] == 0 & !varsel[["exclude_det"]])) {
      warning("Using BVS with an uninformative prior is not recommended.")
    }
  }
  
  if (use_ssvs & use_bvs) {
    stop("SSVS and BVS cannot be applied at the same time.")
  }
  
  if (error_prior == "wishart" & (use_ssvs_error | use_bvs_error)) {
    stop("Wishart prior not allowed when BVS or SSVS are applied to covariance matrix.")
  }
  
  varsel_covar <- use_ssvs_error | use_bvs_error
  
  # Generate priors ----
  
  # Get model specs to obtain total number of coeffs
  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  
  if (k == 1 & (use_ssvs_error | use_bvs_error)) {
    stop("BVS or SSVS cannot be applied to covariance matrix when there is only one endogenous variable.")
  } 
  
  use_exo <- FALSE
  if (object[["model"]][["m"]] > 0) {
    use_exo <- TRUE
    m <- object[["model"]][["m"]]
    s <- object[["model"]][["s"]]
  } else {
    s <- 0
    m <- 0
  }
  if (use_exo) {
    s <- s + 1
  }
  
  # Total # of non-deterministic coefficients
  n_a <- k * (k * p)
  n_b <- k * (m * s)
  
  # Add number of non-cointegration deterministic terms
  n_det <- 0
  if (object[["model"]][["n"]] > 0){
    n_det <- object[["model"]][["n"]] * k
  }
  
  tot_par <- n_a + n_b + n_det
  
  covar <- object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar")
  structural <- object[["model"]][["structural"]]
  if (covar & structural) {
    stop("Error covariances and structural coefficients cannot be estimated at the same time.")
  }
  n_struct <- 0
  if (structural & k > 1) {
    n_struct <- (k - 1) * k / 2
    tot_par <- tot_par + n_struct
  }
  
  # Additional input check
  if (!is.null(object[["data"]][["train"]][["z"]])) {
    if (tot_par != ncol(object[["data"]][["train"]][["z"]])) {
      stop("Model specifications are not consistent with data matrix 'object$data$train$z'.")
    } 
  }
  
  sv <- error_prior == "sv"
  
  # Priors ----
  ## Coefficients ----
  if (tot_par > 0) {
    
    ### Prior means ----
    
    mu <- matrix(rep(0, tot_par - n_struct), k)
    
    # Add 1 to first own lags for cointegrated VAR model
    if (coint_var & p > 0) {
      mu[1:k, 1:k] <- diag(1, k)
    }
    
    # Prior for intercept terms
    if (n_det > 0) {
      
      if (!is.null(coef[["const"]]))  {
        
        pos <- which(dimnames(object[["data"]][["train"]][["x"]])[[2]] == "const")
        
        if (length(pos) == 1) {
          if ("character" %in% class(coef[["const"]])) {
            if (coef[["const"]] == "first") {
              mu[, pos] <- object[["data"]][["train"]][["y"]][1, ]
            }
            if (coef[["const"]] == "mean") {
              mu[, pos] <- colMeans(object[["data"]][["train"]][["y"]])
            }
          }
          if ("numeric" %in% class(coef[["const"]])) {
            if (length(coef[["const"]]) == 1 | length(coef[["const"]]) == k) {
              mu[, pos] <- coef[["const"]]  
            } else {
              stop("When a numeric is provided in argument 'coef$const', it must be either a single number or a vector of the same length as the number of endogenous varibles in the model.")
            }
          } 
        }
      }
    }
    
    mu <- matrix(mu)
    
    if (structural) {
      mu <- rbind(mu, matrix(0, n_struct))
    }
    
    object[["priors"]][["a"]] <- list(type = "normal",
                                      mu = mu)
    
    ### Prior covariances ----
    
    if (minnesota) {
      #### Minnesota prior ----
      minn <- minnesota_prior(object = object,
                              kappa1 = coef[["minnesota"]][["kappa1"]],
                              kappa2 = coef[["minnesota"]][["kappa2"]],
                              kappa3 = coef[["minnesota"]][["kappa3"]],
                              kappa4 = coef[["minnesota"]][["kappa4"]],
                              max_var = NULL,
                              coint_var = FALSE,
                              sigma = "AR")
      
      object[["priors"]][["a"]][["v_inv"]] <- minn[["v_inv"]]
    }
    
    #### SSVS prior ----
    if (use_ssvs) {
      
      if (object[["model"]][["tvp"]]) {
        stop("SSVS is not supported for TVP models.")
      }
      
      if (sv) {
        stop("Not allowed to use SSVS with stochastic volatility models.")
      }
      
      ssvs_temp <- ssvs_prior(object, tau = varsel[["tau"]], semiautomatic = varsel[["semiautomatic"]])
      temp <- inclusion_prior(object, prob = varsel[["inprior"]], exclude_deterministics = varsel[["exclude_det"]],
                              minnesota_like = !is.null(varsel[["minnesota"]]),
                              kappa1 = varsel[["minnesota"]][1],
                              kappa2 = varsel[["minnesota"]][2],
                              kappa3 = varsel[["minnesota"]][3],
                              kappa4 = varsel[["minnesota"]][4])
      object[["priors"]][["a"]][["v_inv"]] <- diag(1 / ssvs_temp[["tau1"]][, 1]^2, tot_par)
      object[["priors"]][["a"]][["inprior"]] <- temp[["prior"]]
      object[["priors"]][["a"]][["include"]] <- temp[["include"]]
      object[["priors"]][["a"]][["tau0"]] <- ssvs_temp[["tau0"]]
      object[["priors"]][["a"]][["tau1"]] <- ssvs_temp[["tau1"]]
      rm(temp)
      rm(ssvs_temp)
    }
    
    #### Regular prior ----
    if (!minnesota & !use_ssvs) {
      v_i <- diag(coef[["v_i"]], tot_par)
      # Add priors for deterministic terms if they were specified
      if (n_det > 0 & !is.null(coef[["v_i_det"]])) {
        diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef[["v_i_det"]]
      }
      object[["priors"]][["a"]][["v_inv"]] <- v_i
    }
    
    #### BVS prior ----
    if (use_bvs) {
      temp <- inclusion_prior(object, prob = varsel[["inprior"]], exclude_deterministics = varsel[["exclude_det"]],
                              minnesota_like = !is.null(varsel[["minnesota"]]),
                              kappa1 = varsel[["minnesota"]][1],
                              kappa2 = varsel[["minnesota"]][2],
                              kappa3 = varsel[["minnesota"]][3],
                              kappa4 = varsel[["minnesota"]][4])
      object[["priors"]][["a"]][["inprior"]] <- temp[["prior"]]
      object[["priors"]][["a"]][["include"]] <- temp[["include"]]
      rm(temp)
    }
    
    ### TVP prior ----
    if (object[["model"]][["tvp"]]) {
      object[["priors"]][["a"]][["shape"]] <- matrix(coef[["shape"]], tot_par)
      object[["priors"]][["a"]][["rate"]] <- matrix(coef[["rate"]], tot_par)
      if (n_det > 0 & !is.null(coef[["rate_det"]])) {
        object[["priors"]][["a"]][["rate"]][tot_par - n_struct - n_det + 1:n_det, ] <- coef[["rate_det"]]
      }
    }
  }
  
  ## Covar priors ----
  
  if (!structural & covar & k > 1) {
    
    n_covar <- k * (k - 1) / 2
    object[["priors"]][["psi"]][["type"]] <- "normal"
    object[["priors"]][["psi"]][["mu"]] <- matrix(0, n_covar)
    object[["priors"]][["psi"]][["v_inv"]] <- diag(coef[["v_i"]], n_covar)
    if (object[["model"]][["tvp"]]) {
      object[["priors"]][["psi"]][["shape"]] <- matrix(coef[["shape"]], n_covar)
      object[["priors"]][["psi"]][["rate"]] <- matrix(coef[["rate"]], n_covar) 
    }
    
    # Variable selection
    object[["priors"]][["psi"]][["varsel"]] <- "none"
    
    # SSVS priors
    if (use_ssvs_error) {
      object[["priors"]][["psi"]][["varsel"]] <- "ssvs"
      object[["priors"]][["psi"]][["inprior"]] <- matrix(varsel[["inprior"]], n_covar)
      object[["priors"]][["psi"]][["include"]] <- matrix(1:n_covar)
      object[["priors"]][["psi"]][["tau0"]] <- matrix(varsel[["tau"]][1], n_covar)
      object[["priors"]][["psi"]][["tau1"]] <- matrix(varsel[["tau"]][2], n_covar)
    }
    
    # BVS priors
    if (use_bvs_error) {
      object[["priors"]][["psi"]][["varsel"]] <- "bvs"
      object[["priors"]][["psi"]][["inprior"]] <- matrix(varsel[["inprior"]], n_covar)
      object[["priors"]][["psi"]][["include"]] <- matrix(1:n_covar)
    }
  }
  
  ## Error term ----
  if (error_prior == "sv") {
    
    object <- .add_priors_sv_helper(object, sigma, k)
    
  } else {
    
    if (error_prior == "wishart") {
      object[["priors"]][["u_sigma"]][["type"]] <- "wishart"
      help_df <- sigma[["df"]]
      object[["priors"]][["u_sigma"]][["df"]] <- NA_real_
      object[["priors"]][["u_sigma"]][["scale"]] = diag(sigma[["scale"]], k)
    }
    
    if (error_prior == "gamma") {
      object[["priors"]][["u_sigma"]][["type"]] <- "gamma"
      help_df <- sigma[["shape"]]
      object[["priors"]][["u_sigma"]][["shape"]] <- NA_real_
      object[["priors"]][["u_sigma"]][["rate"]] = matrix(sigma[["rate"]], k)
    }
    
    if (minnesota & !is.null(object[["data"]][["train"]][["x"]])) {
      # Store LS estimate of variance covariance matrix for analytical solution
      object[["priors"]][["u_sigma"]][["u_sigma_inv"]] = minn[["sigma_inv"]]
    }
    
    if ("character" %in% class(help_df)) {
      if (grepl("k", help_df)) {
        # Transform character specification to expression and evaluate
        help_df <- eval(parse(text = help_df))
      } else {
        stop("Use no other letter than 'k' in 'sigma$df' to indicate the number of endogenous variables.")
      }
    }
    
    if (help_df < 0) {
      stop("Current specification implies a negative prior degree of\nfreedom or shape parameter of the error term.")
    }
    
    if (error_prior == "wishart") {
      object[["priors"]][["u_sigma"]][["df"]] <- as.integer(help_df)
    }
    if (error_prior == "gamma") {
      object[["priors"]][["u_sigma"]][["shape"]] <- matrix(help_df, k)
    } 
  }
  
  return(object)
}