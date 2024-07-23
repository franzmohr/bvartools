#' Add Priors to Bayesian Models
#'
#' Adds prior specifications to a list of models, which was produced by
#' function \code{\link{create_vec_model}}.
#'
#' @param object a list of class 'bvecmodel'.
#' @param coef a named list of prior specifications for coefficients that do not 
#' determine the cointegration space. For the default specification all prior means are set to zero
#' and the diagonal elements of the inverse prior variance-covariance matrix are set to 1
#' for coefficients corresponding to non-deterministic terms. For deterministic coefficients the prior
#' variances are set to 10 by \code{v_i_det = 0.1}. The variances need to be specified as precisions,
#' i.e. as inverses of the variances. For further specifications such as the Minnesota prior see 'Details'.
#' @param coint a named list of prior specifications for coefficients determining the
#' cointegration space of VEC models. See 'Details'.
#' @param sigma a named list of prior specifications for the error variance-covariance matrix.
#' For the default specification of an inverse Wishart distribution the prior degrees
#' of freedom are set to the number of endogenous variables plus the rank of
#' the cointegration matrix. The prior variance is to 1. See 'Details'.
#' @param ssvs optional; a named list of prior specifications for the SSVS algorithm.
#' Not allowed for TVP models. See 'Details'.
#' @param bvs optional; a named list of prior specifications for the BVS algorithm. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The arguments of the function require named lists. Possible
#' specifications are described in the following. Note that it is important to specify the
#' priors in the correct list. Otherwise, the provided specification will be disregarded
#' and default values will be used.
#' 
#' Argument \code{coef} contains the following elements
#' \describe{
#'   \item{\code{v_i}}{a numeric specifying the prior precision of the coefficients. Default is 1.
#'   Will be ignored if \code{minnesota} is specified.}
#'   \item{\code{v_i_det}}{a numeric specifying the prior precision of coefficients
#'   corresponding to deterministic terms. Default is 0.1. Will be ignored if \code{minnesota} is specified.}
#'   \item{\code{const}}{a numeric or character specifying the prior mean of coefficients, which correspond
#'   to the intercept. If a numeric is provided, all prior means are set to this value.
#'   If \code{coef$const = "mean"}, the mean of the respective endogenous variable is used as prior mean.
#'   If \code{coef$const = "first"}, the first values of the respective endogenous variable is used as prior mean.}
#'   \item{\code{minnesota}}{a list of length 4 containing parameters for the calculation of
#'   the Minnesota prior, where the element names must be \code{kappa0}, \code{kappa1}, \code{kappa2} and \code{kappa3}.
#'   For the endogenous variable \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#'   \deqn{ \frac{\kappa_{0}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#'   \deqn{ \frac{\kappa_{0} \kappa_{1}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{0} \kappa_{2}}{(l+1)^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for exogenous variables,}}
#'   \deqn{ \kappa_{0} \kappa_{3} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#'   where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#'   LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#'   The function only provides priors for the non-cointegration part of the model. However,
#'   the residual standard errors \eqn{\sigma_i} are based on an unrestricted LS regression of the
#'   endogenous variables on the error correction term and the non-cointegration regressors.}
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
#' 
#' Argument \code{coint} can contain the following elements:
#' \describe{
#'   \item{\code{v_i}}{numeric between 0 and 1 specifying the shrinkage of the cointegration space prior. Default is 0.}
#'   \item{\code{p_tau_i}}{numeric of the diagonal elements of the inverse prior matrix of
#'   the central location of the cointegration space \eqn{sp(\beta)}. Default is 1.}
#'   \item{\code{rho}}{a numeric specifying the autocorrelation coefficient
#'   of the state equation of \eqn{\beta}. It must be smaller than 1. Default is 0.999.
#'   Note that in contrast to Koop et al. (2011) \eqn{\rho} is not drawn in the Gibbs sampler of
#'   this package yet.}
#' }
#' 
#' Argument \code{sigma} can contain the following elements:
#' \describe{
#'   \item{\code{df}}{an integer or character specifying the prior degrees of freedom of the error term. Only
#'   used, if the prior is inverse Wishart.
#'   Default is \code{"k"}, which indicates the amount of endogenous variables in the respective model.
#'   \code{"k + 3"} can be used to set the prior to the amount of endogenous variables plus 3. If an integer
#'   is provided, the degrees of freedom are set to this value in all models.
#'   In all cases the rank \eqn{r} of the cointegration matrix is automatically added.}
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
#' Argument \code{ssvs} can contain the following elements:
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model.}
#'   \item{\code{tau}}{a numeric vector of two elements containing the prior standard errors
#'   of restricted variables (\eqn{\tau_0}) as its first element and unrestricted variables (\eqn{\tau_1})
#'   as its second.}
#'   \item{\code{semiautomatic}}{an numeric vector of two elements containing the
#'   factors by which the standard errors associated with an unconstrained least squares
#'   estimate of the model are multiplied to obtain the prior standard errors
#'   of restricted (\eqn{\tau_0}) and unrestricted (\eqn{\tau_1}) variables, respectively.
#'   This is the semiautomatic approach described in George et al. (2008).}
#'   \item{\code{covar}}{logical indicating if SSVS should also be applied to the
#'   error covariance matrix.}
#'   \item{\code{exclude_det}}{logical indicating if deterministic terms should
#'   be excluded from the SSVS algorithm.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters
#'   for the calculation of the Minnesota-like inclusion priors. See below.}
#' }
#' Either \code{tau} or \code{semiautomatic} must be specified.
#' 
#' The argument \code{bvs} can contain the following elements
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model.}
#'   \item{\code{covar}}{logical indicating if BVS should also be applied to the
#'   error covariance matrix.}
#'   \item{\code{exclude_det}}{logical indicating if deterministic terms should
#'   be excluded from the BVS algorithm.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters
#'   for the calculation of the Minnesota-like inclusion priors. See below.}
#' }
#' 
#' If either \code{ssvs$minnesota} or \code{bvs$minnesota} is specified, prior
#' inclusion probabilities are calculated in a Minnesota-like fashion as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_1}{l}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_2}{l}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_3}{1 + l}} \tab for exogenous variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables, 
#' }
#' for lag \eqn{l} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3},
#' \eqn{\kappa_4} as the first, second, third and forth element in
#' \code{ssvs$minnesota} or \code{bvs$minnesota}, respectively.
#' 
#' @return A list of models.
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
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference in
#' a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
#' \doi{10.1016/j.jeconom.2011.07.007}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @examples 
#' 
#' # Load data
#' data("e6")
#' 
#' # Create model
#' model <- create_vec_model(e6, p = 4, r = 1,
#'                           const = "unrestricted", seasonal = "unrestricted",
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model)
#' 
#' @export
add_priors.bvecmodel <- function(object,
                                 coef = list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 0.0001, rate_det = 0.01),
                                 coint = list(v_i = 0, p_tau_i = 1, rho = .999),
                                 sigma = list(df = "k", scale = 1, mu = 0, v_i = 0.01, state_variance = 0.05, offset = 0.0001),
                                 ssvs = NULL,
                                 bvs = NULL,
                                 ...){
  
  # Input checks ----
  
  ## coefficients ----
  if (!is.null(coef)) {
    .add_priors_check_coef(object, coef)
  }
  
  ## cointegration ----
  if (object[["model"]][["tvp"]]) {
    if (!"rho" %in% names(coint)) {
      stop("Argument 'coint$rho' must be specified for VEC models with time varying cointegration parameters.")
    }
    if (coint[["rho"]] >= 1) {
      stop("Argument 'coint$rho' must be smaller than 1.")
    }
    if (coint[["rho"]] < .8) {
      warning("Value of argument 'coint$rho' appears rather small.")
    }
  } else {
    if (!"v_i" %in% names(coint)) {
      stop("Argument 'coint$v_i' must be specified for VEC models with constant cointegration parameters.")
    }
    if (!"p_tau_i" %in% names(coint)) {
      stop("Argument 'coint$p_tau_i' must be specified for VEC models with constant cointegration parameters.")
    }
  }
  
  ## sigma ----
  error_prior <- .add_priors_check_sigma(object, sigma)
  
  ## Minnesota ----
  minnesota <- FALSE # Minnesota prior?
  if (!is.null(coef[["minnesota"]])) {
    minnesota <- TRUE
  }
  
  ## SSVS ----
  use_ssvs <- FALSE
  use_ssvs_error <- FALSE
  use_ssvs_semi <- FALSE
  if (!is.null(ssvs)) {
    
    .add_priors_check_ssvs(object, ssvs)
    
    if (!is.null(ssvs[["covar"]])) {
      if (ssvs[["covar"]]) {
        use_ssvs_error <- TRUE 
      }
    }
    
    if (is.null(ssvs[["exclude_det"]])) {
      ssvs[["exclude_det"]] <- FALSE
    }
    # In case ssvs is specified, check if the semi-automatic approach of 
    # George et al. (2008) should be used
    if (!is.null(ssvs[["semiautomatic"]])) {
      use_ssvs_semi <- TRUE
    }
    
    use_ssvs <- TRUE
    if (minnesota) {
      minnesota <- FALSE
      warning("Minnesota prior specification overwritten by SSVS.")
    }
  }
  
  ## BVS ----
  use_bvs <- FALSE
  use_bvs_error <- FALSE
  if (!is.null(bvs)) {
    use_bvs <- TRUE
    if (is.null(bvs$inprior)) {
      stop("If BVS should be applied, argument 'bvs$inprior' must be specified.")
    }
    if (is.null(bvs$exclude_det)) {
      bvs$exclude_det <- FALSE
    }
    if (!is.null(bvs$covar)) {
      if (bvs$covar) {
        use_bvs_error <- TRUE 
      }
    }
    if (coef[["v_i"]] == 0 | (coef[["v_i_det"]] == 0 & !bvs[["exclude_det"]])) {
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
  r <- object[["model"]][["rank"]]
  
  if (k == 1 & (use_ssvs_error | use_bvs_error)) {
    stop("BVS or SSVS cannot be applied to covarianc matrix when there is only one endogenous variable.")
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
  
  # Substract lag from domestic model for VEC
  p <- p - 1
  
  # Total # of non-deterministic coefficients
  n_alpha <- k * r
  n_gamma <- k * (k * p)
  n_upsilon <- k * (m * s)
  
  # Add number of unrestricted deterministic terms
  n_det <- object[["model"]][["n_unrestricted"]] * k
  
  tot_par <- n_alpha + n_gamma + n_upsilon + n_det
  
  covar <- object[["model"]][["error"]] %in% c("gamma-covar", "sv-covar")
  structural <- object[["model"]][["structural"]]
  if (covar & structural) {
    stop("Error covariances and structural coefficients cannot be estimated at the same time.")
  }
  sv <- object[["model"]][["error"]] %in% c("sv", "sv-covar")
  n_struct <- 0
  n_z <- NCOL(object[["data"]][["z"]])
  if (object[["model"]][["rank"]] > 0) {
    n_w <- NCOL(object[["data"]][["w"]])
  }
  
  if (structural & k > 1) {
    n_struct <- (k - 1) * k / 2
    tot_par <- tot_par + n_struct
  }
  
  #### Cointegration (constant) ----
  n_ect <- k * (k + m)
  if (object[["model"]][["n_restricted"]] > 0) {
    n_ect <- n_ect + object[["model"]][["n_restricted"]] * k
  }
  
  r <- object[["model"]][["rank"]]
  if (r > 0) {
    
    n_alpha <- r * k
    n_beta <- r * n_ect / k
    
    if (object[["model"]][["tvp"]]) {
      
      object[["priors"]][["cointegration"]] <- list(type = "cointspace",
                                                    rho = coint[["rho"]],
                                                    mu = matrix(0, n_beta),
                                                    v_i = Matrix::Diagonal(n_beta, 1 - coint[["rho"]]^2))
    } else {
      object[["priors"]][["cointegration"]] <- list(type = "cointspace",
                                                    v_i = coint[["v_i"]],
                                                    p_tau_i = diag(coint[["p_tau_i"]], n_ect / k)) 
    }
  }
  
  # Non-cointegration (constant) ----
  # Generate prior matrices ----
  if (tot_par > 0) {
    
    # Zero prior means
    mu <- matrix(rep(0, tot_par - n_struct), k)
    
    # Prior for intercept terms
    if (n_det > 0) {
      
      if (!is.null(coef[["const"]]))  {
        
        pos <- which(dimnames(object$data$x)[[2]] == "const") + r
        
        if (length(pos) == 1) {
          if ("character" %in% class(coef[["const"]])) {
            if (coef$const == "first") {
              mu[, r + pos] <- object$data$y[1, ]
            }
            if (coef$const == "mean") {
              mu[, r + pos] <- colMeans(object$data$y)
            }
          }
          if ("numeric" %in% class(coef[["const"]])) {
            mu[, r + pos] <- coef$const
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
    
    # Prior covariances
    if (minnesota) {
      # Minnesota prior ----
      minn <- minnesota_prior(object = object,
                              kappa0 = coef[["minnesota"]][["kappa0"]],
                              kappa1 = coef[["minnesota"]][["kappa1"]],
                              kappa2 = coef[["minnesota"]][["kappa2"]],
                              kappa3 = coef[["minnesota"]][["kappa3"]],
                              max_var = NULL,
                              sigma = "AR")
      
      object[["priors"]][["a"]][["v_i"]] <- minn[["v_i"]]
    }
    
    # SSVS prior ----
    if (use_ssvs) {
      
      if (object[["model"]][["tvp"]]) {
        stop("SSVS is not supported for TVP models.")
      }
      
      ssvs_temp <- ssvs_prior(object, tau = ssvs$tau, semiautomatic = ssvs$semiautomatic)
      temp <- inclusion_prior(object, prob = ssvs$inprior, exclude_deterministics = ssvs$exclude_det,
                              minnesota_like = !is.null(ssvs$minnesota), kappa = ssvs$minnesota)
      if (!is.null(temp)) {
        object$model$varselect <- "SSVS"
        object[["priors"]][["a"]][["v_i"]] <- diag(1 / ssvs_temp$tau1[, 1]^2, tot_par)
        object[["priors"]][["a"]][["ssvs"]][["inprior"]] <- temp$prior
        object[["priors"]][["a"]][["ssvs"]][["include"]] <- temp$include
        object[["priors"]][["a"]][["ssvs"]][["tau0"]] <- ssvs_temp$tau0
        object[["priors"]][["a"]][["ssvs"]][["tau1"]] <- ssvs_temp$tau1
      } else {
        object[["priors"]][["a"]][["v_i"]] <- diag(coef[["v_i"]], tot_par)
      }
      rm(temp)
      rm(ssvs_temp)
    }
    
    # Regular prior ----
    if (!minnesota & !use_ssvs) {
      v_i <- diag(coef[["v_i"]], tot_par)
      # Add priors alpha
      if (r > 0 & object[["model"]][["tvp"]]) {
        v_i[1:n_alpha, 1:n_alpha] <- diag(1 / (1 - coint[["rho"]] * coint[["rho"]]), n_alpha)
      }
      # Add priors for deterministic terms if they were specified
      if (n_det > 0 & !is.null(coef[["v_i_det"]])) {
        diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef[["v_i_det"]]
      }
      
      object[["priors"]][["a"]][["v_i"]] <- v_i
    }
    
    if (use_bvs) {
      temp <- inclusion_prior(object, prob = bvs[["inprior"]], exclude_deterministics = bvs[["exclude_det"]],
                              minnesota_like = !is.null(bvs[["minnesota"]]), kappa = bvs[["minnesota"]])
      
      if (!is.null(temp)) {
        object$model$varselect <- "BVS"
        object[["priors"]][["a"]][["bvs"]][["inprior"]] <- temp[["prior"]]
        object[["priors"]][["a"]][["bvs"]][["include"]] <- temp[["include"]] 
      }
    }
    
    
    # TVP prior ----
    if (object[["model"]][["tvp"]]) {
      object[["priors"]][["a"]][["v_i"]] <- Matrix::Matrix(object[["priors"]][["a"]][["v_i"]]) 
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
    object[["priors"]][["psi"]][["v_i"]] <- diag(coef[["v_i"]], n_covar)
    if (object[["model"]][["tvp"]]) {
      object[["priors"]][["psi"]][["v_i"]] <- Matrix::Matrix(object[["priors"]][["psi"]][["v_i"]]) 
      object[["priors"]][["psi"]][["shape"]] <- matrix(coef[["shape"]], n_covar)
      object[["priors"]][["psi"]][["rate"]] <- matrix(coef[["rate"]], n_covar) 
    }
    
    if (use_ssvs_error) {
      object$model$varselect <- "SSVS"
      object[["priors"]][["psi"]][["ssvs"]][["inprior"]] <- matrix(ssvs[["inprior"]], n_covar)
      object[["priors"]][["psi"]][["ssvs"]][["include"]] <- matrix(1:n_covar)
      object[["priors"]][["psi"]][["ssvs"]][["tau0"]] <- matrix(ssvs[["tau"]][1], n_covar)
      object[["priors"]][["psi"]][["ssvs"]][["tau1"]] <- matrix(ssvs[["tau"]][2], n_covar)
    }
    
    if (use_bvs_error) {
      object$model$varselect <- "BVS"
      object[["priors"]][["psi"]][["bvs"]][["inprior"]] <- matrix(bvs[["inprior"]], n_covar)
      object[["priors"]][["psi"]][["bvs"]][["include"]] <- matrix(1:n_covar)
    }
    
  }
  
  # Error term ----
  if (sv) {
    
    object <- .add_priors_sv_helper(object, sigma, k)
    
  } else {
    if (error_prior == "wishart") {
      object[["priors"]][["sigma"]][["type"]] <- "wishart"
      help_df <- sigma[["df"]]
      object[["priors"]][["sigma"]][["df"]] <- NA
      object[["priors"]][["sigma"]][["scale"]] = diag(sigma[["scale"]], k)
    }
    if (error_prior == "gamma") {
      object[["priors"]][["sigma"]][["type"]] <- "gamma"
      help_df <- sigma[["shape"]]
      object[["priors"]][["sigma"]][["shape"]] <- NA
      object[["priors"]][["sigma"]][["rate"]] = matrix(sigma[["rate"]], k)
    }
    
    if (minnesota & !is.null(object[["data"]][["x"]])) {
      # Store LS estimate of variance coviariance matrix for analytical solution
      object[["priors"]][["sigma"]][["sigma_i"]] = minn[["sigma_i"]]
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
    
    # Add rank to degrees of freedom for cointegration model
    if (!is.na(object[["model"]][["rank"]])) {
      help_df <- help_df + r
    }
    
    if (error_prior == "wishart") {
      object[["priors"]][["sigma"]][["df"]] <- help_df
    }
    if (error_prior == "gamma") {
      object[["priors"]][["sigma"]][["shape"]] <- matrix(help_df, k)
    }
  }
  
  return(object)
}