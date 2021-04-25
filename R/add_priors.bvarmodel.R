#' Add Priors to Model
#'
#' Adds prior specifications to a list of models, which was produced by
#' function \code{\link{gen_var}} or \code{\link{gen_vec}}.
#'
#' @param object a list, usually, the output of a call to \code{\link{gen_var}} or \code{\link{gen_vec}}.
#' @param coef a named list of prior specifications for the coefficients of the
#' models. For the default specification all prior means are set to zero and the diagonal elements of
#' the inverse prior variance-covariance matrix are set to 1 for coefficients corresponding to non-deterministic
#' terms. For deterministic coefficients the prior variances are set to 10 via \code{v_i_det = 0.1}.
#' The variances need to be specified as precisions, i.e. as inverses of the variances.
#' For further specifications such as the Minnesota prior see 'Details'.
#' @param coint a named list of prior specifications for cointegration coefficients of
#' country-specific VEC models. See 'Details'.
#' @param sigma a named list of prior specifications for the error variance-covariance matrix
#' of the country models. For the default specification of an inverse Wishart distribution
#' the prior degrees of freedom are set to the number of endogenous variables and
#' the prior variances to 1. See 'Details'.
#' @param ssvs optional; a named list of prior specifications for the SSVS algorithm. See 'Details'.
#' @param bvs optional; a named list of prior specifications for the BVS algorithm. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details Argument \code{coef} can contain the following elements
#' \describe{
#'   \item{\code{v_i}}{a numeric specifying the prior precision of the coefficients. Default is 1.}
#'   \item{\code{v_i_det}}{a numeric specifying the prior precision of coefficients corresponding to deterministic terms. Default is 0.1.}
#'   \item{\code{coint_var}}{a logical specifying whether the prior mean of the first own lag of an
#'   endogenous variable in a VAR model should be set to 1. Default is \code{FALSE}.}
#'   \item{\code{const}}{a numeric or character specifying the prior mean of coefficients, which correspond
#'   to the intercept. If a numeric is provided, all prior means are set to this value.
#'   If \code{const = "mean"}, the means of the series of endogenous variables are used as prior means.
#'   If \code{const = "first"}, the first values of the series of endogenous variables are used as prior means.}
#'   \item{\code{minnesota}}{a list of length 4 containing parameters for the calculation of
#'   the Minnesota prior, where the element names are \code{kappa0}, \code{kappa1}, \code{kappa2} and \code{kappa3}.
#'   For the endogenous variable \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#'   \deqn{ \frac{\kappa_{0}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#'   \deqn{ \frac{\kappa_{0} \kappa_{1}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{0} \kappa_{2}}{(l+1)^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for exogenous variables,}}
#'   \deqn{ \kappa_{0} \kappa_{3} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#'   where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#'   LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#'   
#'   For VEC models the function only provides priors for the non-cointegration part of the model. The
#'   residual standard errors \eqn{\sigma_i} are based on an unrestricted LS regression of the
#'   endogenous variables on the error correction term and the non-cointegration regressors.}
#'   \item{\code{max_var}}{a numeric specifying the maximum prior variance that is allowed for
#'   non-deterministic coefficients.}
#' }
#' If \code{minnesota} is specified, \code{v_i} and \code{v_i_det} are ignored.
#' 
#' Argument \code{coint} can contain the following elements:
#' \describe{
#'   \item{\code{coint_v_i}}{numeric between 0 and 1 specifying the shrinkage of the cointegration space prior. Default is 0.}
#'   \item{\code{coint_p_tau_i}}{numeric of the diagonal elements of the inverse prior matrix of
#' the central location of the cointegration space \eqn{sp(\beta)}. Default is 1.}
#' }
#' 
#' Argument \code{sigma} can contain the following elements:
#' \describe{
#'   \item{\code{df}}{an integer or character specifying the prior degrees of freedom of the error term. Only
#'   used, if the prior is invese Wishart.
#'   Default is \code{"k"}, which indicates the amout of endogenous variables in the respective country model.
#'   \code{"k + 3"} can be used to set the prior to the amount of endogenous variables plus 3. If an integer
#'   is provided, the degrees of freedom are set to this value in all models.
#'   If a VEC model is estimated, the rank \eqn{r} of the cointegration matrix
#'   is automatically added.}
#'   \item{\code{scale}}{a numeric specifying the prior error variance of endogenous variables.
#'   Default is 1.}
#'   \item{\code{shape}}{a numeric or character specifying the prior shape parameter of the error term. Only
#'   used, if the prior is inverse gamma.
#'   Default is \code{"k"}, which indicates the amout of endogenous variables in the respective country model.
#'   \code{"k + 3"} can be used to set the prior to the amount of endogenous variables plus 3. If a numeric
#'   is provided, the shape parameters are set to this value in all models.
#'   If a VEC model is estimated, the rank \eqn{r} of the cointegration matrix
#'   is automatically added.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error term. Only used, if the
#'   prior is inverse gamma.}
#' }
#' \code{df} and \code{scale} must be specified for an inverse Wishart prior. \code{shape} and \code{rate}
#' are required for an inverse gamma prior. For structural models only a gamma prior specification
#' is allowed.
#' 
#' Argument \code{ssvs} can contain the following elements:
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model. Default is 0.5.}
#'   \item{\code{tau}}{a numeric vector of two elements containing the prior standard errors
#'   of restricted variables (\eqn{\tau_0}) as its first element and unrestricted variables (\eqn{\tau_1})
#'   as its second. Default is \code{c(0.05, 10)}.}
#'   \item{\code{semiautomatic}}{an numeric vector of two elements containing the
#'   factors by which the standard errors associated with an unconstrained least squares
#'   estimate of the VAR model are multiplied to obtain the prior standard errors
#'   of restricted (\eqn{\tau_0}) and unrestricted (\eqn{\tau_1}) variables, respectively.
#'   This is the semiautomatic approach described in George et al. (2008).}
#'   \item{\code{covar}}{logical indicating if SSVS should also be applied to the error covariance matrix
#'   as in George et al. (2008).}
#'   \item{\code{exclude_det}}{logical indicating if deterministic terms should be exempted from the SSVS algorithm.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters for the calculation of
#'   the Minnesota-like inclusion priors. See below.}
#' }
#' Either \code{tau} or \code{semiautomatic} must be specified.
#' 
#' The argument \code{bvs} can contain the following elements
#' \describe{
#'   \item{\code{inprior}}{a numeric between 0 and 1 specifying the prior probability
#'   of a variable to be included in the model.}
#'   \item{\code{covar}}{logical indicating if BVS should also be applied to the error covariance matrix.}
#'   \item{\code{exclude_det}}{logical indicating if deterministic terms should be exempted from the BVS algorithm.}
#'   \item{\code{minnesota}}{a numeric vector of length 4 containing parameters for the calculation of
#'   the Minnesota-like inclusion priors. See below.}
#' }
#' 
#' If either \code{ssvs$minnesota} or \code{bvs$minnesota} is specified, prior inclusion probabilites
#' are calculated in a Minnesota-like fashion as
#' \tabular{cl}{
#' \eqn{\frac{\kappa_1}{l}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_2}{l}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_3}{1 + l}} \tab for exogenous variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables, 
#' }
#' for lag \eqn{l} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3}, \eqn{\kappa_4} as the first, second,
#' third and forth element in \code{ssvs$minnesota} or \code{bvs$minnesota}, respectively.
#' 
#' @return A list of country models.
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
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' model <- gen_var(e1, p = 2, deterministic = 2,
#'                  iterations = 100, burnin = 10)
#' 
#' model <- add_priors(model)
#' 
#' @export
add_priors.bvarmodel <- function(object,
                                 coef = list(v_i = 1, v_i_det = 0.1),
                                 coint = list(v_i = 0, p_tau_i = 1),
                                 sigma = list(df = "k", scale = 1),
                                 ssvs = NULL,
                                 bvs = NULL,
                                 ...){
  
  # rm(list = ls()[-which(ls() == "object")]); coef = list(v_i = 1, v_i_det = 0.1); coint = list(v_i = 0, p_tau_i = 1); sigma = list(df = "k", scale = 1); ssvs = NULL; bvs = NULL
  
  only_one_model <- FALSE
  # If only one model is provided, make it compatible with the rest
  if ("data" %in% names(object)) {
    object <- list(object)
    only_one_model <- TRUE
  }
  
  # Checks - Coefficient priors ----
  if (!is.null(coef)) {
    if (!is.null(coef$v_i)) {
      if (coef$v_i < 0) {
        stop("Argument 'v_i' must be at least 0.")
      }  
    } else {
      if (!any(c("minnesota", "ssvs") %in% names(coef))) {
        stop("If 'coef$v_i' is not specified, at least 'coef$minnesota' or 'coef$ssvs' must be specified.")
      }
    }
  }
  
  if (!is.null(coef$const)) {
    if (class(coef$const) == "character") {
      if (!coef$const %in% c("first", "mean")) {
        stop("Invalid specificatin of coef$const.")
      }
    }
  }
  
  if (length(coint) < 2) {
    stop("Argument 'coint' must be at least of length 2.")
  }
  
  # Checks - Error priors ----
  if (length(sigma) < 2) {
    stop("Argument 'sigma' must be at least of length 2.")
  } else {
    error_prior <- NULL
    if (all(c("shape", "rate") %in% names(sigma))) {
      error_prior <- "gamma"
    }
    if (all(c("df", "scale") %in% names(sigma))) {
      error_prior <- "wishart"
    }
    if (is.null(error_prior)) {
      stop("Invalid specification for argument 'sigma'.")
    }
    if (error_prior == "wishart" & any(unlist(lapply(object, function(x) {x$model$structural})))) {
      stop("Structural models may not use a Wishart prior. Consider specifying argument 'sigma$shape' instead.")
    }
    
    if (error_prior == "wishart") {
      if (sigma$df < 0) {
        stop("Argument 'sigma$df' must be at least 0.")
      }
      if (sigma$scale <= 0) {
        stop("Argument 'sigma$scale' must be larger than 0.")
      } 
    }
    if (error_prior == "gamma") {
      if (sigma$shape < 0) {
        stop("Argument 'sigma$shape' must be at least 0.")
      }
      if (sigma$rate <= 0) {
        stop("Argument 'sigma$rate' must be larger than 0.")
      } 
    }
  }
  
  # Check Minnesota ----
  minnesota <- FALSE # Minnesota prior?
  if (!is.null(coef$minnesota)) {
    minnesota <- TRUE
  }
  
  # Check coint VAR ----
  coint_var <- FALSE # Cointegrated VAR?
  if (!is.null(coef$coint_var)) {
    if (coef$coint_var) {
      coint_var <- TRUE 
    }
  }
  
  # Check SSVS ----
  use_ssvs <- FALSE
  use_ssvs_error <- FALSE
  use_ssvs_semi <- FALSE
  if (!is.null(ssvs)) {
    if (is.null(ssvs$inprior)) {
      stop("Argument 'ssvs$inprior' must be specified for SSVS.")
    }
    if (is.null(ssvs$tau) & is.null(ssvs$semiautomatic)) {
      stop("Either argument 'ssvs$tau' or 'ssvs$semiautomatic' must be specified for SSVS.")
    }
    if (is.null(ssvs$exclude_det)) {
      ssvs$exclude_det <- FALSE
    }
    # In case ssvs is specified, check if the semi-automatic approch of 
    # George et al. (2008) should be used
    if (!is.null(ssvs$semiautomatic)) {
      use_ssvs_semi <- TRUE
    }
    
    use_ssvs <- TRUE
    if (minnesota) {
      minnesota <- FALSE
      warning("Minnesota prior specification overwritten by SSVS.")
    }
    
    if (!is.null(ssvs$covar)) {
      if (ssvs$covar) {
        if (error_prior == "wishart") {
          stop("If SSVS should be applied to error covariances, argument 'sigma$shape' must be specified.")
        }
        use_ssvs_error <- TRUE 
      }
      if (is.null(ssvs$tau)) {
        stop("If SSVS should be applied to error covariances, argument 'ssvs$tau' must be specified.")
      }
    }
  }
  
  # Prior a la Korobilis 2013
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
        if (error_prior == "wishart") {
          stop("If BVS should be applied to error covariances, argument 'sigma$shape' must be specified.")
        }
        use_bvs_error <- TRUE 
      }
    }
    if (coef$v_i == 0 | (coef$v_i_det == 0 & !bvs$exclude_det)) {
      warning("Using BVS with an uninformative prior is not recommended.")
    }
  }
  
  if (use_ssvs & use_bvs) {
    stop("SSVS and BVS cannot be applied at the same time.")
  }
  
  if (error_prior == "wishart" & (use_ssvs_error | use_bvs_error)) {
    stop("Wishart prior not allowed when BVS or SSVS are applied to covariance matrix.")
  }
  
  # Generate priors for each country ----
  for (i in 1:length(object)) {
    
    # Get model specs to obtain total number of coeffs
    k <- length(object[[i]]$model$endogen$variables)
    p <- object[[i]]$model$endogen$lags
    
    if (k == 1 & (use_ssvs_error | use_bvs_error)) {
      stop("BVS or SSVS cannot be applied to covarianc matrix when there is only one endogenous variable.")
    } 
    
    use_exo <- FALSE
    if (!is.null(object[[i]]$model$exogen)) {
      use_exo <- TRUE
      m <- length(object[[i]]$model$exogen$variables)
      s <- object[[i]]$model$exogen$lags
    } else {
      s <- 0
      m <- 0
    }
    
    if (object[[i]]$model$type == "VAR") {
      # Add a lag to exogenous vars to include
      # contemporary variables
      if (use_exo) {
        s <- s + 1
      }
    }
    if (object[[i]]$model$type == "VEC") {
      # Substract lag from domestic model for VEC
      p <- p - 1
    }
    
    # Total # of non-deterministic coefficients
    tot_par <- k * (k * p + m * s)
    
    structural <- object[[i]]$model$structural
    n_struct <- 0
    if (structural & k > 1) {
      n_struct <- (k - 1) * k / 2
      tot_par <- tot_par + n_struct
    }
    
    # Add number of non-cointegration deterministic terms
    n_det <- 0
    if (!is.null(object[[i]]$model$deterministic)){
      if (object[[i]]$model$type == "VAR") {
        n_det <- length(object[[i]]$model$deterministic) * k
      }
      if (object[[i]]$model$type == "VEC") {
        if (!is.null(object[[i]]$model$deterministic$unrestricted)) {
          n_det <- length(object[[i]]$model$deterministic$unrestricted) * k
        }
      }
      tot_par <- tot_par + n_det
    }
    
    #### Cointegration (constant) ----
    if (object[[i]]$model$type == "VEC") {
      n_ect <- k * (k + m)
      if (!is.null(object[[i]]$model$deterministic$restricted)) {
        n_ect <- n_ect + k * length(object[[i]]$model$deterministic$restricted)
      }
      
      r_temp <- object[[i]]$model$rank
      if (r_temp > 0) {
        object[[i]]$priors$cointegration <- list(v_i = coint$v_i,
                                                 p_tau_i = diag(coint$p_tau_i, n_ect / k))
      } else {
        if (r_temp == 0 & tot_par == 0) {
          warning("Model with zero cointegration rank and no non-cointegration regressors is skipped.")
        }
      }
    }
    
    #### Non-cointegration (constant) ####
    # Generate prior matrices ----
    if (tot_par > 0) {
      
      # Zero prior means
      mu <- matrix(rep(0, tot_par - n_struct), k)
      
      # Add 1 to first own lags for coint VAR
      if (coint_var & object[[i]]$model$type == "VAR") {
        mu[1:k, 1:k] <- diag(1, k)
      }
      
      # Prior for intercept terms
      if (n_det > 0) {
        
        if (!is.null(coef$const))  {
          
          if (object[[i]]$model$type == "VAR") {
            pos <- which(dimnames(object[[i]]$data$Z)[[2]] == "const")
          }
          if (object[[i]]$model$type == "VEC") {
            pos <- which(dimnames(object[[i]]$data$X)[[2]] == "const")
          }
          
          if (length(pos) == 1) {
            if (class(coef$const) == "character") {
              if (coef$const == "first") {
                mu[, pos] <- object[[i]]$data$Y[1, ]
              }
              if (coef$const == "mean") {
                mu[, pos] <- colMeans(object[[i]]$data$Y)
              }
            }
            if (class(coef$const) == "numeric") {
              mu[, pos] <- coef$const
            } 
          }
        }
      }
      
      mu <- matrix(mu)
      if (structural) {
        mu <- rbind(mu, matrix(0, n_struct))
      }
      
      if (object[[i]]$model$type == "VAR") {
        object[[i]]$priors$coefficients <- list(mu = mu)
      }
      if (object[[i]]$model$type == "VEC") {
        object[[i]]$priors$noncointegration <- list(mu = mu)
      }
      
      # Prior covariances
      if (minnesota) {
        # Minnesota prior ----
        minn <- minnesota_prior(object = object[[i]],
                                kappa0 = coef$minnesota$kappa0,
                                kappa1 = coef$minnesota$kappa1,
                                kappa2 = coef$minnesota$kappa2,
                                kappa3 = coef$minnesota$kappa3,
                                max_var = NULL,
                                coint_var = FALSE,
                                sigma = "AR")
        
        if (object[[i]]$model$type == "VAR") {
          object[[i]]$priors$coefficients$v_i <- minn$v_i 
        }
        if (object[[i]]$model$type == "VEC") {
          object[[i]]$priors$noncointegration$v_i <- minn$v_i 
        }
      }
      
      # SSVS prior ----
      if (use_ssvs) {
        
        ssvs_temp <- ssvs_prior(object[[i]], tau = ssvs$tau, semiautomatic = ssvs$semiautomatic)
        temp <- inclusion_prior(object[[i]], prob = ssvs$inprior, exclude_deterministics = ssvs$exclude_det,
                                minnesota_like = !is.null(ssvs$minnesota), kappa = ssvs$minnesota)
        object[[i]]$model$varselect <- "SSVS"
        
        if (object[[i]]$model$type == "VAR") {
          object[[i]]$priors$coefficients$v_i <- diag(1 / ssvs_temp$tau1[, 1]^2, tot_par)
          object[[i]]$priors$coefficients$ssvs$inprior <- temp$prior
          object[[i]]$priors$coefficients$ssvs$include <- temp$include
          object[[i]]$priors$coefficients$ssvs$tau0 <- ssvs_temp$tau0
          object[[i]]$priors$coefficients$ssvs$tau1 <- ssvs_temp$tau1
        }
        
        if (object[[i]]$model$type == "VEC") {
          if (r_temp > 0) {
            # Drop ECT output from inclusion and SSVS output
            pos_omit <- 1:(k * NCOL(object[[i]]$data$W))
            ssvs_temp$tau0 <- matrix(ssvs_temp$tau0[-pos_omit,])
            ssvs_temp$tau1 <- matrix(ssvs_temp$tau1[-pos_omit,])
            temp$prior <- matrix(temp$prior[-pos_omit,])
            pos_incl <- which(temp$include %in% pos_omit)
            if (length(pos_incl) > 0) {
              temp$include <- matrix(temp$include[-pos_incl])
            }
            temp$include <- temp$include - (temp$include[1] - 1)
            rm(pos_omit)
            rm(pos_incl)
          }
          
          object[[i]]$priors$noncointegration$v_i <- diag(1 / ssvs_temp$tau1[, 1]^2, tot_par)
          object[[i]]$priors$noncointegration$ssvs$inprior <- temp$prior
          object[[i]]$priors$noncointegration$ssvs$include <- temp$include
          object[[i]]$priors$noncointegration$ssvs$tau0 <- ssvs_temp$tau0
          object[[i]]$priors$noncointegration$ssvs$tau1 <- ssvs_temp$tau1
        }
        rm(temp)
      }
      
      # Regular prior ----
      if (!minnesota & !use_ssvs) {
        v_i <- diag(coef$v_i, tot_par)
        # Add priors for deterministic terms if they were specified
        if (n_det > 0 & !is.null(coef$v_i_det)) {
          diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef$v_i_det
        }
        if (object[[i]]$model$type == "VAR") {
          object[[i]]$priors$coefficients$v_i <- v_i   
        }
        if (object[[i]]$model$type == "VEC") {
          object[[i]]$priors$noncointegration$v_i <- v_i
        }
      }
      
      if (use_bvs) {
        object[[i]]$model$varselect <- "BVS"
        temp <- inclusion_prior(object[[i]], prob = bvs$inprior, exclude_deterministics = bvs$exclude_det,
                                minnesota_like = !is.null(bvs$minnesota), kappa = bvs$minnesota)
        
        if (object[[i]]$model$type == "VAR") {
          object[[i]]$priors$coefficients$bvs$inprior <- temp$prior
          object[[i]]$priors$coefficients$bvs$include <- temp$include
        }
        if (object[[i]]$model$type == "VEC") {
          # Drop ECT output from inclusion and SSVS output
          if (r_temp > 0) {
            pos_omit <- 1:(k * NCOL(object[[i]]$data$W)) # Positions of Pi values
            temp$prior <- matrix(temp$prior[-pos_omit,]) # Incl prior w/o Pi values
            pos_incl <- which(temp$include %in% pos_omit) # Drop Pi values
            if (length(pos_incl) > 0) {
              temp$include <- matrix(temp$include[-pos_incl])
            }
            temp$include <- temp$include - (temp$include[1] - 1)
            rm(pos_omit)
            rm(pos_incl) 
          }
          object[[i]]$priors$noncointegration$bvs$inprior <- temp$prior
          object[[i]]$priors$noncointegration$bvs$include <- temp$include
        }
      }
    }
    
    # Error term ----
    if (error_prior == "wishart") {
      object[[i]]$priors$sigma$type <- "wishart"
      help_df <- sigma$df
      object[[i]]$priors$sigma$df <- NA
      object[[i]]$priors$sigma$scale = diag(sigma$scale, k)
    }
    if (error_prior == "gamma") {
      object[[i]]$priors$sigma$type <- "gamma"
      help_df <- sigma$shape
      object[[i]]$priors$sigma$shape <- NA
      object[[i]]$priors$sigma$rate = diag(sigma$rate, k)
    }
    if (minnesota) {
      # Store LS estimate of variance coviariance matrix for analytical solution
      object[[i]]$priors$sigma$sigma_i = minn$sigma_i
    }
    
    if (class(help_df) == "character") {
      if (grepl("k", help_df)) {
        # Transform character specification to expression and evaluate
        help_df <- eval(parse(text = help_df))
      } else {
        stop("Use letter 'k' in sigma$df to indicate the number of endogenous variables.")
      }
    }
    
    if (help_df < 0) {
      stop("Current specification implies a negative prior degree of\nfreedom or shape parameter of the error term.")
    }
    
    if (object[[i]]$model$type == "VEC") {
      # Add rank to degrees of freedom for cointegration model
      if (!is.na(object[[i]]$model$rank)) {
        help_df <- help_df + r_temp
      }
    }
    
    if (error_prior == "wishart") {
      object[[i]]$priors$sigma$df <- help_df
    }
    if (error_prior == "gamma") {
      object[[i]]$priors$sigma$shape <- help_df
      if (use_ssvs_error) {
        if (object[[i]]$model$structural) {
          stop("Not allowed to apply SSVS to covariances when a structural model is estimated.")
        }
        object[[i]]$priors$sigma$v_i <- diag(1 / ssvs$tau[2]^2, k * (k - 1) / 2)
        object[[i]]$priors$sigma$ssvs$inprior <- matrix(ssvs$inprior, k * (k - 1) / 2)
        object[[i]]$priors$sigma$ssvs$tau0 <- matrix(ssvs$tau[1], k * (k - 1) / 2)
        object[[i]]$priors$sigma$ssvs$tau1 <- matrix(ssvs$tau[2], k * (k - 1) / 2)
      }
      
      if (use_bvs_error) {
        if (object[[i]]$model$structural) {
          stop("Not allowed to apply BVS to covariances when a structural model is estimated.")
        }
        object[[i]]$priors$sigma$mu <- matrix(0, k * (k - 1) / 2)
        object[[i]]$priors$sigma$v_i <- diag(coef$v_i, k * (k - 1) / 2)
        object[[i]]$priors$sigma$bvs$inprior <- matrix(bvs$inprior, k * (k - 1) / 2)
      }
    }
  }
  
  if (only_one_model) {
    object <- object[[1]]
  }
  
  return(object)
}