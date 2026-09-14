#' Add Priors to Bayesian Models
#'
#' Adds prior specifications to a list of models, which was produced by
#' function \code{\link{create_bvecmodel}}.
#'
#' @param object a list of class 'bvecmodel'.
#' @param coef a named list of prior specifications for coefficients that do not
#' determine the cointegration space. It has no default and must contain at least
#' \code{v_i} or \code{minnesota}. Variances are specified as precisions, i.e. as
#' inverses of the variances. See 'Details'.
#' @param coint a named list of prior specifications for coefficients determining the
#' cointegration space of VEC models. It has no default. See 'Details'.
#' @param sigma a named list of prior specifications for the error term. It has
#' no default, and the elements it must contain depend on argument \code{error}
#' of \code{\link{create_bvecmodel}}. See 'Details'.
#' @param varsel a named list of prior specifications for the variable
#' selection algorithm. Required if the model was created with
#' \code{varsel = "ssvs"} or \code{"bvs"}, and not allowed otherwise. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details None of the arguments \code{coef}, \code{coint}, \code{sigma} and
#' \code{varsel} provides default hyperparameters: every value that a model needs
#' must be given in the list it belongs to. A missing required element raises an
#' error, as does an element that is not listed below.
#'
#' Argument \code{coef} can contain the following elements:
#' \describe{
#'   \item{\code{v_i}}{a non-negative numeric specifying the prior precision of the coefficients,
#'   where 0 gives an uninformative prior. Required unless \code{minnesota} is given. It is also
#'   required together with \code{minnesota} if \code{error} is \code{"gamma+covar"} or
#'   \code{"sv+covar"}. The precisions of the other coefficients are taken from \code{minnesota}
#'   if it is given, and from \code{varsel} for SSVS.}
#'   \item{\code{v_i_det}}{a numeric specifying the prior precision of coefficients
#'   corresponding to deterministic terms. If it is not given, \code{v_i} is used.
#'   Not used if \code{minnesota} is given or SSVS is applied.}
#'   \item{\code{const}}{a numeric or character specifying the prior mean of coefficients, which correspond
#'   to the intercept. If a numeric is provided, all prior means are set to this value.
#'   If \code{coef$const = "mean"}, the mean of the respective endogenous variable is used as prior mean.
#'   If \code{coef$const = "first"}, the first values of the respective endogenous variable is used as prior mean.}
#'   \item{\code{minnesota}}{a named list containing the parameters for the calculation of
#'   the Minnesota prior. It must contain \code{kappa1}, \code{kappa2} and \code{kappa4}, and
#'   \code{kappa3} if the model has exogenous variables.
#'   For the endogenous variable \eqn{i} the prior variance of the \eqn{l}th lag of regressor \eqn{j} is obtained as
#'   \deqn{ \frac{\kappa_{1}}{l^2} \textrm{ for own lags of endogenous variables,}} 
#'   \deqn{ \frac{\kappa_{1} \kappa_{2}}{l^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for endogenous variables other than own lags,}}
#'   \deqn{ \frac{\kappa_{1} \kappa_{3}}{(l+1)^2} \frac{\sigma_{i}^2}{\sigma_{j}^2} \textrm{ for exogenous variables,}}
#'   \deqn{ \kappa_{1} \kappa_{4} \sigma_{i}^2 \textrm{ for deterministic terms,}}
#'   where \eqn{\sigma_{i}} is the residual standard deviation of variable \eqn{i} of an unrestricted
#'   LS estimate. For exogenous variables \eqn{\sigma_{i}} is the sample standard deviation.
#'   If the model does not contain exogenous variables, \code{kappa3} will be ignored.
#'   The function only provides priors for the non-cointegration part of the model. However,
#'   the residual standard errors \eqn{\sigma_i} are based on an unrestricted LS regression of the
#'   endogenous variables on the error correction term and the non-cointegration regressors.}
#'   \item{\code{max_var}}{a positive numeric specifying the maximum prior variance of the
#'   coefficients of non-deterministic variables in the Minnesota prior. Larger prior variances
#'   are set to this value. Only used if \code{minnesota} is given.}
#'   \item{\code{shape}}{a numeric specifying the prior shape parameter of the error variances of the
#'   state equation. Required for models with time varying parameters and not used otherwise.}
#'   \item{\code{rate}}{a numeric specifying the prior rate parameter of the error variances of the
#'   state equation. Required for models with time varying parameters and not used otherwise.}
#'   \item{\code{rate_det}}{a numeric specifying the prior rate parameter of the error variances of the
#'   state equation for coefficients, which correspond to deterministic terms. If it is not given,
#'   \code{rate} is used. Only used for models with time varying parameters.}
#' }
#' 
#' Argument \code{coint} specifies the prior on the cointegration space. Its
#' elements, and the priors they give, are described in section 'Prior on the
#' cointegration space' below, which is shared with \code{\link{cointspace_prior}}.
#'
#' Argument \code{sigma} must contain the elements that belong to the
#' \code{error} of the model:
#' \itemize{
#'   \item \code{"wishart"}: \code{df} and \code{scale}. Not available for structural models.
#'   \item \code{"gamma"} and \code{"gamma+covar"}: \code{shape} and \code{rate}.
#'   \item \code{"sv"} and \code{"sv+covar"}: \code{mu}, \code{v_i}, \code{shape}, \code{rate},
#'   \code{state_variance} and \code{offset}.
#' }
#' The elements are
#' \describe{
#'   \item{\code{df}}{a non-negative integer, or a character expression in \code{k}, the number of
#'   endogenous variables, such as \code{"k"} or \code{"k + 3"}, specifying the prior degrees of
#'   freedom of the inverse Wishart prior. The samplers add the rank \eqn{r} of the cointegration
#'   matrix to the posterior degrees of freedom, as the prior of the loadings requires.}
#'   \item{\code{scale}}{a positive numeric specifying the prior error variance of the endogenous
#'   variables in the inverse Wishart prior.}
#'   \item{\code{shape}}{for \code{"gamma"} and \code{"gamma+covar"} a non-negative numeric, or a
#'   character expression in \code{k} as for \code{df}, specifying the prior shape parameter of the
#'   error variances. For models with stochastic
#'   volatility a numeric specifying the prior shape parameter of the error variance of the state
#'   equation of the log-volatilities.}
#'   \item{\code{rate}}{a positive numeric specifying the prior rate parameter that corresponds to
#'   \code{shape}.}
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
#' \eqn{\frac{\kappa_1}{l}} \tab for own lags of endogenous variables, \cr
#' \eqn{\frac{\kappa_2}{l}} \tab for other endogenous variables, \cr
#' \eqn{\frac{\kappa_3}{1 + l}} \tab for exogenous variables, \cr
#' \eqn{\kappa_{4}} \tab for deterministic variables, 
#' }
#' for lag \eqn{l} with \eqn{\kappa_1}, \eqn{\kappa_2}, \eqn{\kappa_3},
#' \eqn{\kappa_4} as the first, second, third and forth element in
#' \code{varsel$minnesota}, respectively.
#' 
#' @inheritSection cointspace_prior Prior on the cointegration space
#'
#' @return The object in \code{object} with the element \code{priors} added, a list with
#' \describe{
#'   \item{\code{beta}}{the prior of the cointegration space with \code{type}
#'   \code{"cointspace"}: \code{v_inv} and \code{p_tau_inv} for constant cointegration
#'   parameters, or \code{rho}, \code{mu} and \code{v_inv} of the state equation for time
#'   varying ones, together with the elements added by \code{p_tau_i = "ml"} or a
#'   uniform prior on \eqn{\rho}.}
#'   \item{\code{a}}{the prior of the loadings and the remaining coefficients, with the
#'   same elements as for a VAR model in \code{\link{add_priors.bvarmodel}}.}
#'   \item{\code{psi}, \code{u_sigma}}{the priors of the error covariance coefficients
#'   and error variances, as for a VAR model.}
#' }
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
#' Johansen, S. (1995). \emph{Likelihood-based inference in cointegrated vector
#' autoregressive models}. Oxford: Oxford University Press.
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
#' e6 <- e6 * 100
#' 
#' # Generate model
#' model <- create_bvecmodel(e6, p = 1, r = 1, const = "restricted",
#'                           iterations = 10, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     coint = list(v_i = 0, p_tau_i = 1),
#'                     sigma = list(df = "k", scale = 1))
#'
#' # The same model with time varying parameters. The cointegration vectors are
#' # then a state path rather than a draw from a cointegration space prior, so
#' # 'coint' takes the autocorrelation of that path instead of its shrinkage and
#' # central location, and 'coef' takes the prior of the state error variances.
#' model <- create_bvecmodel(e6, p = 4, r = 1, tvp = TRUE,
#'                           const = "unrestricted",
#'                           seasonal = "unrestricted",
#'                           iterations = 10, burnin = 10)
#'
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10,
#'                                 shape = 3, rate = 0.0001),
#'                     coint = list(rho = 0.999),
#'                     sigma = list(df = "k", scale = 1))
#'
#' @family model set-up
#' @export
#' @method add_priors bvecmodel
add_priors.bvecmodel <- function(object,
                                 coef,
                                 coint,
                                 sigma,
                                 varsel = NULL,
                                 ...){
  
  # Input checks ----
  
  ## coefficients ----
  if (!is.null(coef)) {
    .add_priors_check_coef(object, coef)

    if (is.null(coef[["v_i_det"]])) {
      coef[["v_i_det"]] <- coef[["v_i"]]
    }
  }

  ## cointegration ----
  # Checked and built by cointspace_prior(), which packages with VEC models of
  # their own layout call as well. It is stored further down, with the others.
  beta_prior <- cointspace_prior(object, coint)
  
  ## sigma ----
  error_prior <- .add_priors_check_sigma(object, sigma)
  
  ## Minnesota ----
  minnesota <- FALSE # Minnesota prior?
  if (!is.null(coef[["minnesota"]])) {
    minnesota <- TRUE
  }
  
  if (!is.null(varsel) & object[["model"]][["varsel"]] == "none") {
    stop("Argument 'varsel' was provided, but according to argument 'object$model$varsel' no variable selection algorithm should be used.")
  }
  
  ## SSVS ----
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
    # In case varsel is specified, check if the semi-automatic approach of 
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
  
  ## BVS ----
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
    # A Minnesota prior is informative and needs no v_i
    if (!minnesota && (coef[["v_i"]] == 0 | (coef[["v_i_det"]] == 0 & !varsel[["exclude_det"]]))) {
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
  
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  use_exo <- m > 0
  
  # Substract lag from domestic model for VEC
  p <- p - 1
  
  # Total # of non-deterministic coefficients
  n_alpha <- k * r
  n_gamma <- k * (k * p)
  n_upsilon <- k * (m * s)
  
  # Add number of unrestricted deterministic terms
  n_det <- object[["model"]][["n"]] * k
  
  tot_par <- n_alpha + n_gamma + n_upsilon + n_det
  
  covar <- object[["model"]][["error"]] %in% c("gamma+covar", "sv+covar")
  structural <- object[["model"]][["structural"]]
  if (covar & structural) {
    stop("Error covariances and structural coefficients cannot be estimated at the same time.")
  }

  # A constant coefficient sampler selects over one set of coefficients or over
  # both: it reads a single selection scheme for the whole model, so a
  # covariance block it is given goes into the selection with the rest. Only the
  # time varying samplers take the covariance block's scheme separately, which
  # is why the same call is allowed there. Left to run, this combination fails
  # inside the sampler on a prior it was never given.
  if ((use_ssvs | use_bvs) & covar & !varsel_covar &
      !object[["model"]][["tvp"]] & k > 1) {
    stop("Variable selection cannot be restricted to the coefficients when the ",
         "model has an error covariance block and constant coefficients: this ",
         "sampler applies one selection scheme to both. Set 'varsel$covar' to ",
         "TRUE to select over the covariances as well, drop the covariances with ",
         "an 'error' of \"gamma\" or \"sv\", or use a time varying model, where ",
         "the two blocks can differ.")
  }

  sv <- object[["model"]][["error"]] %in% c("sv", "sv+covar")
  n_struct <- 0
  n_z <- NCOL(object[["data"]][["train"]][["z"]])
  if (object[["model"]][["rank"]] > 0) {
    n_w <- NCOL(object[["data"]][["train"]][["w"]])
  }
  
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
  
  #### Cointegration ----
  # Built by cointspace_prior() with the checks above, and stored here so that
  # the elements of 'priors' keep their order.
  if (r > 0) {
    object[["priors"]][["beta"]] <- beta_prior
  }
  
  # Non-cointegration ----
  # Generate prior matrices ----
  if (tot_par > 0) {
    
    # Zero prior means
    mu <- matrix(rep(0, tot_par - n_struct), k)
    
    # Prior for intercept terms
    if (n_det > 0) {
      
      if (!is.null(coef[["const"]]))  {
        
        # The columns of 'mu' are the r columns of alpha followed by the
        # regressors in 'x', so the offset belongs on the position and must not
        # be added a second time when indexing.
        pos <- which(dimnames(object[["data"]][["train"]][["x"]])[[2]] == "const") + r

        if (length(pos) == 1) {
          if ("character" %in% class(coef[["const"]])) {
            if (coef[["const"]] == "first") {
              mu[, pos] <- object[["data"]][["train"]][["y"]][1,]
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
    
    # Prior covariances
    if (minnesota) {
      # Minnesota prior ----
      minn <- minnesota_prior(object = object,
                              kappa1 = coef[["minnesota"]][["kappa1"]],
                              kappa2 = coef[["minnesota"]][["kappa2"]],
                              kappa3 = coef[["minnesota"]][["kappa3"]],
                              kappa4 = coef[["minnesota"]][["kappa4"]],
                              max_var = coef[["max_var"]],
                              sigma = "AR")
      
      object[["priors"]][["a"]][["v_inv"]] <- minn[["v_i"]]
    }
    
    # SSVS prior ----
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
    
    # Regular prior ----
    if (!minnesota & !use_ssvs) {
      
      if (object[["model"]][["tvp"]]) {

        # For a TVP model this is the prior precision of the state before the
        # sample rather than of a constant coefficient.
        v_i <- diag(coef[["v_i"]], tot_par)

        # The loadings carry the compensating scale of the cointegration space.
        # beta_t has stationary variance 1 / (1 - rho^2), which for a rho just
        # below one is large, and only the product alpha beta' is identified --
        # so alpha's prior variance is shrunk by the same factor, leaving the
        # product on the scale coef$v_i asks for.
        if (r > 0) {
          diag(v_i)[1:n_alpha] <- 1 / (1 - coint[["rho"]] * coint[["rho"]])
        }
        if (n_det > 0 & !is.null(coef[["v_i_det"]])) {
          diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef[["v_i_det"]]
        }
        object[["priors"]][["a"]][["shape"]] <- matrix(coef[["shape"]], tot_par)
        object[["priors"]][["a"]][["rate"]] <- matrix(coef[["rate"]], tot_par)
        if (n_det > 0 & !is.null(coef[["rate_det"]])) {
          object[["priors"]][["a"]][["rate"]][tot_par - n_struct - n_det + 1:n_det, ] <- coef[["rate_det"]]
        }
      } else {
        v_i <- diag(coef[["v_i"]], tot_par)
        # Add priors for deterministic terms if they were specified
        if (n_det > 0 & !is.null(coef[["v_i_det"]])) {
          diag(v_i)[tot_par - n_struct - n_det + 1:n_det] <- coef[["v_i_det"]]
        }
      }
      object[["priors"]][["a"]][["v_inv"]] <- v_i
    }
    
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
  
  # Error term ----
  if (sv) {
    
    object <- .add_priors_sv_helper(object, sigma, k)
    
  } else {
    pos <- "u_sigma"
    if (error_prior == "wishart") {
      object[["priors"]][[pos]][["type"]] <- "wishart"
      help_df <- sigma[["df"]]
      object[["priors"]][[pos]][["df"]] <- NA_real_
      object[["priors"]][[pos]][["scale"]] = diag(sigma[["scale"]], k)
    }
    if (error_prior == "gamma") {
      object[["priors"]][[pos]][["type"]] <- "gamma"
      help_df <- sigma[["shape"]]
      object[["priors"]][[pos]][["shape"]] <- NA_real_
      object[["priors"]][[pos]][["rate"]] = matrix(sigma[["rate"]], k)
    }
    
    if (minnesota & !is.null(object[["data"]][["train"]][["x"]])) {
      # Store LS estimate of variance coviariance matrix for analytical solution
      object[["priors"]][[pos]][["u_sigma_inv"]] = minn[["sigma_inv"]]
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
    
    # Stored as given. The constant VECs with a Wishart prior add the rank to
    # the posterior degrees of freedom themselves, for the prior of the loadings
    # given the error covariance (Koop et al., 2010, eq. 8); adding it here too
    # counted it twice, and neither the gamma shape nor the time varying models
    # call for it.
    if (error_prior == "wishart") {
      object[["priors"]][[pos]][["df"]] <- help_df
    }
    if (error_prior == "gamma") {
      object[["priors"]][[pos]][["shape"]] <- matrix(help_df, k)
    }
  }
  
  return(object)
}
