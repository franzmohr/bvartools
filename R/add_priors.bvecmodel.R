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
#' error, as does an element of \code{coef} or \code{varsel} that is not listed
#' below. Elements of \code{coint} and \code{sigma} that are not listed below are
#' ignored without a message, so a misspelt name there leaves the setting it was
#' meant for unapplied.
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
#' Argument \code{coint} can contain the following elements. Which of them are
#' required depends on whether the cointegration vectors are constant or time
#' varying, since the two are given different kinds of prior: a matric-variate
#' prior on the cointegration space in the first case, and a state equation in
#' the second.
#' \describe{
#'   \item{\code{v_i}}{non-negative numeric specifying the shrinkage of the cointegration space prior,
#'   or \code{"ml"}. See below. Required for models with constant cointegration parameters and not
#'   used otherwise.}
#'   \item{\code{p_tau_i}}{the inverse of the matrix \eqn{P_\tau}, which determines the
#'   central location of the cointegration space \eqn{sp(\beta)}. Either a numeric of
#'   its diagonal elements, a full symmetric matrix, or \code{"ml"}. See below.
#'   Required for models with constant cointegration parameters and not used otherwise.}
#'   \item{\code{weight}}{positive numeric specifying the weight of a prior centred on the
#'   maximum likelihood estimate in units of the information of the sample. Default is 1.
#'   Only used if \code{p_tau_i = "ml"}.}
#'   \item{\code{rho}}{a numeric specifying the autocorrelation coefficient
#'   of the state equation of \eqn{\beta}. It must be smaller than 1.
#'   Required for models with time varying cointegration parameters and not used otherwise.
#'   If \code{rho_min} and \code{rho_max} are given as well, this is the value the
#'   chain starts \eqn{\rho} at rather than the value it keeps, and it must lie
#'   between them.}
#'   \item{\code{rho_min}, \code{rho_max}}{numerics specifying the support of a
#'   uniform prior on \eqn{\rho}, which makes it a drawn parameter rather than a
#'   fixed hyperparameter. Optional, and either both or neither; they must
#'   satisfy \eqn{0 < }\code{rho_min}\eqn{ < }\code{rho_max}\eqn{ \le 1}.
#'   Koop et al. (2011) use \eqn{(0.999, 1)}.
#'   Only used for models with time varying cointegration parameters.}
#' }
#' For a model with constant cointegration parameters the prior is that of
#' Koop et al. (2010). The sampler uses \code{v_i} and \code{p_tau_i} only through
#' their product, so with \code{v_i = 0} the prior on the cointegration space is
#' uniform whatever \code{p_tau_i} is. An informative prior on the space
#' therefore needs a positive \code{v_i}, which also shrinks the loadings: for
#' \eqn{\beta} close to the centre of the space they have prior
#' \eqn{N(0, \Sigma / v)}.
#'
#' With \code{p_tau_i = "ml"} the prior is centred on the space spanned by
#' Johansen's (1995) maximum likelihood estimate \eqn{\hat{\beta}}, computed from
#' the error correction term as it is stored in the model -- so call
#' \code{\link{scale_error_correction}} first if the series should be scaled.
#' Let \eqn{H} be an orthonormal basis of that space, \eqn{H_\perp} one of its
#' orthogonal complement and \eqn{\beta = H + H_\perp \delta}. The matrix
#' \deqn{P_\tau^{-1} = H H^{\prime} + H_\perp T^{-1} H_\perp^{\prime}, \quad
#' T = \frac{v}{w} (H_\perp^{\prime} S_{11} H_\perp)^{-1},}
#' where \eqn{S_{11}} is the cross product of the residuals of a regression of the
#' error correction term on the short-run regressors and \eqn{w} is \code{weight},
#' makes the prior of \eqn{\delta} given \eqn{\alpha}
#' \eqn{N(0, (\alpha^{\prime} \Sigma^{-1} \alpha)^{-1} \otimes w^{-1} (H_\perp^{\prime} S_{11} H_\perp)^{-1})}.
#' This is the asymptotic distribution of the maximum likelihood estimator
#' with its precision multiplied by \eqn{w}: \code{weight = 1} adds as much
#' information about the space as the sample itself holds. Since the centre is
#' estimated from the same sample, the posterior then overstates the precision
#' of \eqn{sp(\beta)}. Eigenvalues of \eqn{T} are capped at one, which is the
#' uniform prior.
#'
#' With \code{v_i = "ml"} the shrinkage is chosen so that the prior mean of
#' \eqn{tr(\alpha^{\prime} \Sigma^{-1} \alpha)} is 100 times its maximum
#' likelihood estimate. It can be combined with a numeric \code{p_tau_i}.
#'
#' For a model with time varying cointegration parameters, \code{p_tau_i = "ml"}
#' centres the marginal prior of the cointegration space on Johansen's estimate
#' as well, using the informative marginal prior of Koop et al. (2011, working
#' paper version). The state equation becomes
#' \eqn{\beta_t = \rho (I_r \otimes P_\tau) \beta_{t-1} + \eta_t} with
#' \eqn{P_\tau = H H^{\prime} + H_\perp T H_\perp^{\prime}}: the part of
#' \eqn{\beta_t} along \eqn{sp(H)} keeps \eqn{\rho}, the part off it decays
#' faster, and the mode of the marginal distribution of \eqn{sp(\beta_t)} is
#' \eqn{sp(H)} in every period. \eqn{T} is chosen so that the prior spread of the
#' tilt of \eqn{\beta_t} away from \eqn{sp(H)} in a period, approximately
#' \eqn{T^* = (1 - \rho^2)(I - \rho^2 T^2)^{-1}}, is the sampling variance of
#' Johansen's estimator with its precision multiplied by \code{weight}. The
#' state before the sample is given the stationary distribution the transition
#' implies. The transition is stored in \code{object$priors$beta$p_tau}.
#'
#' \eqn{\rho} limits how informative this prior can be: even \eqn{T = 0} leaves
#' the tilt a spread of \eqn{1 - \rho^2} per period, and a \code{weight} asking
#' for more is floored there with a warning.
#'
#' A larger \code{weight} does not always give a tighter posterior. It narrows
#' the prior spread of the tilt in each period by lowering \eqn{T}, but \eqn{T}
#' is also how much of the tilt carries over from one period to the next. As
#' \eqn{T} approaches zero the tilt of each period becomes independent of the
#' last, the data of neighbouring periods stop informing it, and the posterior
#' spread of \eqn{sp(\beta_t)} can widen again. On data set E6 with
#' \eqn{\rho = 0.999}, for example, \code{weight = 1}, which gives
#' \eqn{T \approx 0.44}, produced a tighter posterior than \code{weight = 100},
#' which reaches \eqn{T = 0}. The useful range of \code{weight} is the one that
#' keeps \eqn{T} clearly above zero, and a \eqn{\rho} closer to one widens it.
#'
#' A \code{weight} so small that
#' \eqn{T} would be the identity in every direction leaves the noninformative
#' prior unchanged. Both \eqn{P_\tau} and the prior on the state before the
#' sample are computed at \code{coint$rho}, and do not follow the draw if
#' \eqn{\rho} is drawn. \code{coint$v_i} has no counterpart for these models.
#'
#' For a model with time varying cointegration parameters the state equation is
#' \eqn{\beta_t = \rho \beta_{t-1} + \eta_t} with \eqn{\eta_t \sim N(0, I)}, and
#' the prior on the state before the sample is that equation's own stationary
#' distribution, \eqn{N(0, I / (1 - \rho^2))}. This is what makes the prior
#' proper, so a \eqn{\rho} close to one is intended and one further from it draws
#' a warning. The loadings carry the compensating scale: only the product
#' \eqn{\alpha \beta^{\prime}} is identified, so their prior variance is shrunk
#' by \eqn{1 - \rho^2}, leaving the product on the scale \code{coef$v_i} asks
#' for.
#' 
#' When \eqn{\rho} is drawn, those last two are computed from \code{coint$rho}
#' once and do not follow the draw: the state before the sample keeps the normal
#' prior built here, and the loadings keep the shrinkage built here. This is a
#' deliberate departure from Koop et al. (2011), in whom the state before the
#' sample is the stationary distribution of whatever \eqn{\rho} currently is.
#' It is also what makes the draw an exact Gibbs block rather than their
#' Metropolis-within-Gibbs step: with the state before the sample free of
#' \eqn{\rho}, the conditional posterior of \eqn{\rho} is a normal truncated to
#' the prior's support. Draws of \eqn{\rho} are returned in
#' \code{object$posterior$beta$rho}.
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
#'   freedom of the inverse Wishart prior. The rank \eqn{r} of the cointegration matrix is added
#'   to the value given.}
#'   \item{\code{scale}}{a positive numeric specifying the prior error variance of the endogenous
#'   variables in the inverse Wishart prior.}
#'   \item{\code{shape}}{for \code{"gamma"} and \code{"gamma+covar"} a non-negative numeric, or a
#'   character expression in \code{k} as for \code{df}, specifying the prior shape parameter of the
#'   error variances, to which the rank \eqn{r} is added as well. For models with stochastic
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

    # rho is drawn only if the support of its prior is given, and then both ends
    # of it are needed: one alone would leave the sampler to invent the other,
    # and which end is missing changes the model rather than a detail of it.
    has_rho_min <- "rho_min" %in% names(coint)
    has_rho_max <- "rho_max" %in% names(coint)
    if (xor(has_rho_min, has_rho_max)) {
      stop("Arguments 'coint$rho_min' and 'coint$rho_max' must be specified together. ",
           "Leave both out to hold rho fixed at 'coint$rho'.")
    }
    if (has_rho_min) {
      if (coint[["rho_min"]] <= 0 | coint[["rho_max"]] > 1 |
          coint[["rho_min"]] >= coint[["rho_max"]]) {
        stop("Arguments 'coint$rho_min' and 'coint$rho_max' must satisfy ",
             "0 < rho_min < rho_max <= 1.")
      }
      if (coint[["rho"]] < coint[["rho_min"]] | coint[["rho"]] > coint[["rho_max"]]) {
        stop("Argument 'coint$rho' must lie between 'coint$rho_min' and 'coint$rho_max': ",
             "when rho is drawn it is the value the chain starts at.")
      }
    }

    # A time varying space has no matric-variate prior for a numeric p_tau_i to
    # describe; the only thing p_tau_i can ask for here is the transition
    # centred on the ML estimate.
    coint_p_tau_ml <- identical(coint[["p_tau_i"]], "ml")
    if (!is.null(coint[["p_tau_i"]]) && !coint_p_tau_ml) {
      warning("Argument 'coint$p_tau_i' is only used with the value \"ml\" for VEC models ",
              "with time varying cointegration parameters and is ignored.")
    }
  } else {
    if (!"v_i" %in% names(coint)) {
      stop("Argument 'coint$v_i' must be specified for VEC models with constant cointegration parameters.")
    }
    if (!"p_tau_i" %in% names(coint)) {
      stop("Argument 'coint$p_tau_i' must be specified for VEC models with constant cointegration parameters.")
    }
    coint_v_ml <- identical(coint[["v_i"]], "ml")
    coint_p_tau_ml <- identical(coint[["p_tau_i"]], "ml")
    if (!coint_v_ml && !(is.numeric(coint[["v_i"]]) && length(coint[["v_i"]]) == 1 &&
                         isTRUE(coint[["v_i"]] >= 0))) {
      stop("Argument 'coint$v_i' must be a non-negative number or \"ml\".")
    }
    if (!coint_p_tau_ml && !is.numeric(coint[["p_tau_i"]])) {
      stop("Argument 'coint$p_tau_i' must be numeric or \"ml\".")
    }
    # The sampler only ever uses the product v_i * p_tau_i, so under zero
    # shrinkage any p_tau_i is the uniform prior on the space and a centre
    # estimated for it would be silently thrown away.
    if (coint_p_tau_ml && !coint_v_ml && coint[["v_i"]] == 0) {
      stop("Argument 'coint$v_i' must be positive for 'coint$p_tau_i = \"ml\"': ",
           "with zero shrinkage the prior on the cointegration space is uniform ",
           "whatever 'coint$p_tau_i' is.")
    }
  }
  if ("weight" %in% names(coint)) {
    if (!(is.numeric(coint[["weight"]]) && length(coint[["weight"]]) == 1 &&
          isTRUE(coint[["weight"]] > 0))) {
      stop("Argument 'coint$weight' must be a positive number.")
    }
    if (!coint_p_tau_ml) {
      warning("Argument 'coint$weight' is only used with 'coint$p_tau_i = \"ml\"' and is ignored.")
    }
  }
  
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
  if (r > 0) {
    
    n_ect <- k * (k + m)
    if (object[["model"]][["n_restricted"]] > 0) {
      n_ect <- n_ect + object[["model"]][["n_restricted"]] * k
    }
    
    n_alpha <- r * k
    n_beta <- r * n_ect / k
    
    if (object[["model"]][["tvp"]]) {

      # The cointegration vectors are a state path of their own, with
      # beta_t = rho beta_{t-1} + eta_t and eta_t ~ N(0, I). rho belongs with the
      # prior either way: fixed, it is the value the sampler keeps, and drawn, it
      # is the value the chain starts at.
      #
      # The prior on the state before the sample is that path's own stationary
      # distribution, N(0, I / (1 - rho^2)), whose precision is (1 - rho^2) I.
      # That is what makes it proper: at rho = 1 the state equation is a random
      # walk whose variance grows without bound, and beta, being identified only
      # up to scale, has nothing to pull it back.
      #
      # The shrinkage of the constant model's cointegration space, coint$v_i, has
      # no counterpart here: the space is not drawn from a matric-variate prior
      # but followed period by period. Its central location does, as the
      # transition below when coint$p_tau_i = "ml".
      object[["priors"]][["beta"]] <- list("type" = "cointspace",
                                           "rho" = coint[["rho"]],
                                           "mu" = matrix(0, n_beta),
                                           "v_inv" = diag(1 - coint[["rho"]]^2, n_beta))

      # The support of the uniform prior on rho, which is what turns it from a
      # fixed hyperparameter into a drawn one.
      #
      # Note what stays fixed when it is drawn: v_inv just above, and the
      # shrinkage of the loadings further down, are both computed from
      # coint$rho once and do not follow the draw. Under Koop et al. (2011) the
      # state before the sample is the stationary distribution of whatever rho
      # currently is, so those two would move with it; here they are an ordinary
      # normal prior, pinned at the value rho starts from. That is a deliberate
      # difference and the reason the draw is a plain Gibbs block rather than
      # their Metropolis-within-Gibbs step -- see the vendored
      # draw_coint_rho() in src/core/models/vec_support.h.
      if (has_rho_min) {
        object[["priors"]][["beta"]][["rho_min"]] <- coint[["rho_min"]]
        object[["priors"]][["beta"]][["rho_max"]] <- coint[["rho_max"]]
      }

      # The informative marginal prior of Koop et al. (2011, working paper
      # version, eq. 12): the transition becomes rho (I_r kron P_tau) with
      # P_tau = H H' + H_perp T H_perp', H an orthonormal basis of Johansen's
      # estimate of the space. The part of beta along sp(H) keeps rho, the part
      # off it decays at rho T, so the space at every t has its mode at sp(H).
      #
      # T is set from the spread the prior gives the tilt delta of beta_t away
      # from sp(H). The transition is stationary with H_perp' beta_t having
      # variance (I - rho^2 T^2)^-1 against 1 / (1 - rho^2) along H, so for a
      # beta_t of typical length the tilt has spread
      # T* = (1 - rho^2) (I - rho^2 T^2)^-1. That is matched to Johansen's
      # sampling variance of the tilt, (H_perp' S11 H_perp)^-1 divided by the
      # smallest eigenvalue of alpha_h' Omega^-1 alpha_h -- exact for rank one,
      # the more cautious direction otherwise -- with its precision multiplied
      # by 'weight'. Eigenvalue by eigenvalue, tau^2 = (1 - (1 - rho^2) / t*) / rho^2.
      #
      # t* is bounded on both sides. Above one tau would exceed one and the
      # prior would push away from sp(H); it is capped there, which is the
      # identity. Below 1 - rho^2 there is no tau: even T = 0 leaves the tilt
      # that much room per period, so rho limits how informative the prior can
      # be, and the request is floored at T = 0 with a warning.
      #
      # The floor is the tightest *prior*, not the tightest posterior. T is also
      # the persistence of the tilt, so near zero each period's tilt is informed
      # by that period's data alone and the posterior can widen again -- on e6 at
      # rho = 0.999, T = 0.44 left a tighter posterior than T = 0. The roxygen
      # details say so for the user.
      #
      # The state before the sample gets the stationary distribution the
      # transition implies, N(0, I_r kron P_tau* / (1 - rho^2)) with
      # P_tau* = H H' + H_perp T* H_perp'. Both are computed at coint$rho and, if
      # rho is drawn, do not follow the draw -- as with v_inv above.
      if (coint_p_tau_ml) {
        k_beta <- n_ect / k
        if (NROW(object[["data"]][["train"]][["y"]]) <=
            NCOL(object[["data"]][["train"]][["x"]]) + k_beta) {
          stop("Not enough observations for the maximum likelihood estimate that ",
               "'coint$p_tau_i = \"ml\"' is based on.")
        }

        rho <- coint[["rho"]]
        weight <- if (is.null(coint[["weight"]])) 1 else coint[["weight"]]

        ml <- .coint_ml(object)
        h <- ml[["beta"]] %*% solve(.mroot(crossprod(ml[["beta"]])))
        h_perp <- qr.Q(qr(h), complete = TRUE)[, -(1:r), drop = FALSE]
        alpha_h <- ml[["alpha"]] %*% t(crossprod(h, ml[["beta"]]))

        info_alpha <- min(eigen(t(alpha_h) %*% solve(ml[["omega"]]) %*% alpha_h,
                                symmetric = TRUE)$values)
        target <- solve(t(h_perp) %*% tcrossprod(ml[["r1"]]) %*% h_perp) /
          (weight * info_alpha)
        target <- eigen((target + t(target)) / 2, symmetric = TRUE)

        floor_t <- 1 - rho^2
        if (any(target[["values"]] < floor_t)) {
          warning("Argument 'coint$weight' asks for a tighter prior on the cointegration ",
                  "space than 'coint$rho' = ", rho, " allows; the transition is set to ",
                  "T = 0 in the directions concerned. A value of 'coint$rho' closer to one ",
                  "leaves room for more.")
        }
        t_star <- pmin(pmax(target[["values"]], floor_t), 1)

        # Capped at one in every direction is no prior on the direction at all,
        # and the noninformative prior built above stays as it is.
        if (any(t_star < 1)) {
          tau <- pmin(sqrt(pmax(0, 1 - floor_t / t_star)) / rho, 1)
          basis <- h_perp %*% target[["vectors"]]
          p_tau <- tcrossprod(h) + basis %*% diag(tau, nrow = length(tau)) %*% t(basis)
          p_tau_star_inv <- tcrossprod(h) +
            basis %*% diag(1 / t_star, nrow = length(t_star)) %*% t(basis)

          object[["priors"]][["beta"]][["p_tau"]] <- (p_tau + t(p_tau)) / 2
          object[["priors"]][["beta"]][["v_inv"]] <-
            kronecker(diag(1, r), floor_t * (p_tau_star_inv + t(p_tau_star_inv)) / 2)
        }
      }
    } else {

      k_beta <- n_ect / k
      coint_v_inv <- coint[["v_i"]]

      if (coint_v_ml | coint_p_tau_ml) {
        if (NROW(object[["data"]][["train"]][["y"]]) <=
            NCOL(object[["data"]][["train"]][["x"]]) + k_beta) {
          stop("Not enough observations for the maximum likelihood estimate that ",
               "'coint$v_i = \"ml\"' or 'coint$p_tau_i = \"ml\"' is based on.")
        }

        # Johansen's estimate, re-expressed for the orthonormal basis H of the
        # space it spans: alpha beta' = alpha_h H'.
        ml <- .coint_ml(object)
        h <- ml[["beta"]] %*% solve(.mroot(crossprod(ml[["beta"]])))
        h_perp <- qr.Q(qr(h), complete = TRUE)[, -(1:r), drop = FALSE]
        alpha_h <- ml[["alpha"]] %*% t(crossprod(h, ml[["beta"]]))
      }

      # With beta close to H the loadings have prior N(0, Sigma / v), so the
      # prior mean of tr(alpha' Sigma^-1 alpha) is r k / v. It is set to 100 times
      # its ML estimate, which keeps the loadings weakly shrunk on the scale the
      # data put them on.
      if (coint_v_ml) {
        coint_v_inv <- r * k /
          (100 * sum(diag(t(alpha_h) %*% solve(ml[["omega"]]) %*% alpha_h)))
      }

      if (coint_p_tau_ml) {

        # Write beta = H + H_perp delta. Given alpha the prior of Koop et al.
        # (2010) implies vec(delta) ~ N(0, (alpha' Sigma^-1 alpha)^-1 kron T / v)
        # for P_tau^-1 = H H' + H_perp T^-1 H_perp'. Johansen's estimator has the
        # same Kronecker form, vec(delta) ~ N(0, (alpha' Omega^-1 alpha)^-1 kron
        # (H_perp' S11 H_perp)^-1) with S11 the sum of the residual cross
        # products, so T = v / weight * (H_perp' S11 H_perp)^-1 makes the prior on
        # the space that sampling distribution, worth 'weight' samples.
        #
        # Eigenvalues of T above one are capped at one. One is the uniform prior
        # on the space; beyond it the prior would favour the complement of the
        # estimate rather than being weakly informative about it.
        weight <- if (is.null(coint[["weight"]])) 1 else coint[["weight"]]
        tilt <- solve(t(h_perp) %*% tcrossprod(ml[["r1"]]) %*% h_perp)
        tilt <- eigen(coint_v_inv / weight * (tilt + t(tilt)) / 2, symmetric = TRUE)
        tau <- pmin(tilt[["values"]], 1)
        t_inv <- tilt[["vectors"]] %*% diag(1 / tau, nrow = length(tau)) %*% t(tilt[["vectors"]])
        p_tau_inv <- tcrossprod(h) + h_perp %*% t_inv %*% t(h_perp)
        p_tau_inv <- (p_tau_inv + t(p_tau_inv)) / 2
      } else if (is.matrix(coint[["p_tau_i"]])) {
        p_tau_inv <- coint[["p_tau_i"]]
        if (!all(dim(p_tau_inv) == k_beta)) {
          stop("Argument 'coint$p_tau_i' must be a ", k_beta, " x ", k_beta,
               " matrix, one row and column per series in the error correction term.")
        }
        if (!isSymmetric(unname(p_tau_inv))) {
          stop("Argument 'coint$p_tau_i' must be a symmetric matrix.")
        }
      } else {
        if (!length(coint[["p_tau_i"]]) %in% c(1, k_beta)) {
          stop("Argument 'coint$p_tau_i' must have one element or one per series ",
               "in the error correction term.")
        }
        p_tau_inv <- diag(coint[["p_tau_i"]], k_beta)
      }

      object[["priors"]][["beta"]] <- list("type" = "cointspace",
                                           "v_inv" = coint_v_inv,
                                           "p_tau_inv" = p_tau_inv)
    }
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
    
    # Add rank to degrees of freedom for cointegration model
    if (!is.na(object[["model"]][["rank"]])) {
      help_df <- help_df + r
    }
    
    if (error_prior == "wishart") {
      object[["priors"]][[pos]][["df"]] <- help_df
    }
    if (error_prior == "gamma") {
      object[["priors"]][[pos]][["shape"]] <- matrix(help_df, k)
    }
  }
  
  return(object)
}
