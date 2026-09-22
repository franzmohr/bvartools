#' Build the Prior on the Cointegration Space
#'
#' Checks the specification of the prior on the cointegration space of a VEC
#' model and builds it, in the form in which \code{\link{add_priors}} stores it
#' as element \code{priors$beta}.
#'
#' @param object a list of class 'bvecmodel', or any model with the same layout:
#' elements \code{model$rank} and \code{model$tvp}, and the regressors
#' \code{data$train$y}, \code{data$train$w} and \code{data$train$x}.
#' @param coint a named list of prior specifications for the coefficients
#' determining the cointegration space. It has no default. See section
#' 'Prior on the cointegration space'.
#'
#' @details
#' \code{\link{add_priors}} calls this function for VEC models, so it rarely
#' needs to be called directly. It is exported for packages that build VEC
#' models of their own layout, such as the sub-models of a global VEC model,
#' whose error correction term also holds weakly exogenous and global
#' variables. Calling it gives them the same checks and the same prior.
#'
#' All dimensions are taken from the regressors, not from the model
#' specification: the prior has one row and column for each series in the error
#' correction term \code{data$train$w}, and one block for each of the
#' \code{model$rank} cointegration vectors.
#'
#' @section Prior on the cointegration space:
#' Argument \code{coint} can contain the following elements. Which of them are
#' required depends on whether the cointegration vectors are constant or time
#' varying, since the two are given different kinds of prior: a matric-variate
#' prior on the cointegration space in the first case, and a state equation in
#' the second. Any other element raises an error.
#' \describe{
#'   \item{\code{v_i}}{non-negative numeric specifying the shrinkage of the cointegration space prior,
#'   or \code{"ml"}. See below. Required for models with constant cointegration parameters and not
#'   used otherwise.}
#'   \item{\code{p_tau_i}}{the inverse of the matrix \eqn{P_\tau}, which determines the
#'   central location of the cointegration space \eqn{sp(\beta)}. Either a numeric of
#'   its diagonal elements, a full symmetric matrix, or \code{"ml"}, with one row and
#'   column per series in the error correction term \code{data$train$w}. See below.
#'   Required for models with constant cointegration parameters. For models with time
#'   varying cointegration parameters only \code{"ml"} is used.}
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
#'   \item{\code{g_i}}{the inverse of the matrix \eqn{G} that scales the prior of the
#'   loadings, for models with constant cointegration parameters and stochastic volatility,
#'   \code{error = "sv"} or \code{"sv+covar"}, and refused for every other model. Either a
#'   numeric of its diagonal elements, a full symmetric positive definite matrix with one row
#'   and column per endogenous variable, or \code{"ml"} for the inverse of the maximum
#'   likelihood estimate of the error covariance. Optional. See below.}
#' }
#' For a model with constant cointegration parameters the prior is that of
#' Koop et al. (2010). The sampler uses \code{v_i} and \code{p_tau_i} only through
#' their product, so with \code{v_i = 0} the prior on the cointegration space is
#' uniform whatever \code{p_tau_i} is. An informative prior on the space
#' therefore needs a positive \code{v_i}, which also shrinks the loadings: for
#' \eqn{\beta} close to the centre of the space they have prior
#' \eqn{N(0, \Sigma / v)}.
#'
#' In Koop et al. (2010) the loadings' prior is scaled by a matrix \eqn{G}, which
#' may be the error covariance \eqn{\Sigma} or any fixed, known matrix. The models
#' with a constant error covariance take \eqn{G = \Sigma}. With stochastic
#' volatility the covariance differs from period to period, so \eqn{G} is fixed
#' for the whole run instead: \eqn{G^{-1}} is \code{g_i} if it is given, and
#' otherwise the error precision implied by the starting values of the
#' log-volatilities, averaged over the sample once before the first draw. That
#' fallback depends on \code{\link{add_initial_values}}, so giving \code{g_i}
#' makes the prior independent of how the chain is started. \code{g_i = "ml"}
#' uses Johansen's (1995) estimate of the error covariance, the same one
#' \code{v_i = "ml"} is based on.
#'
#' With \code{p_tau_i = "ml"} the prior is centred on the space spanned by
#' Johansen's (1995) maximum likelihood estimate \eqn{\hat{\beta}}, computed from
#' the error correction term as it is stored in the model -- so call
#' \code{\link{scale_error_correction}} first if the series should be scaled;
#' scaling is refused afterwards.
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
#' For a model with time varying cointegration parameters the state equation is
#' \eqn{\beta_t = \rho \beta_{t-1} + \eta_t} with \eqn{\eta_t \sim N(0, I)}, and
#' the prior on the state before the sample is that equation's own stationary
#' distribution, \eqn{N(0, I / (1 - \rho^2))}. This is what makes the prior
#' proper, so a \eqn{\rho} close to one is intended and one further from it draws
#' a warning. \code{\link{add_priors}} gives the loadings the compensating scale:
#' only the product \eqn{\alpha \beta^{\prime}} is identified, so their prior
#' variance is shrunk by \eqn{1 - \rho^2}, leaving the product on the scale
#' \code{coef$v_i} asks for.
#'
#' When \eqn{\rho} is drawn, those last two are computed from \code{coint$rho}
#' once and do not follow the draw: the state before the sample keeps the normal
#' prior built here, and the loadings keep the shrinkage of
#' \code{\link{add_priors}}. This is a deliberate departure from Koop et al.
#' (2011), in whom the state before the sample is the stationary distribution of
#' whatever \eqn{\rho} currently is. It is also what makes the draw an exact
#' Gibbs block rather than their Metropolis-within-Gibbs step: with the state
#' before the sample free of \eqn{\rho}, the conditional posterior of \eqn{\rho}
#' is a normal truncated to the prior's support. Draws of \eqn{\rho} are
#' returned in \code{object$posterior$beta$rho}.
#'
#' How far the space moves is not set by the prior on any variance. The steps
#' \eqn{\eta_t} have unit variance whatever \code{coef$rate} and
#' \code{coef$rate_alpha} are, and how far a step moves
#' \eqn{\Pi_t = \alpha_t \beta_t^{\prime}} depends on the scale of \eqn{\beta_t}.
#' Only the product is identified, so the posterior can pair small loadings with
#' large cointegration vectors, which a unit step barely turns, or large loadings
#' with small ones, which it turns a lot. Two things decide which. One is the
#' prior precision of the loadings, \code{coef$v_i} \eqn{/ (1 - \rho^2)}, which
#' restrains the drift only if it is tight relative to the scale of \eqn{\Pi} the
#' data call for, and that scale depends on the units of the data. The other is
#' the series in the error correction term: a step moves
#' \eqn{\beta_t^{\prime} w_t} by \eqn{\eta_t^{\prime} w_t}, which grows with the
#' levels of the series and not only with their variation. For log levels far
#' from zero the shift is of the order of the levels themselves, and the sampler
#' either lets it act as a random walk intercept, which no prior on the
#' deterministic terms controls and which can absorb the residuals of an
#' equation, or draws loadings close to zero, which switches the error
#' correction term off.
#'
#' Both have been seen. In the US sub-model of the \pkg{bgvars} data set
#' \code{gvar2023}, with the oil price among its endogenous variables and a rank
#' of one, a prior that pinned the state variances of all coefficients, the
#' loadings and the constant included, at 1e-14 still left the residual standard
#' deviation of one equation as low as 0.1 to 0.5 of the maximum likelihood one
#' in some chains, with which equation and how low depending on the chain, and a
#' saved chain reported a posterior mean log-likelihood several hundred above the
#' maximum of the constant coefficient model. The chains do not move between
#' these outcomes within a few thousand draws, so a single chain can report any
#' of them. The Austrian sub-model under the same prior lost fit that the
#' constant coefficient model has, and the chains that were inspected drew
#' loadings close to zero. On data simulated from a constant coefficient VEC
#' model, whose posterior under such a prior should reproduce the maximum
#' likelihood fit, single chains did so only with the series centred and
#' \code{coef$v_i = 10}, or around zero and \code{coef$v_i = 100}. Around 12,
#' where log levels usually are, a chain drew loadings close to zero, and around
#' zero with \code{coef$v_i = 1} one fitted part of the residuals.
#'
#' The same mechanism makes \code{\link{scale_error_correction}} risky in a model
#' whose coefficients vary over time. It divides the series by the standard
#' deviation of their differences, which makes them, and with them both a drift
#' in a loading and a step of \eqn{\beta_t^{\prime} w_t}, larger by the same
#' factor. The coefficient paths can then absorb the residuals of an equation
#' almost entirely: its residual variance -- or, under stochastic volatility, its
#' volatility path -- is reported far below that of a least squares fit of the
#' same regressors and nearly constant over the sample. In the model of the
#' vignette on time varying parameters and stochastic volatility in error
#' correction models, \code{vignette("tvp-sv-vec", package = "bvartools")}, a
#' scaled error correction term did that to one equation in each of four chains
#' with different seeds, and the unscaled one with a small \code{coef$rate} in
#' none.
#'
#' It is safer to leave the series unscaled, to set \code{coef$v_i} with the
#' scale of \eqn{\Pi} in mind as well as choosing small rates, and to compare the
#' residual variances with those of a least squares fit of the same regressors,
#' in chains with different seeds, for a model that is to be used. Small rates
#' alone do not guarantee that a model passes that check, and passing it is what
#' says that neither the coefficients nor the cointegration vectors fit the
#' residuals. \code{scale_error_correction(object, scale = FALSE, centre = TRUE)}
#' centres the series before posterior simulation, which removes the part of a
#' step of \eqn{\beta_t^{\prime} w_t} that comes from the levels of the series
#' rather than from their variation, and \code{\link{rescale_error_correction}}
#' writes the draws back in terms of the series as they are, so that
#' \code{\link{vec_to_var}} and the forecasts can use them. It does not reach the
#' other channel, the scale of the loadings, and it is not a full remedy: in the
#' US sub-model above, centring alone roughly halved how far the error correction
#' term and the constant moved over the sample and still left the oil price
#' equation at about 0.8 of the maximum likelihood residual standard deviation in
#' two chains. Whether a centred model is right is again what the comparison with
#' least squares shows.
#'
#' For a model with time varying cointegration parameters
#' \code{p_tau_i = "ml"} centres the marginal prior of the
#' cointegration space on Johansen's estimate as well, using the informative
#' marginal prior of Koop et al. (2011, working paper version). The state
#' equation becomes
#' \eqn{\beta_t = \rho (I_r \otimes P_\tau) \beta_{t-1} + \eta_t} with
#' \eqn{P_\tau = H H^{\prime} + H_\perp T H_\perp^{\prime}}: the part of
#' \eqn{\beta_t} along \eqn{sp(H)} keeps \eqn{\rho}, the part off it decays
#' faster, and the mode of the marginal distribution of \eqn{sp(\beta_t)} is
#' \eqn{sp(H)} in every period. \eqn{T} is chosen so that the prior spread of the
#' tilt of \eqn{\beta_t} away from \eqn{sp(H)} in a period, approximately
#' \eqn{T^* = (1 - \rho^2)(I - \rho^2 T^2)^{-1}}, is the sampling variance of
#' Johansen's estimator with its precision multiplied by \code{weight}. The
#' state before the sample is given the stationary distribution the transition
#' implies. The transition is stored as element \code{p_tau}.
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
#' @return \code{NULL} for a model without cointegration, \code{model$rank = 0},
#' whose \code{coint} is checked all the same. Otherwise a list with
#' \code{type = "cointspace"} and, for constant cointegration parameters,
#' \code{v_inv} and \code{p_tau_inv}, and \code{g_inv} if \code{g_i} was given,
#' or, for time varying ones, \code{rho},
#' \code{mu} and \code{v_inv} of the state equation, together with
#' \code{rho_min} and \code{rho_max} for a uniform prior on \eqn{\rho} and the
#' transition \code{p_tau} added by \code{p_tau_i = "ml"}.
#'
#' @references
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
#' @examples
#'
#' # Load data
#' data("e6")
#' e6 <- e6 * 100
#'
#' # Generate model
#' model <- create_bvecmodel(e6, p = 2, r = 1,
#'                           const = "unrestricted",
#'                           iterations = 10, burnin = 10)
#'
#' # A prior centred on the maximum likelihood estimate of the space
#' prior <- cointspace_prior(model, coint = list(v_i = 0.01, p_tau_i = "ml"))
#' prior$p_tau_inv
#'
#' @export
cointspace_prior <- function(object, coint) {

  # Checks ----

  allowed_coint_arguments <- c("v_i", "p_tau_i", "weight", "rho", "rho_min", "rho_max", "g_i")
  for (i in names(coint)) {
    if (!i %in% allowed_coint_arguments) {
      stop(paste0("Element '", i, "' in argument 'coint' is not recognised."))
    }
  }

  # G is the error covariance in every other model, so only the constant VEC
  # with stochastic volatility has a G to give. The core refuses one anywhere
  # else as well; refusing it here says so in the terms of this function.
  coint_g_ml <- identical(coint[["g_i"]], "ml")
  if (!is.null(coint[["g_i"]])) {
    if (isTRUE(object[["model"]][["tvp"]]) ||
        !isTRUE(object[["model"]][["error"]] %in% c("sv", "sv+covar"))) {
      stop("Argument 'coint$g_i' is only used for VEC models with constant cointegration ",
           "parameters and stochastic volatility, error = \"sv\" or \"sv+covar\". Every ",
           "other model scales the prior of the loadings by its error covariance.")
    }
    if (!coint_g_ml && !is.numeric(coint[["g_i"]])) {
      stop("Argument 'coint$g_i' must be numeric or \"ml\".")
    }
  }

  if (object[["model"]][["tvp"]]) {
    if (!"rho" %in% names(coint)) {
      stop("Argument 'coint$rho' must be specified for VEC models with time varying cointegration parameters.")
    }
    # The stationary distribution of the state, which is the prior of the state
    # before the sample, has variance 1 / (1 - rho^2): below -1 that is
    # negative, and the sampler was handed a prior that is not one.
    if (!is.numeric(coint[["rho"]]) || length(coint[["rho"]]) != 1 || is.na(coint[["rho"]]) ||
        coint[["rho"]] >= 1 || coint[["rho"]] <= -1) {
      stop("Argument 'coint$rho' must be a number larger than -1 and smaller than 1.")
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

  # Prior ----

  r <- object[["model"]][["rank"]]
  if (!(r > 0)) {
    return(NULL)
  }

  # The dimensions come from the regressors rather than from the model
  # specification. For a model of create_bvecmodel() the error correction term
  # has k + m + n_restricted columns and the two agree; a model of another layout
  # -- a sub-model of a global VEC model, whose error correction term also holds
  # weakly exogenous and global variables, with 'm' counting their lags -- would
  # otherwise get a prior of the wrong size.
  if (is.null(object[["data"]][["train"]][["w"]])) {
    stop("Argument 'object' has a cointegration rank of ", r, " but no error ",
         "correction term 'data$train$w'.")
  }
  k <- NCOL(object[["data"]][["train"]][["y"]])
  k_beta <- NCOL(object[["data"]][["train"]][["w"]])
  n_beta <- r * k_beta

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
    prior <- list("type" = "cointspace",
                  "rho" = coint[["rho"]],
                  "mu" = matrix(0, n_beta),
                  "v_inv" = diag(1 - coint[["rho"]]^2, n_beta))

    # The support of the uniform prior on rho, which is what turns it from a
    # fixed hyperparameter into a drawn one.
    #
    # Note what stays fixed when it is drawn: v_inv just above, and the
    # shrinkage of the loadings in add_priors(), are both computed from
    # coint$rho once and do not follow the draw. Under Koop et al. (2011) the
    # state before the sample is the stationary distribution of whatever rho
    # currently is, so those two would move with it; here they are an ordinary
    # normal prior, pinned at the value rho starts from. That is a deliberate
    # difference and the reason the draw is a plain Gibbs block rather than
    # their Metropolis-within-Gibbs step -- see the vendored
    # draw_coint_rho() in src/core/models/vec_support.h.
    if (has_rho_min) {
      prior[["rho_min"]] <- coint[["rho_min"]]
      prior[["rho_max"]] <- coint[["rho_max"]]
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
    # section says so for the user.
    #
    # The state before the sample gets the stationary distribution the
    # transition implies, N(0, I_r kron P_tau* / (1 - rho^2)) with
    # P_tau* = H H' + H_perp T* H_perp'. Both are computed at coint$rho and, if
    # rho is drawn, do not follow the draw -- as with v_inv above.
    if (coint_p_tau_ml) {
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

        prior[["p_tau"]] <- (p_tau + t(p_tau)) / 2
        prior[["v_inv"]] <-
          kronecker(diag(1, r), floor_t * (p_tau_star_inv + t(p_tau_star_inv)) / 2)
      }
    }
  } else {

    coint_v_inv <- coint[["v_i"]]

    if (coint_v_ml | coint_p_tau_ml | coint_g_ml) {
      if (NROW(object[["data"]][["train"]][["y"]]) <=
          NCOL(object[["data"]][["train"]][["x"]]) + k_beta) {
        stop("Not enough observations for the maximum likelihood estimate that ",
             "'coint$v_i = \"ml\"', 'coint$p_tau_i = \"ml\"' or 'coint$g_i = \"ml\"' ",
             "is based on.")
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

    prior <- list("type" = "cointspace",
                  "v_inv" = coint_v_inv,
                  "p_tau_inv" = p_tau_inv)

    # G^-1, which the stochastic volatility VEC scales the loadings' prior by
    # for the whole run. Left out, the core takes the precision the starting
    # log-volatilities imply, averaged over the sample.
    if (!is.null(coint[["g_i"]])) {
      if (coint_g_ml) {
        g_inv <- solve(ml[["omega"]])
      } else if (is.matrix(coint[["g_i"]])) {
        g_inv <- coint[["g_i"]]
        if (!all(dim(g_inv) == k)) {
          stop("Argument 'coint$g_i' must be a ", k, " x ", k,
               " matrix, one row and column per endogenous variable.")
        }
        if (!isSymmetric(unname(g_inv))) {
          stop("Argument 'coint$g_i' must be a symmetric matrix.")
        }
      } else {
        if (!length(coint[["g_i"]]) %in% c(1, k)) {
          stop("Argument 'coint$g_i' must have one element or one per endogenous variable.")
        }
        g_inv <- diag(coint[["g_i"]], k)
      }
      g_inv <- (g_inv + t(g_inv)) / 2
      if (any(!is.finite(g_inv)) ||
          min(eigen(g_inv, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
        stop("Argument 'coint$g_i' must be positive definite.")
      }
      prior[["g_inv"]] <- g_inv
    }
  }

  return(prior)
}
