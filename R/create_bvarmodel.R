#' Create a Vector Autoregressive Model
#' 
#' Produces the input for the estimation of a vector autoregressive (VAR) model.
#' 
#' @param data a time-series object of endogenous variables.
#' @param p an integer vector of the lag order (default is \code{p = 2}).
#' @param exogen an optional time-series object of external regressors.
#' @param s an integer vector of the lag order of the external regressors (default is \code{s = 2}).
#' Ignored if \code{exogen} is \code{NULL}.
#' @param deterministic a character specifying which deterministic terms should
#' be included. Available values are \code{"none"}, \code{"const"} (default) for an intercept,
#' \code{"trend"} for a linear trend, and \code{"both"} for an intercept with a linear trend.
#' @param seasonal logical. If \code{TRUE}, seasonal dummy variables are
#' generated as additional deterministic terms. The amount of dummies depends on the frequency of the
#' time-series object provided in \code{data}. Defaults to \code{FALSE}.
#' @param iid an optional character vector naming endogenous variables whose
#' equations carry no coefficients at all -- no lags, no deterministic terms,
#' nothing. They are white noise, and reach the rest of the model only through
#' the error covariance. They must be the first columns of \code{data}. See
#' 'Details'.
#' @param structural logical indicating whether data should be prepared for the estimation of a
#' structural VAR model. Defaults to \code{FALSE}.
#' @param tvp logical indicating whether the model parameters are time varying.
#' Defaults to \code{FALSE}.
#' @param error character specifying the model that should be used for the estimation
#' of the covariance matrix of the error term. Default is \code{"wishart"}. See 'Details'.
#' @param quantile a numeric vector of quantiles in the interval \eqn{(0, 1)} that should be
#' estimated. Only used, if \code{error = "ald"}. Defaults to \code{0.5}, the median. One model
#' is created per quantile, so a vector produces a list of models, unless \code{quantile_grid}
#' is \code{TRUE}. See 'Details'.
#' @param varsel character specifying the type of variable selection algorithm
#' that should be employed. Default is \code{"none"}. See 'Details'.
#' @param algorithm algorithm that should be used for posterior simulation. If
#' \code{NULL} (default), the algorithm is named by \code{tvp} and \code{error}.
#' The one non-standard option is \code{"discount"}. See 'Details'.
#' @param delta_beta,delta_sigma numeric discount factors in \eqn{(0, 1]} of the
#' discounted model, the first governing the coefficients and the second the
#' error covariance. Both default to 1, at which the quantity they govern does
#' not move. Ignored unless \code{algorithm = "discount"}, and a vector in
#' either produces one model per value. See 'Details'.
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 20000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 2000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' @param thin an integer thinning interval of the sampler (defaults to 1). After
#' the burn-in the sampler keeps the last of every \code{thin} draws, so it runs
#' \code{burnin + iterations * thin} draws and still keeps \code{iterations}. Unlike
#' \code{\link[=thin.bvarmodel]{thin}}, which thins draws already made, the draws
#' that are not kept are never held in memory.
#' @param missing a character, what to do with periods in which a series of
#' \code{data} is \code{NA}. \code{"omit"} (default) drops them, as
#' \code{\link[stats]{na.omit}} does. \code{"estimate"} keeps them and estimates
#' what was not observed together with the model. See 'Details'.
#' @param aggregate an optional named list, one element per series of
#' \code{data} that is observed at a lower frequency than the model, holding the
#' weights with which it aggregates the periods of the model, oldest first, as
#' \code{\link{aggregation_weights}} returns them. The series is \code{NA} in
#' every period in which it is not observed, and its observations stand in the
#' last period they aggregate. Implies \code{missing = "estimate"}.
#' @param soft an optional character vector of series of \code{data} whose
#' observations hold up to a normal measurement error instead of exactly, each
#' series with an error variance of its own. See 'Details'.
#' @param quantile_grid logical. If \code{TRUE}, the quantiles in \code{quantile} are
#' estimated as one structural quantile VAR rather than as one model each. Needs
#' \code{error = "ald"} and constant coefficients. Defaults to \code{FALSE}. See 'Details'.
#'
#' @details The function produces the data matrices for vector autoregressive (VAR)
#' models, which can also include unmodelled, non-deterministic variables:
#' \deqn{A_0 y_t = \sum_{i=1}^{p} A_i y_{t - i} +
#' \sum_{i=0}^{s} B_i x_{t - i} +
#' C d_t + u_t,}
#' where
#' \eqn{y_t} is a K-dimensional vector of endogenous variables,
#' \eqn{A_0} is a \eqn{K \times K} coefficient matrix of contemporaneous endogenous variables,
#' \eqn{A_i} is a \eqn{K \times K} coefficient matrix of endogenous variables,
#' \eqn{x_t} is an M-dimensional vector of exogenous regressors and
#' \eqn{B_i} its corresponding \eqn{K \times M} coefficient matrix.
#' \eqn{d_t} is an N-dimensional vector of deterministic terms and
#' \eqn{C} its corresponding \eqn{K \times N} coefficient matrix.
#' \eqn{p} is the lag order of endogenous variables, \eqn{s} is the lag
#' order of exogenous variables, and \eqn{u_t} is an error term.
#' 
#' The model can be rewritten as
#' \deqn{A_0 y_t = Z_t a + u_t,}
#' where \eqn{Z_t} is a \eqn{KT \times K * (Kp + M(s + 1) + N)} data matrix and
#' \eqn{a} the corresponding coefficient vector. Unless structural models are
#' be estimated, \eqn{A_0} is assumed to be an identity matrix.
#' 
#' If a vector is provided as argument \code{p} or \code{s}, the function will
#' produce a distinct model for all possible combinations of those specifications.
#' 
#' If \code{structural = TRUE}, the data matrix \eqn{Z_t} is augmented by negative
#' contemporaneous observations of endogenous variables, which correspond to
#' the coefficients in the lower triangular of \eqn{A_0}.
#' 
#' If \code{tvp} is \code{TRUE}, the respective coefficients
#' of the above model are assumed to be time varying. If \code{error} is \code{"sv"} or \code{"sv+covar"},
#' the error covariance matrix is assumed to be time varying.
#' 
#' Argument \code{error} specifies the structure of the covariance matrix of
#' the error term and how it is estimated. Possible specifications are:
#' \itemize{
#'  \item \code{"wishart"}: The covariance is estimated using a Wishart prior.
#'  \item \code{"gamma"}: Only the diagonal elements of the covariance matrix are estimated using a gamma prior.
#' Off-diagonal elements are not estimated and set to zero.
#'  \item \code{"gamma+covar"}: The diagonal elements of the covariance matrix are estimated using a gamma prior.
#' Covariances are estimated based on a triangular decomposition.
#'  \item \code{"sv"}: Only the diagonal elements of the covariance matrix are estimated using a stochastic volatility
#' algorithm. Off-diagonal elements are not estimated and set to zero.
#'  \item \code{"sv+covar"}: Only the diagonal elements of the covariance matrix are estimated using a stochastic volatility
#' algorithm. Covariances are estimated based on a triangular decomposition.
#'  \item \code{"ald"}: The errors are assumed to follow an asymmetric Laplace distribution, which
#' turns the model into a Bayesian quantile regression: the coefficients describe the conditional
#' quantile specified in argument \code{quantile} instead of the conditional mean. Off-diagonal
#' elements of the covariance matrix are not estimated and set to zero. See 'Details'.
#' }
#' 
#' Models with \code{error = "ald"} estimate a conditional quantile after Kozumi and Kobayashi (2011).
#' Minimising the quantile loss at \eqn{q} corresponds to maximising the likelihood of an asymmetric
#' Laplace distribution, which is a scale mixture of normal distributions. Conditional on the latent
#' scales of that mixture every equation is a weighted normal regression, which is what makes the
#' model a Gibbs sampler like the others.
#' 
#' Three properties of these models differ from the rest of the package. Covariances are not
#' estimated, since rotating the equations into each other leaves a residual whose quantile is not
#' the one that was asked for. A single quantile does not forecast, since the \eqn{h} step ahead quantile is
#' not the quantile of the iterated one step ahead quantiles. And the asymmetric Laplace is a working
#' likelihood rather than a claim about the data, so the posterior locates the quantile, but the
#' spread of the draws is not a calibrated credible interval without the adjustment of Yang et al.
#' (2016), which is not applied. Variable selection is available as \code{"bvs"}, not as
#' \code{"ssvs"}.
#'
#' With \code{quantile_grid = TRUE} the quantiles in \code{quantile}, at least two, form one
#' model rather than one model each: the structural quantile VAR of Chavleishvili and Manganelli
#' (2019). Every level is estimated -- each is the chain its single-quantile model would have
#' drawn -- and the draws of \code{posterior$a} and \code{posterior$u_scale} are stacked level by
#' level, the block of the first quantile first. With \code{structural = TRUE}, or a single
#' variable, the grid describes the whole distribution of each variable given the ones ordered
#' before it: the estimated quantiles are sorted (Chernozhukov et al., 2010), interpolated linearly
#' between the levels and continued by exponential tails. That distribution is what
#' \code{\link{add_posterior_forecasts}} simulates from, variable by variable, so a grid does
#' forecast. A scenario given there pins variables in forecast periods, and a
#' \code{forecast_quantile} takes every draw at one level, which gives the quantile paths the
#' impulse responses of the model are differences of. \code{\link{add_posterior_loglik}} scores
#' the density the grid describes. Only models with constant coefficients take a grid. Summaries,
#' plots and impulse responses work on the single levels, which
#' \code{\link{split_quantile_grid}} returns.
#' 
#' Available specifications for argument \code{varsel} are:
#' \itemize{
#'  \item \code{"none"}: No variable selection algorithm is used.
#'  \item \code{"bvs"}: Bayesian variable selection as proposed in Korobilis (2013).
#'  \item \code{"ssvs"}: Stochastic search variable selection as proposed in George et al. (2008).
#' }
#'
#' The one specification for argument \code{algorithm} is \code{"discount"}, the
#' discounted time varying parameter model of West & Harrison (1997, ch. 16)
#' with the discounted Wishart of Uhlig (1997), estimated by
#' \code{VarTvpDiscount}. It is not a sampler: its posterior is closed form --
#' one pass over the sample, no chain -- so \code{burnin} must be 0 and
#' \code{thin} 1, and \code{iterations} says only how many i.i.d. draws a
#' forecast takes from the answer. Its error covariance is the inverse Wishart
#' whole, so \code{error} must be \code{"wishart"} and neither variable
#' selection nor a structural model is available. What it buys is speed and an
#' exact marginal likelihood: the sum of \code{/posterior/loglik} is the log
#' marginal likelihood of the sample given the two discounts, so a grid over
#' them can be compared without a chain being run for any of it.
#'
#' \code{\link{add_posterior_coefficients}} estimates it like any other
#' algorithm, the filter being part of the vendored BayesTS core. It consumes no
#' random numbers, so two runs agree to the bit and a model estimated here and
#' the same model estimated by the \code{bayests} command line over a file
#' written with \code{\link{write_to_hdf5}} give the same numbers rather than
#' merely the same distribution. What comes back is a posterior rather than a
#' chain: one row per period under \code{posterior$a$mean},
#' \code{posterior$a$cov}, \code{posterior$u_sigma$scale} and
#' \code{posterior$df}, and no \code{coeffs} anywhere, because joining one
#' draw per period would look like a sampled path and is not one.
#'
#' Argument \code{iid} restricts the equations of the variables it names to
#' carry no coefficients at all. Such a variable is white noise: nothing
#' forecasts it, it forecasts nothing, and what the model is estimated for is
#' its contemporaneous correlation with the errors of the equations that do have
#' dynamics. That is how a high-frequency surprise becomes a variable of a
#' monthly VAR rather than an instrument beside one, following Jarocinski and
#' Karadi (2020).
#'
#' The restriction is exact rather than a tight prior on those coefficients: the
#' sampler never draws them, so they are zero in every draw. Everything else
#' treats the model as an ordinary VAR -- the restricted variables are still
#' regressors in the other equations, the error covariance still covers them,
#' and \code{\link{irf}} and \code{\link{fevd}} read the draws unchanged.
#'
#' \strong{The variables named in \code{iid} have to be the first columns of
#' \code{data}}, in any order among themselves. The restriction is carried by a
#' variable's position rather than by a list of positions, which is what lets
#' the sampler apply it without being told again which coefficients it owns. A
#' data set in another order is refused with a message naming the columns that
#' are there instead. It is available for models with constant coefficients and
#' cannot be combined with \code{structural} or \code{varsel}.
#'
#' With \code{missing = "estimate"} the panel does not have to be observed
#' whole. Every observation becomes one linear constraint on the panel, stored
#' in \code{data$train$constraints}, and every sweep of the sampler draws the
#' periods that were not observed from their distribution given the
#' constraints, the coefficients and the error covariance before it draws those
#' (Chan, Poon and Zhu 2023). A series observed at a lower frequency is a
#' constraint on several periods, the weights of which \code{aggregate} gives,
#' so that a quarterly series enters a monthly model as the average, sum or
#' growth rate of three months (Schorfheide and Song 2015). Which periods a
#' constraint reaches is decided by the dates of \code{data}, and one that
#' reaches before the estimation sample is left out with a message. The periods
#' before the sample that the first lags reach are taken as observed, with gaps
#' there filled by interpolation. What \code{data$train$y} holds where nothing
#' was observed is a starting value: a linear interpolation of what was.
#' The draws of the completed panel come back as \code{posterior$y}.
#'
#' A series named in \code{soft} holds its constraints up to a normal error
#' whose precision is estimated with a gamma prior, which
#' \code{\link{add_prior_options}} sets. It suits an aggregate whose weights
#' are an approximation, such as the growth rate of an average.
#'
#' Estimating what was not observed is available for \code{error = "wishart"},
#' \code{"gamma"}, \code{"gamma+covar"}, \code{"sv"} and \code{"sv+covar"},
#' with constant or time varying coefficients. It cannot be combined with
#' \code{structural} or the discounted model. It can be combined with
#' \code{iid}, for models with constant coefficients: a variable named there
#' that was not observed in a period -- a high-frequency surprise whose series
#' starts after the sample does, as in Jarocinski and Karadi (2020) -- is drawn
#' from its white noise given the errors of the other equations, and its gaps
#' start at zero rather than at an interpolation. Forecasts start from
#' each draw's completed panel. Forecasts are scored and conditioned on a
#' scenario in \code{data$forecast$constraints} only by models with constant
#' coefficients and \code{error = "wishart"} or \code{"gamma"}.
#'
#' @return An object of class 'bvarmodel' or, if a vector is given in \code{p}, \code{s}
#' or \code{quantile}, a list of class 'modellist' with one such object per
#' specification. A 'bvarmodel' is a list with the elements
#' \describe{
#'   \item{\code{data}}{a list with element \code{original}, which holds the time-series
#'   objects \code{endogen}, \code{exogen} and \code{deterministic}, and element
#'   \code{train}, which holds the estimation sample: \code{y}, a \eqn{T \times K}
#'   time-series object of the endogenous variables, \code{x}, a time-series object of
#'   the regressors, and \code{z}, the corresponding \eqn{TK} row matrix of regressors
#'   in SUR form, which is absent for the discounted model.}
#'   \item{\code{model}}{a list of the specification, including \code{type}
#'   (\code{"VAR"}), \code{algorithm}, the name of the posterior simulation algorithm,
#'   \code{k}, \code{p}, \code{m}, \code{s} and \code{n}, the numbers of endogenous
#'   variables, lags, exogenous variables, their lags and deterministic terms,
#'   \code{endogen}, the names of the endogenous variables, and \code{deterministic},
#'   \code{structural}, \code{error}, \code{varsel}, \code{tvp}, \code{iterations},
#'   \code{burnin}, \code{thin} if it is above 1, and, for \code{error = "ald"},
#'   \code{quantile} as specified.}
#' }
#' The later steps of the workflow add the elements \code{priors}, \code{initial} and
#' \code{posterior}.
#'
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 50, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' @references
#' 
#' Chan, J. C. C., Poon, A., & Zhu, D. (2023). High-dimensional conditionally
#' Gaussian state space models with missing data. \emph{Journal of
#' Econometrics, 236}(1), 105468. \doi{10.1016/j.jeconom.2023.05.005}
#'
#' Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). \emph{Bayesian Econometric Methods}
#' (2nd ed.). Cambridge: University Press.
#'
#' Chavleishvili, S., & Manganelli, S. (2019). Forecasting and stress testing with quantile
#' vector autoregression. \emph{ECB Working Paper}, 2330.
#'
#' Chernozhukov, V., Fernandez-Val, I., & Galichon, A. (2010). Quantile and probability curves
#' without crossing. \emph{Econometrica, 78}(3), 1093--1125.
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' Kozumi, H., & Kobayashi, G. (2011). Gibbs sampling methods for Bayesian quantile regression.
#' \emph{Journal of Statistical Computation and Simulation, 81}(11), 1565--1578.
#' \doi{10.1080/00949655.2010.496117}
#' 
#' Schorfheide, F., & Song, D. (2015). Real-time forecasting with a
#' mixed-frequency VAR. \emph{Journal of Business & Economic Statistics,
#' 33}(3), 366--380. \doi{10.1080/07350015.2014.954707}
#'
#' Lütkepohl, H. (2006). \emph{New Introduction to Multiple Time Series Analysis} (2nd ed.). Berlin: Springer.
#' 
#' Uhlig, H. (1997). Bayesian vector autoregressions with stochastic volatility.
#' \emph{Econometrica, 65}(1), 59--73. \doi{10.2307/2171813}
#' 
#' West, M., & Harrison, J. (1997). \emph{Bayesian forecasting and dynamic models}
#' (2nd ed.). New York: Springer.
#' 
#' @seealso \code{\link{bvartools_model}} describes the object this returns, element by element.
#' @family model set-up
#' @export
create_bvarmodel <- function(data, p = 2,
                             exogen = NULL, s = 2,
                             deterministic = "const",
                             seasonal = FALSE,
                             iid = NULL,
                             structural = FALSE,
                             error = "wishart",
                             quantile = 0.5,
                             tvp = FALSE,
                             varsel = "none",
                             algorithm = NULL,
                             delta_beta = 1,
                             delta_sigma = 1,
                             iterations = 20000,
                             burnin = 2000,
                             thin = 1,
                             missing = "omit",
                             aggregate = NULL,
                             soft = NULL,
                             quantile_grid = FALSE) {
  
  # Input checks ----
  if (!"ts" %in% class(data)) {
    stop("Argument 'data' must be an object of class 'ts'.")
  }
  
  if (!is.null(exogen)) {
    if (!"ts" %in% class(exogen)) {
      stop("Argument 'exogen' must be an object of class 'ts'.")
    }
    if (!is.numeric(s) || length(s) == 0 || anyNA(s) || any(s < 0 | s %% 1 != 0)) {
      stop("Argument 's' must contain non-negative integers when 'exogen' is specified.")
    }
  }
  
  if (seasonal & !deterministic %in% c("const", "both")) {
    stop("Argument 'deterministic' must be either 'const' or 'both' when using 'seasonal = TRUE'.")
  }
  
  if ("character" %in% class(error)) {
    if (!error %in% c("wishart", "gamma", "gamma+covar", "sv", "sv+covar", "ald")) {
      stop("Invalid specification of argument 'error'.")
    }
  } else {
    stop("Argument 'error' must be of class 'character'.")
  }
  
  # The quantile is the one specification that produces several models on its
  # own, so it is checked here and looped over below. It is meaningless for
  # every model that is not an asymmetric Laplace one, which is why an object
  # of theirs does not carry it at all.
  if (!is.logical(quantile_grid) || length(quantile_grid) != 1 || is.na(quantile_grid)) {
    stop("Argument 'quantile_grid' must be TRUE or FALSE.")
  }
  if (quantile_grid) {
    if (error != "ald" || tvp || !is.null(algorithm)) {
      stop("A grid of quantiles needs error = \"ald\" and constant coefficients.")
    }
    if (length(quantile) < 2 || anyDuplicated(quantile) > 0) {
      stop("Argument 'quantile' must contain at least two distinct values when ",
           "'quantile_grid' is TRUE.")
    }
    quantile <- sort(quantile)
  }
  if (error == "ald") {
    if (!"numeric" %in% class(quantile)) {
      stop("Argument 'quantile' must be of class 'numeric'.")
    }
    if (length(quantile) == 0 | any(is.na(quantile))) {
      stop("Argument 'quantile' must contain at least one non-missing value.")
    }
    if (any(quantile <= 0 | quantile >= 1)) {
      stop("Argument 'quantile' must contain values between 0 and 1.")
    }
  }
  
  # Endogenous variables ----
  if (is.null(dimnames(data))) {
    # If 'data' is a simple ts object, transform it into a matrix object
    # to keep variable name information
    tsp_temp <- stats::tsp(data)
    data <- stats::ts(as.matrix(data), class = c("mts", "ts", "matrix"))
    stats::tsp(data) <- tsp_temp
    dimnames(data)[[2]] <- "y"
  }
  
  # Exogenous variables ----
  use_exo <- !is.null(exogen)
  if (use_exo) {
    if (is.null(dimnames(exogen))) {
      # If 'exogen' is a simple ts object, transform it into a matrix object
      # to keep variable name information
      tsp_temp <- stats::tsp(exogen)
      exogen <- stats::ts(as.matrix(exogen), class = c("mts", "ts", "matrix"))
      stats::tsp(exogen) <- tsp_temp
      dimnames(exogen)[[2]] <- "exogen"
    }
  }
  
  if (NCOL(data) == 1 & structural) {
    structural <- FALSE
    if (error == "gamma+covar") {
      error <- "gamma"
    }
    if (error == "sv+covar") {
      error <- "sv"
    }
  }
  
  if (structural & error %in% c("wishart", "gamma+covar", "sv+covar")) {
    stop(paste0("Structural models cannot be estimated with argument 'error' specified as '", error,"'."))
  }
  
  if (!varsel %in% c("none", "bvs", "ssvs")) {
    stop("Specification of argument 'varsel' is not supported.")
  }
  
  # Refused where the specification is made rather than where the priors are
  # added, since a quantile regression model with SSVS is not a model this
  # package has at all.
  if (error == "ald" & varsel == "ssvs") {
    stop("Variable selection algorithm 'ssvs' is not available for a quantile ",
         "regression model. Consider using 'bvs' instead.")
  }
  
  # The panel as it was observed, and the same panel with every gap filled
  # for the regressors to be built from. See .panel_constraints().
  if (!is.character(missing) || length(missing) != 1 || !missing %in% c("omit", "estimate")) {
    stop("Argument 'missing' must be either 'omit' or 'estimate'.")
  }
  if (!is.null(aggregate)) {
    missing <- "estimate"
  }
  if (!is.null(soft) && missing != "estimate") {
    stop("Argument 'soft' needs missing = 'estimate'.")
  }
  observed <- NULL
  if (missing == "estimate") {
    if (error == "ald" || !is.null(algorithm) || structural) {
      stop("Estimating what was not observed is not available for a quantile regression ",
           "model, the discounted model or a structural model.")
    }
    panel <- .check_panel_arguments(data, aggregate, soft)
    aggregate <- panel[["aggregate"]]
    observed <- data
    data <- .fill_panel(data, aggregate, zero = iid)
  }

  # A lag order of 1.5 failed with "subscript out of bounds" and one of -1 was
  # taken as zero.
  if (!is.numeric(p) || length(p) == 0 || anyNA(p) || any(p != round(p)) || any(p < 0)) {
    stop("Argument 'p' must be a vector of non-negative whole numbers.")
  }

  data_name <- dimnames(data)[[2]]
  k <- NCOL(data)
  p_max <- max(p)
  
  model_type <- NULL
  if (tvp) {
    model_type <- paste0(model_type, "Tvp")
  } else {
    model_type <- paste0(model_type, "Normal")
  }
  if (error == "wishart") {
    model_type <- paste0(model_type, "Wishart")
  }
  if (error %in% c("gamma", "gamma+covar")) {
    model_type <- paste0(model_type, "Gamma")
  }
  if (error %in% c("sv", "sv+covar")) {
    model_type <- paste0(model_type, "Stochvol")
  }
  if (error == "ald") {
    model_type <- paste0(model_type, "Ald")
  }
  model_type <- paste0("Var", model_type)

  # The one algorithm outside the naming grammar for a VAR. See
  # R/discount_models.R for what it is and where the rest of the package has to
  # be told about it.
  if (!is.null(algorithm)) {
    if (!identical(algorithm, "discount")) {
      stop("Specified algorithm not recognized.")
    }
    model_type <- "VarTvpDiscount"
    # Refused here rather than by BayesTS on the written file: for a grid of
    # models that is one error instead of one per file.
    .check_discount_specification(error, varsel, structural, burnin, thin)
    .check_discount_deltas(delta_beta, delta_sigma, k)
    # The coefficients of this model drift, which is what `tvp` says. What
    # governs how far they drift is `delta_beta`, and it is a model in its own
    # right at one, where they do not move at all.
    tvp <- TRUE
  }
  use_discount <- identical(model_type, "VarTvpDiscount")

  model <- NULL
  model[["type"]] <- ifelse(length(data_name) == 1, "AR", "VAR")
  model[["algorithm"]] <- model_type
  model[["k"]] <- length(data_name)
  model[["p"]] <- 0L
  model[["m"]] <- 0L
  model[["s"]] <- 0L
  model[["n"]] <- 0L
  model[["varsel"]] <- varsel
  model[["endogen"]] <- dimnames(data)[[2]]
  model[["n_iid"]] <- .check_iid_variables(iid, data_name, structural, varsel, tvp)
  if (model[["n_iid"]] == 0L) {
    model[["n_iid"]] <- NULL
  }
  if (use_exo) {
    model[["exogen"]] <- dimnames(exogen)[[2]]
  }
  
  temp <- data
  temp_name <- data_name
  if (p_max >= 1) {
    # Obtain lags of endogenous variables
    for (i in 1:p_max) {
      temp <- cbind(temp, stats::lag(data, -i))
      if (nchar(p_max) > 2) {
        i_temp <- paste0(c(rep(0, nchar(p_max) - nchar(i)), i), collapse = "")
      } else {
        i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
      }
      temp_name <- c(temp_name, paste0(data_name, ".", i_temp))
    }
  }
  
  # Exogenous variables ---- 
  if (use_exo) {
    exo_name <- model[["exogen"]]
    m <- length(exo_name)
    s_max <- max(s)
    
    temp <- cbind(temp, exogen)
    if (nchar(s_max) > 2) {
      i_temp <- rep(0, nchar(s_max))
    } else {
      i_temp <- rep(0, 2)
    }
    i_temp <- paste0(i_temp, collapse = "")
    temp_name <- c(temp_name, paste0(exo_name, ".l", i_temp))
    if (s_max > 0) {
      for (i in 1:s_max) {
        temp <- cbind(temp, stats::lag(exogen, -i))
        if (nchar(s_max) > 2) {
          i_temp <- paste0(c(rep(0, nchar(s_max) - nchar(i)), i), collapse = "")
        } else {
          i_temp <- paste0(c(rep(0, 2 - nchar(i)), i), collapse = "")
        }
        temp_name <- c(temp_name, paste0(exo_name, ".l", i_temp))
      } 
    }
    
    model[["m"]] <- length(exo_name)
    model[["s"]] <- 0
    
  } else {
    s <- 0
    s_max <- 0
    m <- 0L
  }
  
  tt <- nrow(temp)
  det_data <- NULL
  det_name <- NULL
  det_pos <- ncol(temp)
  
  # Add intercept term
  if (deterministic %in% c("const", "both")) {
    temp <- cbind(temp, 1)
    temp_name <- c(temp_name, "const")
    det_name <- c(det_name, "const")
  }
  
  # Add linear trend
  if (deterministic %in% c("trend", "both")) {
    # One in the first period that is estimated on, which is the first row with
    # every lag available. It used to be the row after max(p, s), which is that
    # period only when the data and the exogenous series start together:
    # 'temp' spans both, so an 'exogen' reaching further back shifted the trend
    # -- and with "both" the meaning of the intercept -- by as many periods.
    first <- which(stats::complete.cases(temp))[1]
    temp <- cbind(temp, seq_len(tt) - first + 1)
    temp_name <- c(temp_name, "trend")
    det_name <- c(det_name, "trend")
  }
  
  # Add seasonal dummies
  if (seasonal) {
    freq <- stats::frequency(data)
    if (freq == 1) {
      warning("The frequency of the provided data is 1. No seasonal dummmies are generated.")
    } else {
      pos <- which(stats::cycle(temp) == 1)[1]
      pos <- rep(1:freq, 2)[pos:(pos + (freq - 2))]
      for (i in 1:(freq - 1)) {
        s_temp <- rep(0, freq)
        s_temp[pos[i]] <- 1
        temp <- cbind(temp, rep(s_temp, length.out = tt))
        temp_name <- c(temp_name, paste("season.", i, sep = ""))
        det_name <- c(det_name, paste("season.", i, sep = ""))
      }
    }
  }
  
  # Update model specs for deterministic terms
  use_det <- FALSE
  if (length(det_name) > 0) {
    model[["n"]] <- length(det_name)
    model[["deterministic"]] <- det_name
    use_det <- TRUE
    det_data <- temp[, det_pos + 1:model[["n"]]]
    
    # If 'det_data' is a simple ts object, transform it into a matrix object
    # to keep variable name information
    tsp_det_data <- stats::tsp(det_data)
    det_data <- stats::ts(as.matrix(det_data), class = c("mts", "ts", "matrix"))
    stats::tsp(det_data) <- tsp_det_data
    dimnames(det_data) <- list(NULL, det_name)
  }
  
  temp <- stats::na.omit(temp)
  
  # Set if the model is structural
  if ("logical" %in% class(structural)) {
    model[["structural"]] <- structural
    if (structural) {
      model[["type"]] <- "SVAR" 
    }
  } else {
    stop("Argument 'structural' must be of class 'logical'.")
  }
  
  ## errors ----
  if ("character" %in% class(error)) {
    if (!error %in% c("wishart", "gamma", "gamma+covar", "sv", "sv+covar", "ald")) {
      stop("Invalid specification of argument 'error'.")
    }
    model[["error"]] <- error
  } else {
    stop("Argument 'error' must be of class 'character'.")
  }
  
  ## tvp ----
  if ("logical" %in% class(tvp)) {
    model[["tvp"]] <- tvp
  } else {
    stop("Argument 'tvp' must be of class 'logical'.")
  }
  
  # Iterations, burnin and thinning ----
  model[["iterations"]] <- as.integer(iterations)
  model[["burnin"]] <- as.integer(burnin)
  # Carried only when it thins, as `quantile` is carried only by the models that
  # read it: a model without it keeps every draw.
  model[["thin"]] <- .check_sampler_thin(thin)
  
  # Data that is equal across models ----
  
  # Endogenous variables y
  y <- stats::ts(as.matrix(temp[, 1:k]), class = c("mts", "ts", "matrix"))
  stats::tsp(y) <- stats::tsp(temp)
  dimnames(y)[[2]] <- temp_name[1:k]
  
  # What was observed of the sample, where it was not observed whole
  constraints <- NULL
  if (!is.null(observed)) {
    constraints <- .panel_constraints(observed, y, aggregate, soft)
  }

  # Structural data
  y_A0 <- NULL
  if (structural & k > 1) {
    y_A0 <- kronecker(-y, diag(1, k))
    pos <- NULL
    for (j in 1:k) {
      pos <- c(pos, (j - 1) * k + 1:j)
    }
    y_A0 <- y_A0[, -pos]
  }
  
  # Create model list ----
  # A vector of quantiles is a list of models, in the same way a grid of lag orders
  # is: one quantile per model, so the grid parallelises without the sampler
  # knowing about it. Models that are not asymmetric Laplace ones pass through
  # this loop once and never see the field.
  # A grid of quantiles is the one exception: its levels are one model, so the
  # loop runs once and the model carries all of them.
  quantiles <- if (error == "ald" && !quantile_grid) quantile else NA_real_
  if (quantile_grid) {
    model[["quantiles"]] <- quantile
  }
  
  # A grid over the two discounts, for the one algorithm that has them. Every
  # other model gets the single specification it always had, so the two inner
  # loops run once and change nothing. See create_bvecmodel() for why a grid
  # over them is worth having.
  if (use_discount) {
    grid_beta <- as.numeric(delta_beta)
    grid_sigma <- as.numeric(delta_sigma)
  } else {
    grid_beta <- NA_real_
    grid_sigma <- NA_real_
  }

  result <- NULL
  for (i in p) { # for each lag p
    for (j in s) { # for each lag s
     for (q in quantiles) { # for each quantile
      for (d_beta in grid_beta) {
       for (d_sigma in grid_sigma) {
      pos <- NULL
      model_i <- model
      if (error == "ald" && !quantile_grid) {
        model_i[["quantile"]] <- q
      }
      if (use_discount) {
        model_i[["delta_beta"]] <- d_beta
        model_i[["delta_sigma"]] <- d_sigma
      }
      if (i >= 1) {
        pos <- c(pos, k + 1:(k * i))
        model_i[["p"]] <- as.integer(i)
      }  
      if (use_exo) {
        pos <- c(pos, k + k * p_max + 1:(m * (j + 1)))
        model_i[["s"]] <- as.integer(j)
      }
      if (use_det) {
        pos <- c(pos, k + k * p_max + m * (s_max + 1) + 1:length(det_name))
      }
      
      x <- NULL
      z <- NULL
      if (length(pos) > 0) {
        # Create data input matrix of the respective model
        x <- stats::ts(as.matrix(temp[, pos]), class = c("mts", "ts", "matrix")) 
        stats::tsp(x) <- stats::tsp(temp)
        dimnames(x)[[2]] <- temp_name[pos]
        z <- kronecker(x, diag(1, k))
      }
      
      # If specified add structural data to SUR form
      if (!is.null(y_A0)) {
        z <- cbind(z, y_A0)
      }
      dimnames(z) <- NULL

      # Not kept for a discounted model, which reads the compact regressors;
      # see create_bvecmodel().
      if (use_discount) {
        z <- NULL
      }
      
      # Create individual model
      result_i <- list("model" = model_i,
                       "data" = list("original" = .drop_null(list("endogen" = if (is.null(observed)) data else observed,
                                                       "exogen" = exogen,
                                                       "deterministic" = det_data)),
                                     "train" = .drop_null(list("y" = y,
                                                    "x" = x,
                                                    "z" = z,
                                                    "constraints" = constraints))))
      
      # Update class of individual model
      class(result_i) <- append("bvarmodel", class(result_i)) 
      
      result <- c(result, list(result_i)) 
      
       }
      }
     }
    }
  }
  
  if (length(result) == 1) {
    result <- result[[1]]
  } else {
    class(result) <- append("modellist", class(result)) 
  }
  
  return(result)
}
