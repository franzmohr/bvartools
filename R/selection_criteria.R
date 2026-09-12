#' Model Selection Criteria
#'
#' Generic function used to calculate selection criteria.
#'
#' @param object an object with suitable input data passed forward to method.
#' @param ... arguments passed forward to method.
#'
#' @details The in-sample criteria are obtained from the draws of the
#' log-likelihood that \code{\link{add_posterior_loglik}} adds to a model. With
#' \eqn{ll_t^{(i)}} the log-likelihood of period \eqn{t} in draw \eqn{i} of
#' \eqn{R} draws, \eqn{LL^{(i)} = \sum_{t = 1}^{T} ll_t^{(i)}} and
#' \eqn{\kappa} the number of estimated parameters, these are
#' \itemize{
#'  \item{\code{"LL"}: the log-likelihood \eqn{LL^{(i)}}, which is summarised by
#' the mean, the median and the bounds of the credible band of its draws.}
#'  \item{\code{"AIC"}: \eqn{D + 2 \kappa};}
#'  \item{\code{"BIC"}: \eqn{D + \ln(T) \kappa};}
#'  \item{\code{"HQ"}: \eqn{D + 2 \ln(\ln(T)) \kappa};}
#'  \item{\code{"WAIC"}: \eqn{-2 \sum_{t = 1}^{T} \left( \ln \left(
#' \frac{1}{R} \sum_{i = 1}^{R} \exp(ll_t^{(i)}) \right) -
#' \mathrm{Var}_i \left[ ll_t^{(i)} \right] \right)};}
#'  \item{\code{"LOOIC"}: \eqn{-2} times the expected log pointwise predictive
#' density of leave-one-out cross validation, obtained by Pareto smoothed
#' importance sampling.}
#' }
#'
#' \eqn{D} is the deviance of the model at its point estimate. Since AIC, BIC
#' and HQ correct the deviance of a fitted model for the optimism of having
#' fitted it, it is the deviance at the point estimate that their penalties
#' belong to. The mean of the deviance over the posterior is the larger
#' quantity, by about the effective number of parameters, so adding a penalty
#' to it would charge the complexity of the model a second time and tilt every
#' comparison towards the smaller model. \eqn{D} is therefore recovered from
#' the draws as their mean deviance less
#' \eqn{\sum_{t = 1}^{T} \mathrm{Var}_i [ll_t^{(i)}]}, the effective number of
#' parameters that WAIC also penalises with. Evaluated at the maximum
#' likelihood estimates AIC reduces to \eqn{T} times the expression in
#' Luetkepohl (2006) plus terms that do not depend on the lag order, so both
#' order a set of lag orders in the same way.
#'
#' AIC, BIC and HQ are point estimates rather than quantities with a posterior
#' distribution, so their credible bands are \code{NA}. WAIC and LOOIC are
#' point estimates as well, and their bands are normal intervals built from
#' their standard errors, which is the usual way they are reported.
#'
#' The choice between the criteria is one of what \eqn{\kappa} means for the
#' model at hand. A count of parameters describes a model with constant
#' coefficients and a weak prior, which is the case AIC, BIC and HQ are derived
#' for and the one in which they reproduce the textbook lag order selection. It
#' does not describe a model whose coefficients or variances follow a state
#' equation, where the prior on the state variances decides how much of the
#' nominal freedom is used, and it does not describe a shrinkage prior, which
#' buys back degrees of freedom that \eqn{\kappa} does not see. WAIC and LOOIC
#' penalise by the flexibility the fit actually used and remain defined in
#' those cases, so they are the criteria to compare constant, time varying and
#' stochastic volatility specifications with.
#'
#' LOOIC estimates the same quantity as WAIC by reweighting the posterior
#' towards the one that has not seen a period, which is the more accurate route
#' when it works and which says when it does not: the shape parameter of the
#' Pareto tail fitted to the importance ratios of a period exceeds its
#' threshold when that period is too influential to be reweighted, and the
#' print methods report how many periods this happened for.
#'
#' All criteria require that the models that are compared were estimated on the
#' same observations. \code{\link{create_bvarmodel}} and
#' \code{\link{create_bvecmodel}} ensure this for a vector of lag orders by
#' trimming the data by the maximum lag, and \code{\link{align_model_obs}}
#' restricts models that were created separately to their common sample.
#'
#' @references
#'
#' Luetkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#'
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using
#' leave-one-out cross-validation and WAIC. \emph{Statistics and Computing, 27}(5), 1413--1432.
#' \doi{10.1007/s11222-016-9696-4}
#'
#' Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024). Pareto smoothed
#' importance sampling. \emph{Journal of Machine Learning Research, 25}(72), 1--58.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable
#' information criterion in singular learning theory. \emph{Journal of Machine Learning
#' Research, 11}, 3571--3594.
#'
#' Zhang, J., & Stephens, M. A. (2009). A new and efficient estimation method for the generalized
#' Pareto distribution. \emph{Technometrics, 51}(3), 316--325. \doi{10.1198/tech.2009.08017}
#'
#' @export
selection_criteria <- function (object, ...) {
 UseMethod("selection_criteria")
}
