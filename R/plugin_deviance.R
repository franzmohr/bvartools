# Deviance of a model at its point estimate, recovered from the posterior
# draws of the log-likelihood.
#
# AIC, BIC and HQ are not penalties bolted onto an arbitrary measure of fit.
# Each estimates how the model would do on data it has not seen, and each
# corrects the deviance of a *fitted* model for the optimism of having fitted
# it. The correction is derived for a likelihood that was maximised, and it is
# the deviance at that maximum -- or, here, at the Bayesian point estimate --
# that the correction belongs to.
#
# The mean of the deviance over the posterior is a different quantity. It is
# not optimistic but pessimistic, because it averages over draws that are away
# from the point estimate, and it exceeds the deviance at the point estimate
# by about the effective number of parameters:
#
#   E[D(theta)] = D(theta_hat) + p_D
#
# Adding the penalty to the posterior mean therefore charges the complexity of
# the model once through the penalty and a second time through the averaging.
# The size of the mistake is not cosmetic: on the data of Luetkepohl (2006) it
# amounts to half the penalty again and reverses the choice of the lag order.
#
# The deviance at the point estimate is recovered by undoing the averaging,
# with the effective number of parameters estimated by the variance of the
# pointwise log-likelihood across draws -- the same quantity WAIC penalises
# with, so no second estimator is introduced. This keeps every criterion a
# function of the stored draws of the log-likelihood alone, which is what
# makes it work unchanged for the error correction, time varying and
# stochastic volatility models, where the residuals cannot be reconstructed
# from a single matrix of regressors.
.plugin_deviance <- function(loglik) {

  loglik <- as.matrix(loglik)
  mean_deviance <- -2 * sum(colMeans(loglik))

  # A single draw is its own point estimate and has no spread to correct for
  if (nrow(loglik) < 2) {
    return(mean_deviance)
  }

  mean_deviance - sum(apply(loglik, 2, stats::var))
}

# The criteria that are point estimates rather than quantities with a
# posterior distribution are reported in the same four columns as the others,
# with the estimate in both the mean and the median and no credible band,
# which they do not have: a criterion is a function of the data and the
# estimator, not a parameter to be uncertain about.
.point_criterion <- function(value) {
  data.frame("mean" = value,
             "median" = value,
             "qlower" = NA_real_,
             "qupper" = NA_real_)
}
