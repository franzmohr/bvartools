# Model Selection Criteria

Generic function used to calculate selection criteria.

## Usage

``` r
selection_criteria(object, ...)

# S3 method for class 'selcrit'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'selcritlist'
print(x, digits = max(3L, getOption("digits") - 3L), relative = 0, ...)
```

## Arguments

- object:

  an object with suitable input data passed forward to method.

- ...:

  arguments passed forward to method.

- x:

  an object used to select a method. Usually, the result of a call to
  `selection_criteria`.

- digits:

  the minimum number of significant digits to be printed in values.

- relative:

  an integer specifying the model that is used as the reference of
  relative forecast performance. Default is `0`, which indicates that
  results are not displayed in relation to each other.

## Value

The value returned by the method for the class of `object`, as described
on the pages of the methods.

## Details

The in-sample criteria are obtained from the draws of the log-likelihood
that
[`add_posterior_loglik`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.md)
adds to a model. With \\ll_t^{(i)}\\ the log-likelihood of period \\t\\
in draw \\i\\ of \\R\\ draws, \\LL^{(i)} = \sum\_{t = 1}^{T}
ll_t^{(i)}\\ and \\\kappa\\ the number of estimated parameters, these
are

- `"LL"`: the log-likelihood \\LL^{(i)}\\, which is summarised by the
  mean, the median and the bounds of the credible band of its draws.

- `"AIC"`: \\D + 2 \kappa\\;

- `"BIC"`: \\D + \ln(T) \kappa\\;

- `"HQ"`: \\D + 2 \ln(\ln(T)) \kappa\\;

- `"WAIC"`: \\-2 \sum\_{t = 1}^{T} \left( \ln \left( \frac{1}{R}
  \sum\_{i = 1}^{R} \exp(ll_t^{(i)}) \right) - \mathrm{Var}\_i \left\[
  ll_t^{(i)} \right\] \right)\\;

- `"LOOIC"`: \\-2\\ times the expected log pointwise predictive density
  of leave-one-out cross validation, obtained by Pareto smoothed
  importance sampling.

\\D\\ is the deviance of the model at its point estimate. Since AIC, BIC
and HQ correct the deviance of a fitted model for the optimism of having
fitted it, it is the deviance at the point estimate that their penalties
belong to. The mean of the deviance over the posterior is the larger
quantity, by about the effective number of parameters, so adding a
penalty to it would charge the complexity of the model a second time and
tilt every comparison towards the smaller model. \\D\\ is therefore the
deviance at the posterior mean of the parameters, obtained by evaluating
the log-likelihood of the model once at that point. For a VEC model the
point is the matrix of the rank of the model that is closest to the
posterior mean of \\\Pi = \alpha \beta^{\prime}\\ in the metric of the
likelihood, i.e. that minimises \\tr(\bar{Q} \Delta W^{\prime} W
\Delta^{\prime})\\, where \\\Delta\\ is its difference to the posterior
mean, \\\bar{Q}\\ the posterior mean of the error precision and \\W\\
the error correction term, since \\\alpha\\ and \\\beta\\ are identified
only up to a rotation and their own means say nothing. The point does
not depend on the scale of the series. Under a flat prior \\D\\ is close
to the deviance at the maximum likelihood estimates. Evaluated at the
maximum likelihood estimates AIC reduces to \\T\\ times the expression
in Luetkepohl (2006) plus terms that do not depend on the lag order, so
both order a set of lag orders in the same way.

AIC, BIC and HQ are point estimates rather than quantities with a
posterior distribution, so their credible bands are `NA`. WAIC and LOOIC
are point estimates as well, and their bands are normal intervals built
from their standard errors, which is the usual way they are reported.

The choice between the criteria is one of what \\\kappa\\ means for the
model at hand. A count of parameters describes a model with constant
coefficients and a weak prior, which is the case AIC, BIC and HQ are
derived for and the one in which they reproduce the textbook lag order
selection. It does not describe a model whose coefficients or variances
follow a state equation, where the prior on the state variances decides
how much of the nominal freedom is used, and it does not describe a
shrinkage prior, which buys back degrees of freedom that \\\kappa\\ does
not see. WAIC and LOOIC penalise by the flexibility the fit actually
used and remain defined in those cases.

LOOIC estimates the same quantity as WAIC by reweighting the posterior
towards the one that has not seen a period, which is the more accurate
route when it works and which says when it does not: the shape parameter
of the Pareto tail fitted to the importance ratios of a period exceeds
its threshold when that period is too influential to be reweighted, and
the print methods report how many periods this happened for.

For a model whose coefficients or variances follow a state equation, the
pointwise log-likelihood that WAIC and LOOIC are built from is evaluated
at the states of each period, which have seen that period's observation.
The state of a period is informed mostly by its own observation, so
leaving the observation out by reweighting fails for exactly the periods
a path bends towards, and neither criterion separates models that differ
in how much of that freedom they have, such as the cointegration ranks
of a time varying VEC model. The criterion for those comparisons is
`"LPL"`, the log predictive likelihood, whose densities condition only
on the data before the period they evaluate. Koop, León-González and
Strachan (2011) choose between time varying cointegration models with
it.

Two things produce those densities and both report them as `"LPL"`.
[`add_predictive_loglik`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md)
takes one per window of an expanding window exercise, carrying the
states forward with their state equations, which is the form the
comparison of ranks above is made in. A single model carries one per
horizon of its forecast in `posterior$forecast$loglik`, written by
BayesTS against the values in `data$test$y`, each conditioning on the
periods realised before it. They are the same quantity computed two ways
and are summarised by the same code, so the same densities give the same
number.

All criteria require that the models that are compared were estimated on
the same observations.
[`create_bvarmodel`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
and
[`create_bvecmodel`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md)
ensure this for a vector of lag orders by trimming the data by the
maximum lag, and
[`align_model_obs`](https://franzmohr.github.io/bvartools/reference/align_model_obs.md)
restricts models that were created separately to their common sample.

## References

Koop, G., León-González, R., & Strachan, R. W. (2011). Bayesian
inference in a time varying cointegration model. *Journal of
Econometrics, 165*(2), 210–220.
[doi:10.1016/j.jeconom.2011.07.007](https://doi.org/10.1016/j.jeconom.2011.07.007)

Luetkepohl, H. (2006). *New introduction to multiple time series
analysis* (2nd ed.). Berlin: Springer.

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing, 27*(5), 1413–1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research, 25*(72), 1–58.

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
*Journal of Machine Learning Research, 11*, 3571–3594.

Zhang, J., & Stephens, M. A. (2009). A new and efficient estimation
method for the generalized Pareto distribution. *Technometrics, 51*(3),
316–325.
[doi:10.1198/tech.2009.08017](https://doi.org/10.1198/tech.2009.08017)

## See also

Methods:
[`selection_criteria.bvarmodel`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md),
[`selection_criteria.bvecmodel`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md),
[`selection_criteria.expandingwindow`](https://franzmohr.github.io/bvartools/reference/selection_criteria.expandingwindow.md),
[`selection_criteria.externalforecast`](https://franzmohr.github.io/bvartools/reference/selection_criteria.externalforecast.md),
[`selection_criteria.modellist`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md).
