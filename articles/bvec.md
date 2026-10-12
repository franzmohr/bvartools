# Bayesian Error Correction Models with Priors on the Cointegration Space

## Introduction

This vignette provides the code to set up and estimate a Bayesian vector
error correction (BVEC) model with the `bvartools` package. The
presented Gibbs sampler is based on the approach of Koop et al. (2010),
who propose a prior on the cointegration space. The estimated model has
the following form

``` math
\Delta y_t = \Pi y_{t - 1} + \sum_{l = 1}^{p - 1} \Gamma_l \Delta y_{t - l} + C d_t + u_t,
```
where $`\Pi = \alpha \beta^{\prime}`$ with cointegration rank $`r`$,
$`u_t \sim N(0, \Sigma)`$ and $`d_t`$ contains deterministic terms. For
an introduction to vector error correction models see
<https://www.r-econometrics.com/timeseries/vecintro/>.

## Data

To illustrate the workflow of the analysis, the Austrian long-term
interest rate (`lr`) and inflation (`Dp`) are taken from data set
`at_macrodata`, both in percent per quarter from 1979Q2 to 2019Q4. The
series are seasonally adjusted, so that $`d_t`$ contains only an
intercept.

``` r

library(bvartools)

data("at_macrodata")

data <- window(at_macrodata[["domestic"]][, c("lr", "Dp")], end = c(2019, 4)) * 100
plot(data) # Plot the series
```

![plot of chunk data](figures/bvec-data-1.png)

plot of chunk data

The `create_bvecmodel` function produces an object of class `bvecmodel`,
where the element `data$train` contains the data matrices `y`, `w` and
`x` for the BVEC estimator, where `y` is the matrix of dependent
variables, `w` is a matrix of potentially cointegrated regressors, and
`x` is the matrix of non-cointegration regressors.

``` r

object <- create_bvecmodel(data, p = 4, r = 1,
                           const = "unrestricted",
                           iterations = 5000, burnin = 1000)
```

Argument `p` represents the lag order of the VAR form of the model and
`r` is the cointegration rank of `Pi`. Function `create_bvecmodel`
requires to specify the inclusion of intercepts, trends and seasonal
dummies separately. This allows to decide on whether they enter the
cointegration term (`"restricted"`) or the non-cointegration part
(`"unrestricted"`) of the model.

Models such as the one proposed in Koop, Leon-Gonzalez and Strachan
(2010) put a restriction on the space spanned by $`\beta`$. This means,
that the scale of the variables in the error correction term becomes
important. A reasonable approach to handle this, is to rescale the
series in the error correction term. Function
`scale_error_correction_term` can be used for this purpose. It rescales
all non-deterministic elements by their respective differenced series’
standard deviation. If `object[["data"]][["train"]][["w"]]` contains a
column named `trend`, its valus are replaced by $`(t - T/2) / T`$.

``` r

object <- scale_error_correction(object)
```

## Priors

Function `add_priors` adds the necessary prior specifications to object
`object`. For the current application non-informative priors are used:

``` r

object <- add_priors(object,
                     coef = list(v_i = 0, v_i_det = 0),
                     coint = list(v_i = 0, p_tau_i = 1),
                     sigma = list(df = 1, scale = .0001))
```

### Adding initial values

Function `add_initial_values` generates of first draws of the Gibbs
sampler for the specified model(s) in object `model` and augments the
object accordingly.

``` r

object <- add_initial_values(object, method = "maxlik")
```

## Estimation

### User-specific algorithm

The following code produces posterior draws using the algorithm
described in Koop et al. (2010).

``` r

# Reset random number generator for reproducibility
set.seed(7654321)

# Obtain data matrices
y <- t(object[["data"]][["train"]][["y"]])
w <- t(object[["data"]][["train"]][["w"]])
x <- t(object[["data"]][["train"]][["x"]])

r <- object[["model"]][["rank"]] # Set rank

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
k_w <- nrow(w) # Number of regressors in error correction term
k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms
k_gamma <- k * k_x # Total number of non-cointegration coefficients

k_alpha <- k * r # Number of elements in alpha
k_beta <- k_w * r # Number of elements in beta

# Priors
a_mu_prior <- object[["priors"]][["a"]][["mu"]] # Prior means
a_v_i_prior <- object[["priors"]][["a"]][["v_inv"]] # Inverse of the prior covariance matrix

v_i <- object[["priors"]][["beta"]][["v_inv"]]
p_tau_i <- object[["priors"]][["beta"]][["p_tau_inv"]]

sigma_df_prior <- object[["priors"]][["u_sigma"]][["df"]] # Prior degrees of freedom
sigma_scale_prior <- object[["priors"]][["u_sigma"]][["scale"]] # Prior covariance matrix
sigma_df_post <- tt + sigma_df_prior + r # Posterior degrees of freedom (Koop et al., 2010, eq. 8)

# Initial values
beta <- object[["initial"]][["beta"]]
sigma_i <- object[["initial"]][["u_sigma_inv"]]

g_i <- sigma_i

iterations <- object[["model"]][["iterations"]] # Number of iterations of the Gibbs sampler
burnin <- object[["model"]][["burnin"]] # Number of burn-in draws
draws <- iterations + burnin # Total number of draws

# Data containers
draws_alpha <- matrix(NA, k_alpha, iterations)
draws_beta <- matrix(NA, k_beta, iterations)
draws_pi <- matrix(NA, k * k_w, iterations)
draws_gamma <- matrix(NA, k_gamma, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
  
  # Draw conditional mean parameters
  temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = sigma_i,
                         v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
                         gamma_mu_prior = a_mu_prior,
                         gamma_v_i_prior = a_v_i_prior)
  alpha <- temp[["alpha"]]
  beta <- temp[["beta"]]
  Pi <- temp[["Pi"]]
  gamma <- temp[["Gamma"]]
  
  # Draw variance-covariance matrix
  u <- y - Pi %*% w - matrix(gamma, k) %*% x
  sigma_scale_post <- solve(tcrossprod(u) + v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
  sigma_i <- matrix(rWishart(1, sigma_df_post, sigma_scale_post)[,, 1], k)
  sigma <- solve(sigma_i)
  
  # Update g_i
  g_i <- sigma_i
  
  # Store draws
  if (draw > burnin) {
    draws_alpha[, draw - burnin] <- alpha
    draws_beta[, draw - burnin] <- beta
    draws_pi[, draw - burnin] <- Pi
    draws_gamma[, draw - burnin] <- gamma
    draws_sigma[, draw - burnin] <- sigma
  }
}
```

The `bvec` function can be used to collect output of the above Gibbs
sampler in a standardized object of class `bvecmodel`, which can be used
further for forecasting, impulse response analysis or forecast error
variance decomposition. The error correction term is described by the
draws of `alpha` and `beta` and not by those of `Pi`, since the latter
cannot be decomposed into the former.

``` r

# Number of non-deterministic coefficients
k_nondet <- (k_x - 1) * k

# Generate bvecmodel object
bvec_est <- bvec(y = object[["data"]][["train"]][["y"]],
                 w = object[["data"]][["train"]][["w"]],
                 x = object[["data"]][["train"]][["x"]][, 1:6],
                 x_d = object[["data"]][["train"]][["x"]][, -(1:6), drop = FALSE],
                 r = 1,
                 alpha = draws_alpha,
                 beta = draws_beta,
                 Gamma = draws_gamma[1:k_nondet,],
                 C = draws_gamma[(k_nondet + 1):nrow(draws_gamma),],
                 Sigma = draws_sigma)
```

The draws were obtained from the rescaled error correction term, so they
refer to the rescaled series in `w`. Function `rescale_error_correction`
undoes the transformation of `scale_error_correction`: it puts the
series in `object$data$train$w` back on the scale of the input data and
multiplies the draws of `beta` by the inverse of the scaling factors.
The draws of `alpha` are not affected, because the error correction term
of the estimated model is $`\alpha \beta^{\prime} D^{-1} w_t`$, where
$`D`$ is the diagonal matrix of scaling factors.

``` r

bvec_est <- rescale_error_correction(bvec_est)
```

This step is required before the model is used any further. Functions
such as `vec_to_var` recover the levels of the endogenous variables from
the differences in `y` and their lags in `w`, which is only possible if
both are on the same scale.

Obtain point estimates of cointegration variables:

``` r

beta <- bvec_est[["posterior"]][["beta"]][["coeffs"]] # Obtain draws of beta
beta <- apply(beta / beta[, 1], 2, median) # Normalise and obtain medians
beta <- matrix(beta, k_w) # Transform vector into a matrix
beta <- round(beta, 3) # Round values
dimnames(beta) <- list(dimnames(w)[[1]], NULL) # Rename matrix dimensions

beta # Print
#>       [,1]
#> l.lr  1.00
#> l.Dp -3.73
```

Posterior draws can be inspected with `plot`. Note that the coefficients
of the error correction term are displayed as draws of `Pi`, i.e. as the
product of `alpha` and `beta`, since the two matrices are only
identified up to a rotation.

``` r

plot(bvec_est)
```

![plot of chunk histograms](figures/bvec-histograms-1.png)

plot of chunk histograms

![plot of chunk histograms](figures/bvec-histograms-2.png)

plot of chunk histograms

![plot of chunk histograms](figures/bvec-histograms-3.png)

plot of chunk histograms

![plot of chunk histograms](figures/bvec-histograms-4.png)

plot of chunk histograms

Obtain summaries of posterior draws

``` r

summary(bvec_est)
#> 
#> Bayesian VEC model with p = 4 
#> 
#> Endogenous variables: lr, Dp
#> 
#> Variable: lr 
#> 
#>                Mean       SD  Naive SD Time-series SD      2.5%        50%
#> l.lr      0.0001826 0.005978 8.454e-05      8.454e-05 -0.012610  0.0002075
#> l.Dp     -0.0016969 0.020973 2.966e-04      2.966e-04 -0.042542 -0.0014739
#> d.lr.l01  0.4577010 0.083269 1.178e-03      1.178e-03  0.294385  0.4578165
#> d.Dp.l01  0.0143547 0.021850 3.090e-04      3.171e-04 -0.028898  0.0141437
#> d.lr.l02 -0.2179819 0.089660 1.268e-03      1.268e-03 -0.395582 -0.2178561
#> d.Dp.l02  0.0068847 0.021430 3.031e-04      3.031e-04 -0.036199  0.0070069
#> d.lr.l03  0.1535708 0.082841 1.172e-03      1.172e-03 -0.009668  0.1529289
#> d.Dp.l03 -0.0092369 0.017375 2.457e-04      2.457e-04 -0.042910 -0.0093650
#> const    -0.0073756 0.009003 1.273e-04      1.308e-04 -0.024748 -0.0075247
#>             97.5%  
#> l.lr      0.01220  
#> l.Dp      0.03909  
#> d.lr.l01  0.62765 *
#> d.Dp.l01  0.05786  
#> d.lr.l02 -0.04196 *
#> d.Dp.l02  0.04837  
#> d.lr.l03  0.31754  
#> d.Dp.l03  0.02397  
#> const     0.01083  
#> 
#> Variable: Dp 
#> 
#>              Mean      SD  Naive SD Time-series SD      2.5%      50%    97.5%
#> l.lr      0.11775 0.05243 0.0007414       0.001363  0.012807  0.11811  0.22234
#> l.Dp     -0.44188 0.10205 0.0014432       0.001511 -0.639957 -0.44277 -0.23932
#> d.lr.l01  0.43813 0.38683 0.0054705       0.005803 -0.336262  0.44031  1.17628
#> d.Dp.l01 -0.28792 0.10423 0.0014741       0.001598 -0.491242 -0.28728 -0.08439
#> d.lr.l02  0.14555 0.41643 0.0058893       0.005889 -0.666667  0.14524  0.96768
#> d.Dp.l02 -0.12033 0.09902 0.0014003       0.001522 -0.316389 -0.11916  0.07470
#> d.lr.l03  0.49342 0.38754 0.0054807       0.005606 -0.295404  0.50143  1.23670
#> d.Dp.l03  0.01055 0.07979 0.0011285       0.001186 -0.145929  0.01018  0.16865
#> const     0.12823 0.06019 0.0008512       0.001835  0.009656  0.12933  0.24556
#>           
#> l.lr     *
#> l.Dp     *
#> d.lr.l01  
#> d.Dp.l01 *
#> d.lr.l02  
#> d.Dp.l02  
#> d.lr.l03  
#> d.Dp.l03  
#> const    *
#> 
#> Variance-covariance matrix:
#> 
#>           Mean        SD  Naive SD Time-series SD     2.5%      50%    97.5%  
#> lr_lr 0.005023 0.0005729 8.101e-06      8.596e-06 0.004016 0.004983 0.006258 *
#> lr_Dp 0.004788 0.0019967 2.824e-05      2.911e-05 0.001005 0.004755 0.008958 *
#> Dp_Dp 0.110004 0.0129794 1.836e-04      1.904e-04 0.087269 0.109124 0.137207 *
```

### Default algorithm

As an alternative to a user-specific algorithm, function
`add_posterior_coefficients` can be used to estimate such a model as
well:

``` r

bvec_est <- add_posterior_coefficients(object)
bvec_est <- rescale_error_correction(bvec_est)
```

### A prior centred on the maximum likelihood estimate

The prior used above is uniform on the cointegration space. That is not
a consequence of `p_tau_i = 1` alone: the sampler uses `v_i` and
`p_tau_i` only through their product, so with `v_i = 0` the prior on the
space is uniform whatever `p_tau_i` is. An informative prior on the
space needs a positive `v_i`.

A natural centre for such a prior is the space spanned by Johansen’s
(1995) maximum likelihood estimate of $`\beta`$. With
`coint = list(p_tau_i = "ml")`, `add_priors` builds $`P_\tau`$ so that,
given $`\alpha`$, the prior on how far $`\beta`$ tilts away from that
space is the sampling distribution of the estimator itself, with its
precision multiplied by `weight`. With `weight = 1` the prior adds as
much information about the space as the sample holds. With `v_i = "ml"`
the shrinkage of the loadings is set from their maximum likelihood
estimate as well, so that it is weak on the scale the data put them on.

Because the prior is built from Johansen’s estimate on the series as
they are stored in the model, it moves with the units of the error
correction term, and scaling the series makes practically no difference
to the posterior. The model below is therefore estimated without
`scale_error_correction`, and consequently without
`rescale_error_correction`. If the series are scaled, that has to happen
before `add_priors`, which is enforced.

``` r

object_ml <- create_bvecmodel(data, p = 4, r = 1,
                              const = "unrestricted",
                              iterations = 5000, burnin = 1000)

object_ml <- add_priors(object_ml,
                        coef = list(v_i = 0, v_i_det = 0),
                        coint = list(v_i = "ml", p_tau_i = "ml", weight = 1),
                        sigma = list(df = 1, scale = .0001))

set.seed(7654321)
object_ml <- add_initial_values(object_ml, method = "maxlik")

bvec_ml <- add_posterior_coefficients(object_ml)
```

The posterior median of the normalised cointegration vector stays close
to that under the uniform prior, while its credible interval narrows:

``` r

beta_ratio <- function(object) {
  draws <- object[["posterior"]][["beta"]][["coeffs"]]
  quantile(draws[, 2] / draws[, 1], c(0.05, 0.25, 0.5, 0.75, 0.95))
}

round(rbind("Uniform prior" = beta_ratio(bvec_est),
            "ML prior, weight = 1" = beta_ratio(bvec_ml)), 3)
#>                          5%    25%   50%    75%    95%
#> Uniform prior        -8.886 -4.798 -3.73 -3.033 -2.337
#> ML prior, weight = 1 -6.538 -4.437 -3.73 -3.213 -2.653
```

Since the centre is estimated from the same sample that the model is
estimated on, the data enter twice, and the posterior overstates the
precision of the cointegration space. Smaller values of `weight` reduce
this. The prior is most useful where the likelihood alone is weakly
informative about the space, for example in short samples or in larger
systems.

## Evaluation

Posterior draws can be thinned with function `thin`:

``` r

bvec_est <- thin(bvec_est, thin = 5)
```

The function `vec_to_var` can be used to transform the VEC model into a
VAR in levels:

``` r

bvar_form <- vec_to_var(bvec_est)
```

The output of `vec_to_var` is an object of class `bvarmodel` for which
summary statistics, forecasts, impulse responses and variance
decompositions can be obtained in the usual manner using `summary`,
`predict`, `irf` and `fevd`, respectively.

``` r

summary(bvar_form)
#> 
#> Bayesian VAR model with p = 4 
#> 
#> Endogenous variables: lr, Dp
#> 
#> Variable: lr 
#> 
#>            Mean       SD  Naive SD Time-series SD     2.5%       50%    97.5%  
#> lr.l1  1.462600 0.080821 0.0025558      0.0025558  1.30481  1.463401  1.61341 *
#> Dp.l1  0.012825 0.017809 0.0005632      0.0005632 -0.02134  0.012691  0.04681  
#> lr.l2 -0.678378 0.140281 0.0044361      0.0044361 -0.95878 -0.677750 -0.41486 *
#> Dp.l2 -0.007967 0.018494 0.0005848      0.0005848 -0.04631 -0.008298  0.02779  
#> lr.l3  0.367034 0.141654 0.0044795      0.0046471  0.09629  0.368067  0.63775 *
#> Dp.l3 -0.016042 0.017226 0.0005447      0.0005447 -0.05128 -0.015946  0.01561  
#> lr.l4 -0.150961 0.085027 0.0026888      0.0026888 -0.31248 -0.154417  0.01747  
#> Dp.l4  0.009305 0.017704 0.0005598      0.0005598 -0.02400  0.009041  0.04341  
#> const -0.007044 0.009282 0.0002935      0.0002935 -0.02570 -0.007319  0.01294  
#> 
#> Variable: Dp 
#> 
#>            Mean      SD Naive SD Time-series SD     2.5%      50%  97.5%  
#> lr.l1  0.552796 0.38736 0.012249       0.012249 -0.23189  0.56319 1.2829  
#> Dp.l1  0.270377 0.08403 0.002657       0.002657  0.10763  0.27196 0.4436 *
#> lr.l2 -0.268056 0.67276 0.021275       0.022296 -1.51349 -0.25655 1.0647  
#> Dp.l2  0.166112 0.08262 0.002613       0.002613  0.01237  0.16384 0.3259 *
#> lr.l3  0.314862 0.66809 0.021127       0.023902 -1.03624  0.31717 1.5361  
#> Dp.l3  0.130134 0.08475 0.002680       0.002680 -0.04382  0.13125 0.2913  
#> lr.l4 -0.481880 0.38308 0.012114       0.012114 -1.18983 -0.48205 0.2782  
#> Dp.l4 -0.006325 0.07951 0.002514       0.002514 -0.16467 -0.00548 0.1422  
#> const  0.126567 0.06054 0.001914       0.002289  0.00541  0.12863 0.2378 *
#> 
#> Variance-covariance matrix:
#> 
#>           Mean        SD  Naive SD Time-series SD     2.5%      50%    97.5%  
#> lr_lr 0.005053 0.0005778 1.827e-05      2.076e-05 0.004044 0.005020 0.006301 *
#> lr_Dp 0.004830 0.0019254 6.089e-05      6.089e-05 0.001248 0.004791 0.008741 *
#> Dp_Dp 0.109786 0.0127505 4.032e-04      4.032e-04 0.086723 0.109157 0.137814 *
```

### Forecasts

``` r

bvar_form <- add_forecast_input(bvar_form, n_ahead = 10)
```

``` r

bvar_form <- add_posterior_forecasts(bvar_form)
```

``` r

# Genrate forecasts
bvar_pred <- predict(bvar_form, n_ahead = 10)

# Plot forecasts
plot(bvar_pred, n_pre = 20)
```

![plot of chunk forecasts](figures/bvec-forecasts-1.png)

plot of chunk forecasts

![plot of chunk forecasts](figures/bvec-forecasts-2.png)

plot of chunk forecasts

### Forecast error impulse response

``` r

FEIR <- irf(bvar_form, impulse = "lr", response = "Dp", n_ahead = 20)

plot(FEIR, main = "Forecast Error Impulse Response", xlab = "Period", ylab = "Response")
```

![plot of chunk feir](figures/bvec-feir-1.png)

plot of chunk feir

## Forecast error variance decomposition

``` r

bvar_fevd_oir <- fevd(bvar_form, response = "Dp", n_ahead = 20)

plot(bvar_fevd_oir, main = "OIR-based FEVD of inflation")
```

![plot of chunk fevd-oir](figures/bvec-fevd-oir-1.png)

plot of chunk fevd-oir

## Citing bvartools

If you use `bvartools` in published work, please cite it.
`citation("bvartools")` prints the reference, and the package has the
DOI [10.5281/zenodo.22736604](https://doi.org/10.5281/zenodo.22736604),
which always resolves to the latest archived version.

## References

Johansen, S. (1995). *Likelihood-based inference in cointegrated vector
autoregressive models*. Oxford: Oxford University Press.

Koop, G., León-González, R., & Strachan R. W. (2010). Efficient
posterior simulation for cointegrated models with priors on the
cointegration space. *Econometric Reviews, 29*(2), 224-242.
<https://doi.org/10.1080/07474930903382208>

Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis
in linear multivariate models, *Economics Letters, 58*, 17-29.
<https://doi.org/10.1016/S0165-1765(97)00214-0>
