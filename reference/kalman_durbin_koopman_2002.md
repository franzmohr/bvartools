# Durbin and Koopman Simulation Smoother

An implementation of the Kalman filter and backward smoothing algorithm
proposed by Durbin and Koopman (2002).

## Usage

``` r
kalman_durbin_koopman_2002(y, z, sigma_u, sigma_v, B, a_init, P_init)
```

## Arguments

- y:

  a \\K \times T\\ matrix of endogenous variables.

- z:

  a \\KT \times M\\ matrix of explanatory variables.

- sigma_u:

  the constant \\K \times K\\ error variance-covariance matrix. For time
  varying variance-covariance matrices a \\KT \times K\\ can be
  specified.

- sigma_v:

  the constant \\M \times M\\ coefficient variance-covariance matrix.
  For time varying variance-covariance matrices a \\MT \times M\\ can be
  specified.

- B:

  the constant \\M \times M\\ autocorrelation matrix of the transition
  equation. For a time varying transition an \\MT \times M\\ matrix can
  be specified.

- a_init:

  an M-dimensional vector of initial states.

- P_init:

  an \\M \times M\\ variance-covariance matrix of the initial states.

## Value

A \\M \times T+1\\ matrix of state vector draws. Column \\i\\ is the
state the observation in column \\i\\ of `y` loads on, for \\i =
1,...,T\\, so a caller wanting one state per observation takes the first
\\T\\ columns. The first column is not a state preceding the sample: it
is \\a_1\\, already conditioned on every observation. The last column is
the transition applied once past the end of the sample and is informed
by no observation; it is returned because the recursions build it on the
way.

## Details

The function uses algorithm 2 from Durbin and Koopman (2002) to produce
a draw of the state vector \\a_t\\ for \\t = 1,...,T\\ for a state space
model with measurement equation \$\$y_t = Z_t a_t + u_t\$\$ and
transition equation \$\$a\_{t + 1} = B_t a\_{t} + v_t,\$\$ where \\u_t
\sim N(0, \Sigma\_{u,t})\\ and \\v_t \sim N(0, \Sigma\_{v,t})\\. \\y_t\\
is a K-dimensional vector of endogenous variables and \\Z_t =
z_t^{\prime} \otimes I_K\\ is a \\K \times M\\ matrix of regressors with
\\z_t\\ as a vector of regressors.

The algorithm takes into account Jarociński (2015), where a possible
missunderstanding in the implementation of the algorithm of Durbin and
Koopman (2002) is pointed out. Following that note the function sets the
mean of the initial state to zero in the first step of the algorithm.

This is the routine the time varying parameter samplers of the package
use. It comes from the vendored BayesTS core, so there is one
implementation of the smoother rather than one for R and one for the
samplers.

The draw uses R's random number generator, so
[`set.seed`](https://rdrr.io/r/base/Random.html) makes it reproducible.

## References

Durbin, J., & Koopman, S. J. (2002). A simple and efficient simulation
smoother for state space time series analysis. *Biometrika, 89*(3),
603–615.

Jarociński, M. (2015). A note on implementing the Durbin and Koopman
simulation smoother. *Computational Statistics and Data Analysis, 91*,
1–3.
[doi:10.1016/j.csda.2015.05.001](https://doi.org/10.1016/j.csda.2015.05.001)

## Examples

``` r

# Load data
data("e1")
data <- diff(log(e1))

# Generate model data
temp <- create_bvarmodel(data = data, p = 2, deterministic = "const",
                         iterations = 1, burnin = 0)
y <- t(temp$data$train$y)
z <- temp$data$train$z
k <- nrow(y)
tt <- ncol(y)
m <- ncol(z)

# Priors
a_mu_prior <- matrix(0, m)
a_v_i_prior <- diag(0.1, m)

a_Q <- diag(.0001, m)

# Initial value of Sigma
sigma <- tcrossprod(y) / tt
sigma_i <- solve(sigma)

# Initial values for Kalman filter
y_init <- y * 0
a_filter <- matrix(0, m, tt + 1)

# Initialise the Kalman filter
for (i in 1:tt) {
  y_init[, i] <- y[, i] - z[(i - 1) * k + 1:k,] %*% a_filter[, i]
}
a_init <- post_normal_sur(y = y_init, z = z, sigma_i = sigma_i,
                          a_prior = a_mu_prior, v_i_prior = a_v_i_prior)
y_filter <- matrix(y) - z %*% a_init
y_filter <- matrix(y_filter, k) # Reshape

# Kalman filter and backward smoother
a_filter <- kalman_durbin_koopman_2002(y = y_filter, z = z, sigma_u = sigma,
                                       sigma_v = a_Q, B = diag(1, m),
                                       a_init = matrix(0, m), P_init = a_Q)

a <- a_filter + matrix(a_init, m, tt + 1)
```
