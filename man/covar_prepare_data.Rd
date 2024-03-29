% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{covar_prepare_data}
\alias{covar_prepare_data}
\title{Covariance: Data Preparation}
\usage{
covar_prepare_data(y, omega_i, k, tt, tvp)
}
\arguments{
\item{y}{a \eqn{KT \times 1} vector of input data.}

\item{omega_i}{a \eqn{K \times K} or \eqn{KT \times KT} matrix of error variances.
The matrix must be sparse.}

\item{k}{an integer of the number of endogenous variables.}

\item{tt}{an integer of the number of observations.}

\item{tvp}{logical indicating if the SUR matrix with the values of regressors
should be prepared for the estimation of constant or time varying parameters.}
}
\value{
A list with three elements:
\item{y}{The prepared vector of endogenous variables.}
\item{z}{The prepared matrix of regressors.}
\item{omega_i}{The prepared diagonal matrix of measurement error variances.}
All matrices are returned as sparse matrices.
}
\description{
Convenience function, which generates the input data for posterior
simulation of covariance parameters.
}
\details{
For the model
\deqn{y_t = Z_{t} a_t + u_t}
with \eqn{u_t \sim N(0, \Psi \Omega_{t} \Psi^{\prime})} and \eqn{\Omega_{t}}
as a diagonal matrix of error variances, the function produces
the input data for the posterior simulation of the lower triangular covariance coefficients
of \eqn{\Psi} as presented in Primiceri (2005).
}
\examples{

# Create artificial data
k <- 3
tt <- 4
u <- matrix(1:(k * tt))
omega_i <- Matrix(diag(1:3, k))

# Generate input data (constant parameters)
covar_prepare_data(u, omega_i, k, tt, FALSE)

# Generate input data (time varying parameters)
covar_prepare_data(u, omega_i, k, tt, TRUE)

}
\references{
Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods} (2nd ed.). Cambridge: Cambridge University Press.

Primiceri, G. E. (2005). Time varying structural vector autoregressions and monetary policy. \emph{The Review of Economic Studies 72}(3), 821--852. \doi{10.1111/j.1467-937X.2005.00353.x}
}
