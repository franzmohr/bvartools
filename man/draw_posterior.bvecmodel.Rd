% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_posterior.bvecmodel.R
\name{draw_posterior.bvecmodel}
\alias{draw_posterior.bvecmodel}
\title{Posterior Simulation for Vector Error Correction Models}
\usage{
\method{draw_posterior}{bvecmodel}(object, FUN = NULL, mc.cores = NULL, ...)
}
\arguments{
\item{object}{a list of model specifications, which should be passed on
to function \code{FUN}. Usually, the output of a call to \code{\link{gen_vec}}
in combination with \code{\link{add_priors}}.}

\item{FUN}{the function to be applied to each list element in argument \code{object}.
If \code{NULL} (default), the internal function \code{\link{bvecpost}} is used.}

\item{mc.cores}{the number of cores to use, i.e. at most how many child
processes will be run simultaneously. The option is initialized from
environment variable MC_CORES if set. Must be at least one, and
parallelization requires at least two cores.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
For multiple models a list of objects of class \code{bvarlist}.
For a single model the object has the class of the output of the applied posterior
simulation function. In case the package's own functions are used, this will
be \code{"bvec"}.
}
\description{
Forwards model input to posterior simulation functions for vector error correction models.
}
\examples{

# Load data 
data("e6")
e6 <- e6 * 100

# Generate model
model <- gen_vec(e6, p = 1, r = 1, const = "restricted",
                 iterations = 10, burnin = 10)
# Chosen number of iterations and burn-in should be much higher.

# Add priors
model <- add_priors(model)

# Obtain posterior draws
object <- draw_posterior(model)

}
\references{
Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
simulation for cointegrated models with priors on the cointegration space.
\emph{Econometric Reviews, 29}(2), 224--242.
\doi{10.1080/07474930903382208}

Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference in
a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
\doi{10.1016/j.jeconom.2011.07.007}
}
