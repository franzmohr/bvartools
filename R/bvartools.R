#' bvartools: Bayesian Inference of Vector Autoregressive Models
#' 
#' A collection of R and C++ functions, which assist in the Bayesian inference of
#' vector autoregressive (VAR) and vector error correction (VEC) models.
#' 
#' The package \code{bvartools} implements some common functions used for Bayesian
#' inference for linear, multivariate time series models. It should give researchers
#' maximum freedom in setting up MCMC algorithms in R and keep calculation time
#' limited at the same time. This is achieved by implementing posterior simulation
#' functions in C++. Its main features are
#' \itemize{
#' \item The \code{bvar} and \code{bvec} functions collect the output of a Gibbs sampler
#' in standardised objects, which can be used for further analyses.
#' \item Further functions such as \code{predict}, \code{irf}, \code{fevd} for forecasting,
#' impulse response analysis and forecast error variance decomposition, respectively.
#' \item Computationally intensive functions - such as for posterior
#' simulation - are written in C++ using the \code{RcppArmadillo} package of Eddelbuettel
#' and Sanderson (2014).
#' \item Posterior simulation functions for commonly used Gibbs sampler algorithms.
#' }
#' 
#' @author Franz X. Mohr
#' @docType package
#' @name bvartools
#' 
#'
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias, J. L. (2019). \emph{Bayesian Econometric Methods}
#' (2nd ed.). Cambridge: University Press.
#' 
#' Durbin, J., & Koopman, S. J. (2002). A simple and efficient simulation smoother for
#' state space time series analysis. \emph{Biometrika, 89}(3), 603--615.
#' 
#' Eddelbuettel, D., & Sanderson C. (2014). RcppArmadillo: Accelerating R with high-performance
#' C++ linear algebra. \emph{Computational Statistics and Data Analysis, 71}, 1054--1063.
#' \doi{10.1016/j.csda.2013.02.005}
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Koop, G, & Korobilis, D. (2010), Bayesian multivariate time series Methods for empirical
#' macroeconomics, \emph{Foundations and Trends in Econometrics, 3}(4), 267--358.
#' \doi{10.1561/0800000013}
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
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230.
#' \doi{10.1002/jae.1271}
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Sanderson, C., & Curtin, R. (2016). Armadillo: a template-based C++ library for linear algebra.
#' \emph{Journal of Open Source Software, 1}(2), 26. \doi{10.21105/joss.00026}
#' 
#' @useDynLib bvartools, .registration = TRUE
#' @importFrom coda thin
#' @importFrom Rcpp sourceCpp
#' @import methods
#' @exportPattern "^[[:alpha:]]+"
NULL