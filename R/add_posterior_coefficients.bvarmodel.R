#' Posterior Simulation of Model Coefficients
#' 
#' Forwards model input to posterior simulation functions for vector autoregressive models.
#' 
#' @param object an object of class 'bvarmodel', usually, a result of a
#' call to \code{\link{create_bvarmodel}} in combination with
#' \code{\link{add_priors}} and \code{\link{add_initial_values}}.
#' @param posterior_function the function to be applied to the model in argument \code{object}.
#' If \code{NULL} (default), internal functions are used.
#' @param chains the number of chains to simulate. If \code{NULL} (default), the value in
#' \code{object$model$chains} is used, and one chain when there is none. Each chain is the
#' same simulation with a seed of its own -- the first with the seed of the model, so that one
#' chain draws what it always has -- and the chains are pooled, one after the other, in the
#' draws of \code{posterior}, so that every later step uses all of them. Their number, when
#' above one, is stored in \code{object$model$chains}, and \code{\link{chain_diagnostics}} compares them. A
#' \code{posterior_function} is called once per chain. Not available for discounted models,
#' whose posterior is not a chain.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details Unless \code{posterior_function} is specified, the function forwards
#' the model input to the package's own posterior functions.
#' 
#' The internal samplers draw with the seed in \code{object$model$seed}, which
#' \code{\link{add_initial_values}} sets and \code{\link{add_seed}} replaces.
#' R's random number generator is set to that seed, with R's default kinds, for
#' the simulation and put back as it was afterwards. A call of \code{set.seed()}
#' between \code{add_initial_values()} and this function therefore does not
#' change the draws. A model without a seed draws from R's generator as it
#' stands. A \code{posterior_function} is called as it is and decides itself
#' what to do with the seed; see \code{\link{bayests_posterior}}.
#'
#' A sampler that cannot run raises its error rather than returning something.
#' The message names what about the input it could not work with. Applied to a
#' list of models -- a 'modellist', an 'expandingwindow' or, in \pkg{bgvars}, a
#' 'gvarmodel' -- that ends the run on the first specification that fails,
#' rather than leaving that one without a posterior and carrying it into
#' whatever reads the results.
#' 
#' @return The object in \code{object} with the element \code{posterior} added. Each of its
#' elements is a list whose element \code{coeffs} holds the draws after burn-in as a
#' \code{\link[coda]{mcmc}} object with one row per draw and one column per parameter:
#' \describe{
#'   \item{\code{a}}{the coefficients, \eqn{M} columns or \eqn{TM} for TVP models, in the
#'   order of the columns of \code{data$train$z} with the contemporaneous coefficients
#'   of structural models last. With variable selection element \code{lambda} holds the
#'   inclusion indicators, and for TVP models element \code{sigma} the state variances.}
#'   \item{\code{u_sigma_inv}}{the inverse error covariance matrix, \eqn{K^2} columns, or
#'   \eqn{TK^2} if the error variances vary over time.}
#'   \item{\code{u_omega_inv}}{for all errors but \code{"wishart"}, the error precisions,
#'   \eqn{K} columns or \eqn{TK}.}
#'   \item{\code{psi}}{for \code{"gamma+covar"} and \code{"sv+covar"}, the error covariance
#'   coefficients.}
#'   \item{\code{u_scale}}{for \code{error = "ald"}, the \eqn{K} scales of the asymmetric
#'   Laplace distribution.}
#' }
#' Elements that do not apply to a model are absent or \code{NULL}. Note that
#' \code{\link{bvar}} expects draws in the transposed orientation.
#'
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' model <- add_posterior_coefficients(model)
#' 
#' @seealso \code{\link{bvartools_model}} describes the object this returns, element by element.
#' @family posterior simulation
#' @export
add_posterior_coefficients.bvarmodel <- function(object, posterior_function = NULL, chains = NULL, ...){

  # Several chains are the same simulation run once per chain, each with a seed
  # of its own, and pooled: see chain_diagnostics().
  chains <- .check_chains(object, chains)
  if (chains > 1) {
    return(.run_chains(object, chains, function(single) {
      add_posterior_coefficients.bvarmodel(single, posterior_function = posterior_function, chains = 1, ...)
    }))
  }
  
  # This  allows to employ the method with other compatible classes
  class_of_object <- class(object)
  
  # A failing simulation raises its error rather than being turned into an
  # object that looks estimated. What this used to do was catch it, print it and
  # return the unestimated model with an 'error' element bolted on -- a marker
  # nothing in this package or in bgvars ever read. The model that came back had
  # the class, the specification and the data of a fitted one and no posterior,
  # so the failure surfaced somewhere else entirely: as a missing element in a
  # later step, or as a model silently ranked against its siblings on a log
  # likelihood it did not have.
  #
  # It costs the batch methods something. add_posterior_coefficients() on a
  # 'modellist', an 'expandingwindow' or a 'gvarmodel' is an lapply over this,
  # so one unusable specification now stops the run instead of leaving a hole in
  # the results. Stopping on it is the lesser harm: a long batch that ends in an
  # error says which model was wrong and why, and one that quietly drops a model
  # does not.
  if (is.null(posterior_function)) {

    # Not a sampler. The filter runs in the vendored core like every algorithm
    # below, but what it returns is a posterior rather than a chain, so it does
    # not fall through to the mcpar the draws of the others are labelled with.
    if (.is_discount(object)) {
      return(.discount_coefficients(object))
    }

    # Check if the input is suitable for the posterior simulation functions
    .check_bvarpost_input(object)

    algorithm <- object[["model"]][["algorithm"]]

    if (algorithm %in% c("VarNormalAld", "VarNormalGamma", "VarNormalStochvol", "VarNormalWishart",
                         "VarTvpAld", "VarTvpGamma", "VarTvpStochvol", "VarTvpWishart")) {
      object <- .with_model_seed(object[["model"]][["seed"]], switch(algorithm,
                       VarNormalAld = .VarNormalAldCoefficients(object),
                       VarNormalGamma = .VarNormalGammaCoefficients(object),
                       VarNormalStochvol = .VarNormalStochvolCoefficients(object),
                       VarNormalWishart = .VarNormalWishartCoefficients(object),
                       VarTvpAld = .VarTvpAldCoefficients(object),
                       VarTvpGamma = .VarTvpGammaCoefficients(object),
                       VarTvpStochvol = .VarTvpStochvolCoefficients(object),
                       VarTvpWishart = .VarTvpWishartCoefficients(object)))
    } else {
      stop("Algorithm '", algorithm, "' not supported.")
    }
    object <- .raise_core_warnings(object)

    # 'omega' and the two 'omega_log_zero' elements are what a block drawn under
    # the non-centred prior 'omega_v' adds beside 'sigma'. They are chains like
    # the others, and the HDF5 writer reads their mcpar.
    for (i in c("a", "psi", "u_sigma_inv", "u_omega_inv", "u_scale")) {
      for (j in c("coeffs", "lambda", "sigma", "omega", "omega_log_zero", "omega_log_zero_joint")) {
        if (!is.null(object[["posterior"]][[i]][[j]])) {
          object[["posterior"]][[i]][[j]] <- .mcmc_draws(object[["model"]], object[["posterior"]][[i]][[j]])
        }
      }
    }

    # The C++ side names every block a model can have and leaves the ones this
    # model does not have as NULL. Such an element cannot be written to a file,
    # so a round trip lost it; it is dropped here, which reads the same.
    object[["posterior"]] <- .drop_null(object[["posterior"]])

  } else {
    # Apply own function
    object <- posterior_function(object)
  }
  
  class(object) <- class_of_object
  
  return(object)
}
