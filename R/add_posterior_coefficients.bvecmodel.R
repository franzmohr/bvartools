#' Posterior Simulation of Model Coefficients
#' 
#' Forwards model input to posterior simulation functions for vector error correction models.
#' 
#' @param object an object of class 'bvecmodel', usually, a result of a
#' call to \code{\link{create_bvecmodel}} in combination with
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
#' @details The function implements commonly used posterior simulation algorithms for Bayesian VAR models with
#' both constant and time varying parameters (TVP) as well as stochastic volatility. It can produce posterior
#' draws for standard BVAR models with independent normal-Wishart priors, which can be augmented by stochastic
#' search variable selection (SSVS) as proposed by Geroge et al. (2008) or Bayesian variable selection (BVS)
#' as proposed in Korobilis (2013). Both SSVS or BVS can also be applied to the covariances of the error term.
#' 
#' The implementation follows the descriptions in Chan et al. (2019), George et al. (2008) and Korobilis (2013).
#' For all approaches the SUR form of a VAR model is used to obtain posterior draws. The algorithm is implemented
#' in C++ to reduce calculation time.
#' 
#' The function also supports structural BVEC models, where the structural coefficients are estimated from
#' contemporary endogenous variables, which corresponds to the so-called (A-model). Currently, only
#' specifications are supported, where the structural matrix contains ones on its diagonal and all lower
#' triangular elements are freely estimated. Since posterior draws are obtained based on the SUR form of
#' the VEC model, the structural coefficients are drawn jointly with the other coefficients.
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
#' 'gvecmodel' -- that ends the run on the first specification that fails,
#' rather than leaving that one without a posterior and carrying it into
#' whatever reads the results.
#'
#' @return The object in \code{object} with the element \code{posterior} added, whose
#' elements hold the draws after burn-in as \code{\link[coda]{mcmc}} objects with one
#' row per draw and one column per parameter, in element \code{coeffs}: \code{beta}, the
#' cointegration coefficients, \eqn{K_\beta r} columns or \eqn{T K_\beta r} for time
#' varying cointegration, \code{a}, the loadings and the remaining coefficients, and the
#' draws of the error term as described in \code{\link{add_posterior_coefficients.bvarmodel}}.
#'
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' George, E. I., Sun, D., & Ni, S. (2008). Bayesian stochastic search for VAR model
#' restrictions. \emph{Journal of Econometrics, 142}(1), 553--580.
#' \doi{10.1016/j.jeconom.2007.08.017}
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2010). Efficient posterior
#' simulation for cointegrated models with priors on the cointegration space.
#' \emph{Econometric Reviews, 29}(2), 224--242.
#' \doi{10.1080/07474930903382208}
#' 
#' Korobilis, D. (2013). VAR forecasting using Bayesian variable selection.
#' \emph{Journal of Applied Econometrics, 28}(2), 204--230. \doi{10.1002/jae.1271}
#' 
#' @examples
#' 
#' # Load data 
#' data("e6")
#' e6 <- e6 * 100
#' 
#' # Generate model
#' model <- create_bvecmodel(e6, p = 1, r = 1, const = "restricted",
#'                           iterations = 10, burnin = 10)
#' # Chosen number of iterations and burn-in should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     coint = list(v_i = 0, p_tau_i = 1),
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
add_posterior_coefficients.bvecmodel <- function(object, posterior_function = NULL, chains = NULL, ...){

  # Several chains are the same simulation run once per chain, each with a seed
  # of its own, and pooled: see chain_diagnostics().
  chains <- .check_chains(object, chains)
  if (chains > 1) {
    return(.run_chains(object, chains, function(single) {
      add_posterior_coefficients.bvecmodel(single, posterior_function = posterior_function, chains = 1, ...)
    }))
  }
  
  class_of_object <- class(object)
  
  # A failing simulation raises its error rather than being turned into an
  # object that looks estimated. See the note in
  # add_posterior_coefficients.bvarmodel(), which this mirrors: the 'error'
  # marker this used to return in place of a posterior was read by nothing, so
  # the failure surfaced later and somewhere else.
  if (is.null(posterior_function)) {

    # Not a sampler. The filter runs in the vendored core like every algorithm
    # below, but what it returns is a posterior rather than a chain, so it does
    # not fall through to the mcpar the draws of the others are labelled with.
    if (.is_discount(object)) {
      return(.discount_coefficients(object))
    }

    # Check if the input is suitable for the posterior simulation functions
    .check_bvecpost_input(object)

    algorithm <- object[["model"]][["algorithm"]]

    if (algorithm %in% c("VecKlgs2010", "VecNormalGamma", "VecNormalWishart",
                         "VecNormalStochvol", "VecTvpGamma", "VecTvpWishart",
                         "VecTvpStochvol")) {
      object <- .with_model_seed(object[["model"]][["seed"]], switch(algorithm,
                       VecKlgs2010 = .VecKlgs2010Coefficients(object),
                       VecNormalGamma = .VecNormalGammaCoefficients(object),
                       VecNormalStochvol = .VecNormalStochvolCoefficients(object),
                       VecNormalWishart = .VecNormalWishartCoefficients(object),
                       VecTvpGamma = .VecTvpGammaCoefficients(object),
                       VecTvpStochvol = .VecTvpStochvolCoefficients(object),
                       VecTvpWishart = .VecTvpWishartCoefficients(object)))
    } else {
      stop("Algorithm '", algorithm, "' not supported.")
    }
    object <- .raise_core_warnings(object)

    for (i in c("a", "beta", "psi", "u_sigma_inv", "u_omega_inv")) {
      if (!is.null(object[["posterior"]][[i]][["coeffs"]])) {
        object[["posterior"]][[i]][["coeffs"]] <- .mcmc_draws(object[["model"]], object[["posterior"]][[i]][["coeffs"]])
      }
      # Only the cointegration block has one of these, and only when the
      # prior made rho a parameter rather than a hyperparameter. NULL
      # everywhere else, which is the same as not having it.
      if (!is.null(object[["posterior"]][[i]][["rho"]])) {
        object[["posterior"]][[i]][["rho"]] <- .mcmc_draws(object[["model"]], object[["posterior"]][[i]][["rho"]])
      }
      if (!is.null(object[["posterior"]][[i]][["lambda"]])) {
        object[["posterior"]][[i]][["lambda"]] <- .mcmc_draws(object[["model"]], object[["posterior"]][[i]][["lambda"]])
      }
      if (!is.null(object[["posterior"]][[i]][["sigma"]])) {
        object[["posterior"]][[i]][["sigma"]] <- .mcmc_draws(object[["model"]], object[["posterior"]][[i]][["sigma"]])
      }
      # What a block drawn under the non-centred prior 'omega_v' adds beside
      # 'sigma'. Chains like the others, and the HDF5 writer reads their mcpar.
      for (j in c("omega", "omega_log_zero", "omega_log_zero_joint")) {
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
