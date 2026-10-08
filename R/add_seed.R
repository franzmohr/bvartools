
#' Seed of the Posterior Simulation of a VAR or VEC Model
#'
#' The methods of \code{\link[bayests]{add_seed}} for VAR and VEC models: they
#' set the seed with which the posterior is simulated, in element \code{seed} of
#' \code{object$model}.
#'
#' @param object an object of class 'bvarmodel' or 'bvecmodel'.
#' @param seed a non-negative whole number no larger than
#' \code{.Machine$integer.max}.
#' @param ... not used.
#'
#' @details \code{\link[bayests]{write_to_hdf5}} writes the seed as attribute
#' \code{seed} of group \code{/model}, where the BayesTS executable reads it,
#' and the internal samplers of \code{\link[bayests]{add_posterior_coefficients}} draw
#' with it too, so a model with a given seed gives the same draws however R's
#' generator stands. So do the simulations that follow the coefficients:
#' \code{\link[bayests]{add_posterior_forecasts}}, scenarios included,
#' \code{\link[bayests]{add_predictive_loglik}}, \code{\link{add_sign_restrictions}}
#' and \code{\link{add_sign_zero_restrictions}}. Each of these three stages -- the
#' forecast, the predictive score and the rotations -- draws from a seed of its own
#' derived from the model's, as the BayesTS command line seeds its stages from
#' \code{/model/seed}, so that a forecast repeated gives the same draws and does
#' not reuse the random numbers the chain started from. A model without a seed
#' draws from R's generator as it stands throughout. BayesTS and this package use different generators, so the
#' same seed gives different draws in the two. A list of models is numbered by
#' the method of \pkg{bayests}.
#'
#' @return \code{object} with the seed set.
#'
#' @examples
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' # Generate model
#' model <- create_bvarmodel(data = e1, p = 2, iterations = 100, burnin = 50)
#'
#' # Set the seed of the posterior simulation
#' model <- add_seed(model, 20260916)
#' model[["model"]][["seed"]]
#'
#' @family posterior simulation
#' @export
add_seed.bvarmodel <- function(object, seed, ...) {
  object[["model"]][["seed"]] <- .check_seed(seed)
  object
}

#' @rdname add_seed.bvarmodel
#' @export
add_seed.bvecmodel <- function(object, seed, ...) {
  object[["model"]][["seed"]] <- .check_seed(seed)
  object
}

# .with_model_seed(), .stage_seed() and .with_stage_seed() come from bayests,
# which dfmtools takes them from too, so that the models of both packages are
# seeded alike.
