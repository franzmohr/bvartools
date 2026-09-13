#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class 'bvarmodel'.
#' 
#' @param x an object of class 'bvarmodel'.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
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
#' # Thinning
#' model <- thin(model, 2)
#' 
#' @return An object of class 'bvarmodel'.
#' 
#' @export
#' @method thin bvarmodel
thin.bvarmodel <- function(x, thin = 10, ...) {

  draws <- nrow(x[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  .check_thin(thin, draws)
  pos_thin <- seq(from = thin, to = draws, by = thin)
  x[["posterior"]] <- .thin_draws(x[["posterior"]], pos_thin, draws, thin)

  return(x)
}


# Thins every element of a posterior that holds one row per draw, however
# deeply it is nested and whatever it is called.
#
# The methods used to name the elements to thin, and the names fell behind the
# samplers: the state variances of a time varying model, a$sigma and
# psi$sigma, and the draws of rho, beta$rho, kept every draw while the rest of
# the posterior was thinned, so that row i no longer belonged to the same draw
# across blocks. Anything with as many rows as there are draws is a block of
# draws, since the samplers store one row per draw and nothing else.
.thin_draws <- function(posterior, pos_thin, draws, thin) {

  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]

  for (i in names(posterior)) {
    element <- posterior[[i]]
    if (is.null(element)) {
      next
    }
    if (is.list(element) && !inherits(element, "mcmc")) {
      posterior[[i]] <- .thin_draws(element, pos_thin, draws, thin)
    } else if (NROW(element) == draws) {
      posterior[[i]] <- coda::mcmc(.draws_matrix(element)[pos_thin, , drop = FALSE],
                                   start = start, end = end, thin = thin)
    }
  }

  return(posterior)
}


# The draws of an mcmc object as a plain matrix, one row per draw, and with the
# dimnames it had. as.matrix() would do the first but add an empty dimnames list
# to draws that had none, so that thinning by one or cutting a window to the
# whole sample no longer returned what it was given.
.draws_matrix <- function(draws) {
  attr(draws, "mcpar") <- NULL
  class(draws) <- NULL
  if (is.null(dim(draws))) {
    draws <- matrix(draws, ncol = 1)
  }
  return(draws)
}