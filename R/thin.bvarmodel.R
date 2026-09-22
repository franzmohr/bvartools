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

  kept <- .thin_positions(x, thin)
  x[["posterior"]] <- .thin_draws(x[["posterior"]], kept[["positions"]], kept[["draws"]], thin)

  return(x)
}


# Which rows of the pooled draws thinning keeps. Every chain is thinned on its
# own: they are stacked one after the other, and thinning the stack as one
# sequence kept draws 4, 8, ... of the first chain and 3, 7, ... of the next
# whenever a chain's length was not a multiple of 'thin', and left chains of
# unequal length otherwise. A discounted model has no draws at all.
.thin_positions <- function(x, thin) {

  if (.is_discount(x)) {
    stop("A discounted model has a closed-form posterior rather than draws, so there is ",
         "nothing to thin.", call. = FALSE)
  }

  draws <- nrow(x[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  if (is.null(draws)) {
    stop("Argument 'x' has no posterior draws to thin.", call. = FALSE)
  }

  chains <- x[["model"]][["chains"]]
  chains <- if (is.null(chains)) 1L else as.integer(chains)
  if (draws %% chains != 0) {
    stop("The ", draws, " draws do not divide into ", chains, " chains of equal length.",
         call. = FALSE)
  }
  n <- draws %/% chains
  .check_thin(thin, n)

  within <- seq(from = thin, to = n, by = thin)
  positions <- as.vector(outer(within, (seq_len(chains) - 1) * n, "+"))

  list("positions" = positions, "draws" = draws)
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
#
# The labels are counted from the ones the draws already carry. A sampler that
# kept one draw in `before` labels them `before`, `2 before`, ..., and keeping
# every `thin`-th of those keeps one draw in `before * thin` of the chain. Starting
# the labels afresh at one would say the draws came from the start of the chain.
.thin_draws <- function(posterior, pos_thin, draws, thin) {

  for (i in names(posterior)) {
    element <- posterior[[i]]
    if (is.null(element)) {
      next
    }
    if (is.list(element) && !inherits(element, "mcmc")) {
      posterior[[i]] <- .thin_draws(element, pos_thin, draws, thin)
    } else if (NROW(element) == draws) {
      mcpar <- attr(element, "mcpar")
      if (is.null(mcpar)) {
        mcpar <- c(1, draws, 1)
      }
      # The end is counted from the number of draws kept rather than read off
      # the last position: kept per chain, the positions are not evenly spaced
      # across the pooled rows, and coda wants labels that are.
      start <- mcpar[1] + (pos_thin[1] - 1) * mcpar[3]
      posterior[[i]] <- coda::mcmc(.draws_matrix(element)[pos_thin, , drop = FALSE],
                                   start = start,
                                   end = start + (length(pos_thin) - 1) * thin * mcpar[3],
                                   thin = thin * mcpar[3])
    }
  }

  return(posterior)
}


# The thinning interval a model's sampler keeps one draw in, validated. NULL for
# one, which is how a model that keeps every draw says so.
.check_sampler_thin <- function(thin) {
  if (!is.numeric(thin) || length(thin) != 1 || is.na(thin) ||
      thin < 1 || thin != round(thin)) {
    stop("Argument 'thin' must be a single positive integer.")
  }
  if (thin == 1) {
    return(NULL)
  }
  return(as.integer(thin))
}


# The sampler's draws of one block as an mcmc object whose labels say which draws
# were kept. BayesTS keeps the last of every `thin` draws after the burn-in, so the
# kept ones are iterations thin, 2 thin, ... after the burn-in -- the labels
# bayests writes into a model file. Without thinning this is exactly what
# coda::as.mcmc() returns, so an unthinned model is unchanged down to the type of
# its labels.
.mcmc_draws <- function(model, draws) {
  thin <- model[["thin"]]
  if (is.null(thin) || thin == 1) {
    return(coda::as.mcmc(draws))
  }
  thin <- as.numeric(thin)
  return(coda::mcmc(draws, start = thin, end = NROW(draws) * thin, thin = thin))
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