
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
#' generator stands. BayesTS and this package use different generators, so the
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

# Evaluates 'expr' with R's generator set to 'seed', with R's default kinds so
# that a model draws the same on a cluster worker running L'Ecuyer-CMRG as in a
# plain session, and puts the generator back afterwards, kinds included. Without
# a seed 'expr' is evaluated with the generator as it stands. 'expr' is a
# promise and is only forced after set.seed().
.with_model_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  seed <- .check_seed(seed)

  global <- globalenv()
  # Before RNGkind(), which initialises the generator when it has no state yet.
  had_state <- exists(".Random.seed", envir = global, inherits = FALSE)
  if (had_state) {
    old_state <- get(".Random.seed", envir = global, inherits = FALSE)
  }
  old_kind <- RNGkind()
  on.exit({
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (had_state) {
      assign(".Random.seed", old_state, envir = global)
    } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
      rm(".Random.seed", envir = global)
    }
  }, add = TRUE)

  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  expr
}
