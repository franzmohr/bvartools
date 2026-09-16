#' Seed of the Posterior Simulation
#'
#' Sets the seed with which the posterior of a model is simulated.
#'
#' @param object a model, usually after \code{\link{add_initial_values}}: an
#' object of class 'bvarmodel' or 'bvecmodel', or a list of them of class
#' 'modellist' or 'expandingwindow'.
#' @param seed a non-negative whole number no larger than
#' \code{.Machine$integer.max}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' Calling this function is optional. \code{\link{add_initial_values}} already
#' stores a seed, drawn from R's random number generator, in element \code{seed}
#' of \code{object$model}. \code{add_seed} replaces it, for example to give a
#' model a seed that depends neither on the state of R's generator nor on the
#' worker of a cluster that happens to simulate it.
#'
#' The seed is part of the model. \code{\link{write_to_hdf5}} writes it as
#' attribute \code{seed} of group \code{/model}, where the BayesTS executable
#' reads it, and \code{\link{read_model_from_hdf5}} reads it back. The internal
#' samplers of \code{\link{add_posterior_coefficients}} draw with it too, so a
#' model with a given seed gives the same draws however R's generator stands.
#' BayesTS and this package use different generators, so the same seed gives
#' different draws in the two.
#'
#' A list of models gets the seeds \code{seed}, \code{seed + 1}, \ldots, one per
#' model in the order of its elements, counting through nested lists, so that no
#' two of its models draw the same numbers. Elements that are not estimated, such
#' as external forecasts, are left as they are and are not counted.
#'
#' @return \code{object} with the seed set.
#'
#' @examples
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' # Generate model, add priors and initial values
#' model <- create_bvarmodel(data = e1, p = 2, iterations = 100, burnin = 50)
#' model <- add_priors(model)
#' model <- add_initial_values(model)
#'
#' # add_initial_values() has set a seed; replace it
#' model <- add_seed(model, 20260916)
#' model[["model"]][["seed"]]
#'
#' @family posterior simulation
#' @export
add_seed <- function(object, seed, ...) {
  UseMethod("add_seed")
}

#' @rdname add_seed
#' @export
add_seed.bvarmodel <- function(object, seed, ...) {
  object[["model"]][["seed"]] <- .check_seed(seed)
  object
}

#' @rdname add_seed
#' @export
add_seed.bvecmodel <- function(object, seed, ...) {
  object[["model"]][["seed"]] <- .check_seed(seed)
  object
}

#' @rdname add_seed
#' @export
add_seed.modellist <- function(object, seed, ...) {
  .add_seed_to_list(object, seed, ...)
}

#' @rdname add_seed
#' @export
add_seed.expandingwindow <- function(object, seed, ...) {
  .add_seed_to_list(object, seed, ...)
}

# One seed per model, counted through nested lists: use_expanding_window() on a
# 'modellist' returns a 'modellist' of 'expandingwindow' lists, and numbering
# each inner list from 'seed' again would give models of different lists the
# same seed.
.add_seed_to_list <- function(object, seed, ...) {
  seed <- .check_seed(seed)
  n <- 0
  set_seeds <- function(x) {
    if (inherits(x, c("modellist", "expandingwindow"))) {
      for (i in seq_along(x)) {
        x[[i]] <- set_seeds(x[[i]])
      }
    } else if (inherits(x, c("bvarmodel", "bvecmodel"))) {
      x <- add_seed(x, .offset_seed(seed, n), ...)
      n <<- n + 1
    }
    x
  }
  set_seeds(object)
}

# A seed is stored as an integer. A plain 20260916 in R is a double, and the
# BayesTS reader takes a whole number of either kind, but an integer is what
# write_to_hdf5() should hand it.
.check_seed <- function(seed) {
  if (length(seed) != 1 || !is.numeric(seed) || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) {
    stop("Argument 'seed' must be a single whole number between 0 and ",
         .Machine$integer.max, ".", call. = FALSE)
  }
  as.integer(seed)
}

# 'seed + offset', wrapped so that it stays a valid seed.
.offset_seed <- function(seed, offset) {
  as.integer((as.numeric(seed) + offset) %% (as.numeric(.Machine$integer.max) + 1))
}

# The seed add_initial_values() and bayests_posterior() give a model that has
# none, from R's generator as it stands.
.draw_model_seed <- function() {
  sample.int(.Machine$integer.max, 1L)
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
