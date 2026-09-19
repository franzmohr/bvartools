#' @include irf.R
NULL

#' Impulse Responses and Variance Decompositions of a Stored Model
#'
#' Impulse response functions and forecast error variance decompositions of a
#' model in an HDF5 file, computed a piece of the chain at a time.
#'
#' @param x an object of class 'bvarfile', from \code{\link{open_model}}.
#' @param chunk how many draws are read at a time. Defaults to 100.
#' @param ... arguments of \code{\link{irf.bvarmodel}} or
#' \code{\link{fevd.bvarmodel}}, which do the work on each piece.
#'
#' @details
#' Both quantities are sums over the draws: an impulse response is a set of
#' quantiles of the responses of the draws, a variance decomposition their mean.
#' Neither needs the draws together, so the chain is read in pieces of
#' \code{chunk}, each piece contributes what it has, and the pieces are put
#' together at the end. What comes out is what the same call on the model in
#' memory gives.
#'
#' For a variance decomposition the pieces are combined before the shares are
#' normalised and before groups are collapsed, so \code{normalise_gir} and
#' \code{max_groups} describe the whole chain rather than the piece they were
#' applied to.
#'
#' @return What \code{\link{irf.bvarmodel}} and \code{\link{fevd.bvarmodel}}
#' return: an object of class 'bvarirf' or 'bvarfevd'.
#'
#' @examples
#'
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' model <- gen_var(e1, p = 2, deterministic = "const",
#'                  iterations = 100, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     sigma = list(df = 1, scale = .0001))
#' model <- add_posterior_coefficients(add_initial_values(model))
#'
#' file <- file.path(tempdir(), "bvartools-example-irf.h5")
#' unlink(file)
#' write_to_hdf5(model, filename = file)
#' stored <- open_model(file)
#'
#' ir <- irf(stored, impulse = "invest", response = "income", n_ahead = 5,
#'           chunk = 25)
#' shares <- fevd(stored, response = "income", n_ahead = 5, chunk = 25)
#'
#' @name analysis_of_stored_models
#' @family model comparison
NULL

#' @rdname analysis_of_stored_models
#' @export
#' @method irf bvarfile
irf.bvarfile <- function(x, ..., chunk = 100) {

  arguments <- list(...)
  ci <- if (is.null(arguments[["ci"]])) .95 else arguments[["ci"]]
  keep_draws <- isTRUE(arguments[["keep_draws"]])
  arguments[["ci"]] <- NULL
  arguments[["keep_draws"]] <- NULL

  # Every piece hands back its own draws of the responses, which is what the
  # quantiles below are taken over. They are small -- one row per draw and one
  # column per horizon -- so the chain of responses is held whole where the
  # chain of the model is not.
  pieces <- map_draws(x, function(model) {
    do.call(irf, c(list(model), arguments, list(keep_draws = TRUE)))
  }, chunk = chunk)

  result <- do.call("rbind", lapply(pieces, unclass))

  if (!keep_draws) {
    result <- .summarise_irf_draws(result, ci)
  }

  class(result) <- append("bvarirf", class(result))

  result
}

#' @rdname analysis_of_stored_models
#' @export
#' @method fevd bvarfile
fevd.bvarfile <- function(x, ..., chunk = 100) {

  arguments <- list(...)
  normalise_gir <- isTRUE(arguments[["normalise_gir"]])
  max_groups <- arguments[["max_groups"]]
  arguments[["normalise_gir"]] <- NULL
  arguments[["max_groups"]] <- NULL

  # A decomposition is the mean of the decompositions of the draws, so a piece
  # contributes its own mean and its weight. Normalising or collapsing groups
  # first would average quantities that are no longer shares of the same thing,
  # so both are left to the end.
  pieces <- map_draws(x, function(model) {
    shares <- do.call(fevd, c(list(model), arguments,
                              list(normalise_gir = FALSE, max_groups = NULL)))
    list(shares = unclass(shares), draws = nrow(model[["posterior"]][["u_sigma_inv"]][["coeffs"]]))
  }, chunk = chunk)

  weights <- vapply(pieces, function(piece) piece[["draws"]], numeric(1))
  result <- Reduce(`+`, Map(function(piece, weight) piece[["shares"]] * weight,
                            pieces, weights)) / sum(weights)

  type <- if (is.null(arguments[["type"]])) "oir" else arguments[["type"]]
  if (type %in% c("gir", "sgir") && normalise_gir) {
    result <- t(apply(result, 1, function(z) z / sum(z)))
  }
  colnames(result) <- x[["model"]][["endogen"]]

  max_groups <- .check_max_groups(max_groups)
  if (!is.null(max_groups) && max_groups < x[["model"]][["k"]]) {
    result <- .limit_fevd_groups(result, max_groups)
  }

  result <- stats::ts(result, start = 0, frequency = 1)
  class(result) <- append("bvarfevd", class(result))

  result
}
