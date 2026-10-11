#' Impulse Responses and Variance Decompositions of a Stored Model
#'
#' Impulse response functions and forecast error variance decompositions of a
#' model in an HDF5 file, computed a piece of the chain at a time.
#'
#' @param x an object of class 'modelfile', from \code{\link[bayests]{read_model_file}}.
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
#' For a variance decomposition the pieces hand back the shares of each of
#' their draws, and the decomposition is summarised once over all of them, so
#' \code{statistic}, \code{ci}, \code{normalise_gir} and \code{max_groups}
#' describe the whole chain rather than the piece they were applied to.
#'
#' An identification given per draw is read at the draw's position in the
#' whole chain: a list of impact matrices is cut to each piece's draws, and an
#' impact function is called with the draw's number in the chain. A piece in
#' which no draw is identified -- under \code{type = "sign"}, or an impact
#' function that declines all of its draws -- contributes nothing.
#'
#' @return What \code{\link{irf.bvarmodel}} and \code{\link{fevd.bvarmodel}}
#' return: an object of class 'bvarirf' or 'bvarfevd'.
#'
#' @examples
#'
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
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
#' stored <- read_model_file(file)
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
#' @method irf modelfile
irf.modelfile <- function(x, ..., chunk = 100) {

  arguments <- list(...)
  ci <- if (is.null(arguments[["ci"]])) .95 else arguments[["ci"]]
  keep_draws <- isTRUE(arguments[["keep_draws"]])
  arguments[["ci"]] <- NULL
  arguments[["keep_draws"]] <- NULL

  # Every piece hands back its own draws of the responses, which is what the
  # quantiles below are taken over. They are small -- one row per draw and one
  # column per horizon -- so the chain of responses is held whole where the
  # chain of the model is not.
  pieces <- .map_draw_pieces(x, function(model, piece) {
    if (!.has_identified_draws(model, arguments[["type"]])) {
      return(NULL)
    }
    piece_arguments <- .piece_impact(arguments, piece)
    tryCatch(do.call(irf, c(list(model), piece_arguments, list(keep_draws = TRUE))),
             bvartools_no_draws = function(e) NULL)
  }, chunk = chunk)
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0) {
    stop("No draw of the stored model is identified.")
  }

  result <- do.call("rbind", lapply(pieces, unclass))

  if (!keep_draws) {
    result <- .summarise_irf_draws(result, ci)
  }

  class(result) <- append("bvarirf", class(result))

  result
}

#' @rdname analysis_of_stored_models
#' @export
#' @method fevd modelfile
fevd.modelfile <- function(x, ..., chunk = 100) {

  arguments <- list(...)
  normalise_gir <- isTRUE(arguments[["normalise_gir"]])
  type <- if (is.null(arguments[["type"]])) "oir" else arguments[["type"]]
  statistic <- if (is.null(arguments[["statistic"]])) "mean" else arguments[["statistic"]]
  ci <- arguments[["ci"]]
  max_groups <- .check_max_groups(arguments[["max_groups"]])
  if (!statistic %in% c("mean", "median")) {
    stop("Argument 'statistic' must be \"mean\" or \"median\".")
  }
  .check_ci(ci)
  if (!is.null(ci) && !is.null(max_groups)) {
    stop("Arguments 'ci' and 'max_groups' cannot be combined: the bounds of a pooled column ",
         "are not the pooled bounds.")
  }
  for (i in c("normalise_gir", "max_groups", "statistic", "ci")) {
    arguments[[i]] <- NULL
  }

  # Every piece hands back the shares of each of its draws, and the
  # decomposition is summarised once over all of them, as fevd() does for a
  # model in memory. That is exact for every statistic and needs no weights:
  # averaging the pieces' own results was exact for the mean alone, and only
  # where every draw of a piece counted -- not under an identification that
  # keeps some draws and not others.
  pieces <- .map_draw_pieces(x, function(model, piece) {
    if (!.has_identified_draws(model, type)) {
      return(NULL)
    }
    piece_arguments <- .piece_impact(arguments, piece)
    tryCatch(do.call(fevd, c(list(model), piece_arguments,
                             list(normalise_gir = FALSE, max_groups = NULL,
                                  .return_draws = TRUE))),
             bvartools_no_draws = function(e) NULL)
  }, chunk = chunk)
  pieces <- Filter(function(piece) !is.null(piece) && ncol(piece) > 0, pieces)
  if (length(pieces) == 0) {
    stop("No draw of the stored model is identified.")
  }

  draws <- do.call(cbind, pieces)
  varnames <- x[["model"]][["endogen"]]
  n_ahead <- nrow(draws) / length(varnames) - 1
  .summarise_fevd_draws(draws, n_ahead, varnames, type, normalise_gir, statistic, ci,
                        max_groups)
}


# map_draws() for the analysis above, handing each piece the positions in the
# chain of the draws it holds as well. An identification given per draw -- a
# list of impact matrices, or a function of the draw's number -- has to be
# read at those positions: every piece used to be handed the whole list, which
# it refused for its length, and a function saw the piece's own numbering.
.map_draw_pieces <- function(x, f, chunk) {
  if (!inherits(x, "modelfile")) {
    stop("Argument 'x' must be of class 'modelfile'. Use read_model_file().")
  }
  chunk <- max(1L, as.integer(chunk))
  pieces <- split(seq_len(x[["draws"]]), ceiling(seq_len(x[["draws"]]) / chunk))
  lapply(pieces, function(piece) {
    model <- read_model_from_hdf5(filename = x[["filename"]], group = x[["group"]],
                                  draws = piece)
    f(model, piece)
  })
}

# The arguments of a piece, with a per-draw identification cut to its draws.
.piece_impact <- function(arguments, piece) {
  impact <- arguments[["impact"]]
  if (is.function(impact)) {
    arguments[["impact"]] <- function(draw, i) impact(draw, piece[i])
  } else if (is.list(impact) && !is.data.frame(impact)) {
    arguments[["impact"]] <- impact[piece]
  }
  arguments
}

# Stops for an identification that declined every draw, with a class of its
# own, so that a piece of a stored chain in which it did can be skipped rather
# than stopping the analysis of the others.
.stop_no_identified_draws <- function() {
  stop(structure(class = c("bvartools_no_draws", "error", "condition"),
                 list(message = paste0("No draw has an impact matrix: the identification ",
                                       "declined every one."), call = NULL)))
}

# Whether a piece has a draw to analyse. Under a sign identification a piece
# can hold none that was identified, and analysing it stopped the whole call.
.has_identified_draws <- function(model, type) {
  if (!identical(type, "sign")) {
    return(TRUE)
  }
  q <- model[["posterior"]][["q"]][["coeffs"]]
  !is.null(q) && any(!is.na(as.matrix(q)[, 1]))
}
