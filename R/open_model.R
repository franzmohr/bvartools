#' Open a Model Stored in an HDF5 File
#'
#' Opens a model written with \code{\link{write_to_hdf5}} and returns a handle
#' to it, which the analysis functions work through a piece of the chain at a
#' time rather than reading every draw into the session.
#'
#' @param filename path to an HDF5 file holding a model.
#' @param group the group the model's tree hangs under inside its file, as
#' \code{\link{read_model_from_hdf5}} takes it. Defaults to \code{""}, the root.
#'
#' @details
#' The draws of a model are what it is large in. A model with time varying
#' coefficients keeps a path per draw, and a global model solved from many
#' sub-models keeps a square matrix per lag and draw, so a chain that a
#' posterior summary is comfortable with is a chain an analysis cannot hold
#' beside everything else it needs.
#'
#' A handle carries what a model is -- its specification, its data and its
#' priors -- and the length of its chain, but none of the draws. The methods for
#' it read the draws in pieces and keep only what they are asked for: the
#' responses of an impulse response, the shares of a variance decomposition.
#' What they return is what the same call on the model in memory returns.
#'
#' @return A list of class 'bvarfile' with the elements \code{filename},
#' \code{group}, \code{model}, \code{data}, \code{priors} and \code{draws}, the
#' number of draws in the file.
#'
#' @examples
#'
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' model <- add_priors(model,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     sigma = list(df = 1, scale = .0001))
#' model <- add_posterior_coefficients(add_initial_values(model))
#'
#' file <- file.path(tempdir(), "bvartools-example-model.h5")
#' unlink(file)
#' write_to_hdf5(model, filename = file)
#'
#' stored <- open_model(file)
#' stored
#'
#' # The analysis reads the draws in pieces
#' ir <- irf(stored, impulse = "invest", response = "income", n_ahead = 5,
#'           chunk = 25)
#'
#' @family model comparison
#' @export
open_model <- function(filename, group = "") {

  if (!file.exists(filename)) {
    stop("File ", filename, " does not exist.")
  }

  model <- read_model_from_hdf5(filename = filename, group = group,
                                draws = integer(0))
  draws <- .chain_length_in_hdf5(filename, group)

  result <- list(filename = filename, group = group,
                 model = model[["model"]], data = model[["data"]],
                 priors = model[["priors"]], draws = draws)
  class(result) <- c("bvarfile", "list")

  result
}

#' @rdname open_model
#' @param x an object of class 'bvarfile'.
#' @param ... further arguments passed to or from other methods.
#' @export
print.bvarfile <- function(x, ...) {

  cat("Model in", x[["filename"]], "\n")
  cat(x[["draws"]], "draws of a", x[["model"]][["algorithm"]], "model of",
      x[["model"]][["k"]], "variables\n")

  invisible(x)
}

# How many draws a stored model holds, read from the shape of its error
# precision, which every posterior has.
.chain_length_in_hdf5 <- function(filename, group = "") {

  group <- .normalize_hdf5_group(group)
  file <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (file$is_valid) file$close_all(), add = TRUE)
  root <- .hdf5_model_root(file, group)

  if (!"posterior" %in% names(root)) {
    stop("The model in ", filename, " holds no posterior draws.")
  }
  posterior <- root[["posterior"]]
  if (!"u_sigma_inv" %in% names(posterior)) {
    stop("The model in ", filename, " holds no draws of its error precision, ",
         "so the length of its chain cannot be read.")
  }
  draws <- posterior[["u_sigma_inv"]][["coeffs"]]$dims[1]
  file$close_all()

  as.integer(draws)
}

#' Apply a Function to the Draws of a Stored Model
#'
#' Reads the chain of a stored model in pieces and applies a function to the
#' model holding each piece, which is how the analysis of a model too large to
#' hold is done.
#'
#' @param x an object of class 'bvarfile', from \code{\link{open_model}}.
#' @param f a function taking a model and returning whatever the caller wants to
#' keep of that piece of the chain.
#' @param ... further arguments passed to \code{f}.
#' @param chunk how many draws are read at a time. Defaults to 100.
#'
#' @details
#' The model \code{f} is given is the model in the file with the draws of one
#' piece, so anything that works on a model works on it. What it returns is
#' collected in a list, one element per piece, and combining those is the
#' caller's business: the responses of an impulse response are stacked, the
#' shares of a variance decomposition are averaged over the pieces they came
#' from.
#'
#' @return A list with one element per piece of the chain.
#'
#' @examples
#'
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
#' model <- add_priors(model,
#'                     coef = list(v_i = 0, v_i_det = 0),
#'                     sigma = list(df = 1, scale = .0001))
#' model <- add_posterior_coefficients(add_initial_values(model))
#'
#' file <- file.path(tempdir(), "bvartools-example-draws.h5")
#' unlink(file)
#' write_to_hdf5(model, filename = file)
#'
#' stored <- open_model(file)
#'
#' # The mean of every coefficient, without reading the chain at once
#' sums <- map_draws(stored, function(model) {
#'   colSums(unclass(model[["posterior"]][["a"]][["coeffs"]]))
#' }, chunk = 25)
#' means <- Reduce(`+`, sums) / stored[["draws"]]
#'
#' @family model comparison
#' @export
map_draws <- function(x, f, ..., chunk = 100) {

  if (!inherits(x, "bvarfile")) {
    stop("Argument 'x' must be of class 'bvarfile'. Use open_model().")
  }
  f <- match.fun(f)
  chunk <- max(1L, as.integer(chunk))

  pieces <- split(seq_len(x[["draws"]]), ceiling(seq_len(x[["draws"]]) / chunk))

  lapply(pieces, function(piece) {
    model <- read_model_from_hdf5(filename = x[["filename"]], group = x[["group"]],
                                  draws = piece)
    f(model, ...)
  })
}
