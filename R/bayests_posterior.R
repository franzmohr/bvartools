#' Posterior Simulation with the BayesTS Executable
#'
#' Returns a function that simulates the posterior of a model with the standalone
#' BayesTS executable rather than with the samplers compiled into this package.
#'
#' @param executable the path to the BayesTS executable. If \code{NULL}
#' (default), option \code{bvartools.bayests_executable} or, if that is not set,
#' environment variable \code{BAYESTS_EXECUTABLE}.
#' @param library_path a character vector of directories the executable needs on
#' its search path for its runtime libraries, put in front of \code{PATH} while
#' it runs. A packaged BayesTS carries its libraries and does not need it; an
#' executable from a build tree may.
#' @param scratch the directory where the model files are written while BayesTS
#' draws into them. They are removed afterwards.
#'
#' @details
#' The result is meant for argument \code{posterior_function} of
#' \code{\link{add_posterior_coefficients}}:
#' \preformatted{model <- add_posterior_coefficients(model,
#'   posterior_function = bayests_posterior())}
#' It writes the model with \code{\link{write_to_hdf5}}, runs
#' \code{bayests posterior} on the file without its log-likelihood and forecast
#' steps, reads the draws back with \code{\link{read_model_from_hdf5}} and adds
#' them to the model as element \code{posterior}. Nothing else about the model is
#' changed. The draws have the elements, dimensions and thinning of those of the
#' internal samplers, whose C++ code BayesTS shares; a BayesTS build linked
#' against an optimised BLAS draws faster.
#'
#' BayesTS draws with the seed in \code{object$model$seed}, which
#' \code{\link{add_initial_values}} sets and \code{\link{add_seed}} replaces. A
#' model without one is given a seed drawn from R's random number generator,
#' since BayesTS would otherwise start every model from the same state of its
#' own generator. The seed used is kept in the returned model.
#'
#' The executable runs single threaded unless \code{OMP_NUM_THREADS} says
#' otherwise, which suits simulating one model per worker of a cluster.
#'
#' @return A function of one argument, a model of class 'bvarmodel' or
#' 'bvecmodel', that returns that model with its posterior draws.
#'
#' @family posterior simulation
#' @export
bayests_posterior <- function(executable = NULL, library_path = NULL,
                              scratch = tempdir()) {

  if (is.null(executable)) {
    executable <- getOption("bvartools.bayests_executable",
                            Sys.getenv("BAYESTS_EXECUTABLE"))
  }
  if (length(executable) != 1 || is.na(executable) || !nzchar(executable)) {
    stop("No BayesTS executable given. Pass 'executable', or set option ",
         "'bvartools.bayests_executable' or environment variable ",
         "'BAYESTS_EXECUTABLE'.", call. = FALSE)
  }
  if (!file.exists(executable)) {
    stop("BayesTS executable not found: ", executable, call. = FALSE)
  }
  executable <- normalizePath(executable, winslash = "/")
  force(library_path)
  force(scratch)

  quote_type <- if (.Platform$OS.type == "windows") "cmd" else "sh"

  function(object) {

    if (is.null(object[["model"]][["seed"]])) {
      object[["model"]][["seed"]] <- .draw_model_seed()
    }

    dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
    file <- tempfile(pattern = "bayests_", tmpdir = scratch, fileext = ".h5")
    log <- sub("[.]h5$", ".log", file)
    on.exit(unlink(c(file, log)), add = TRUE)

    if (length(library_path) > 0) {
      old_path <- Sys.getenv("PATH")
      on.exit(Sys.setenv(PATH = old_path), add = TRUE)
      Sys.setenv(PATH = paste(c(library_path, old_path), collapse = .Platform$path.sep))
    }

    write_to_hdf5(object, filename = file)
    status <- system2(executable,
                      c("posterior", shQuote(file, type = quote_type),
                        "--no-loglik", "--no-forecasts"),
                      stdout = log, stderr = log)
    if (!identical(as.integer(status), 0L)) {
      reason <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
      reason <- reason[seq_len(length(reason)) > length(reason) - 5]
      stop("BayesTS failed with exit status ", status, ": ",
           paste(reason, collapse = " | "), call. = FALSE)
    }

    draws <- read_model_from_hdf5(filename = file)[["posterior"]]
    if (is.null(draws)) {
      stop("BayesTS wrote no draws for the model.", call. = FALSE)
    }
    object[["posterior"]] <- draws
    object
  }
}
