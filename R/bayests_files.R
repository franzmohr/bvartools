#' Run BayesTS on Stored Models
#'
#' Returns a function that runs the standalone BayesTS executable on model files
#' or on directories of them, several at a time, leaving the results in the
#' files.
#'
#' @param executable the path to the BayesTS executable. If \code{NULL}
#' (default), option \code{bvartools.bayests_executable} or, if that is not set,
#' environment variable \code{BAYESTS_EXECUTABLE}.
#' @param library_path a character vector of directories the executable needs on
#' its search path for its runtime libraries, put in front of \code{PATH} while
#' it runs. A packaged BayesTS carries its libraries and does not need it; an
#' executable from a build tree may.
#'
#' @details
#' \code{\link{bayests_posterior}} draws one model that is in the session: it
#' writes the model to a scratch file, runs BayesTS on it and reads the draws
#' back. For models that are already stored -- a folder written with
#' \code{\link{write_to_hdf5}}, which is how a model too large to hold is worked
#' on -- that is three copies of the draws per model and a session holding them,
#' for no purpose. The function returned here runs BayesTS on the stored files
#' themselves and leaves what it produces where it wrote it, so the draws never
#' pass through R.
#'
#' The returned function takes paths, which are model files or directories of
#' them, the name of a BayesTS command, and how many to run at once:
#' \describe{
#'   \item{\code{"posterior"}}{draws the posterior, and with it the pointwise
#'   log-likelihood unless \code{--no-loglik} is passed in \code{args}.}
#'   \item{\code{"loglik"}}{the log-likelihood of a model that has draws.}
#'   \item{\code{"forecasts"}}{the forecasts of a model that has draws.}
#'   \item{\code{"check"}}{reads and validates the models without running them,
#'   which is what to call before an estimation that will take hours.}
#' }
#' A model is drawn with the seed in its file, so it does not matter whether it
#' is run on its own, as one of a directory, or through
#' \code{\link{bayests_posterior}}: the draws are the same.
#'
#' @section Running several at once:
#' The executable works through a directory one model at a time and gains
#' nothing from \code{OMP_NUM_THREADS} on models of the size a country model of
#' a global VAR has, so the way to use a machine is to give it several paths and
#' a \code{jobs} above one.
#'
#' They are run as processes of the operating system rather than on a cluster of
#' R workers, because an idle R worker is not cheap: on Windows a fresh one
#' commits about 2 GB before it does anything and about 4 GB once it has loaded
#' a package of this kind, which for six workers is more memory than the six
#' BayesTS processes they would be waiting for. Here one R session starts the
#' processes, each writing its output to a log of its own, and waits for them.
#'
#' @return A function of paths, a command, further arguments for the executable
#' and the number of processes to run at once. It returns the paths invisibly
#' and fails with the end of BayesTS's output if any of the runs does.
#'
#' @examples
#'
#' \dontrun{
#' run <- bayests_files(executable = "/opt/bayests/bin/bayests")
#'
#' # One directory of models, drawn in place
#' run("models/gvar/submodels/US")
#'
#' # Every country, six processes at a time
#' run(list.dirs("models/gvar/submodels", recursive = FALSE), jobs = 6)
#'
#' # Validate first, draw without the log-likelihood
#' run("models/gvar/submodels/US", command = "check")
#' run("models/gvar/submodels/US", args = "--no-loglik")
#' }
#'
#' @family posterior simulation
#' @export
bayests_files <- function(executable = NULL, library_path = NULL) {

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

  function(path, command = "posterior", args = character(0), jobs = 1,
           poll = 1) {

    if (length(path) == 0) {
      return(invisible(path))
    }
    if (anyNA(path) || !all(nzchar(path))) {
      stop("Argument 'path' must be files or directories.", call. = FALSE)
    }
    missing <- path[!file.exists(path)]
    if (length(missing) > 0) {
      stop("No such file or directory: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    full <- normalizePath(path, winslash = "/")
    jobs <- max(1L, as.integer(jobs))

    jobs_running <- list()
    queue <- seq_along(full)
    failures <- character(0)

    collect <- function(job) {
      status <- tryCatch(as.integer(trimws(readLines(job[["done"]], warn = FALSE)[1])),
                         warning = function(w) NA_integer_)
      unlink(job[["script"]])
      if (identical(status, 0L)) {
        unlink(c(job[["log"]], job[["done"]]))
        return(NULL)
      }
      reason <- if (file.exists(job[["log"]])) readLines(job[["log"]], warn = FALSE) else character()
      reason <- reason[seq_len(length(reason)) > length(reason) - 5]
      unlink(c(job[["log"]], job[["done"]]))
      paste0(job[["path"]], " (exit status ", status, "): ",
             paste(reason, collapse = " | "))
    }

    on.exit({
      for (job in jobs_running) {
        unlink(c(job[["script"]], job[["log"]], job[["done"]]))
      }
    }, add = TRUE)

    while (length(queue) > 0 || length(jobs_running) > 0) {

      while (length(jobs_running) < jobs && length(queue) > 0) {
        index <- queue[1]
        queue <- queue[-1]
        jobs_running[[length(jobs_running) + 1]] <-
          .bayests_start(executable, library_path, command, full[index], args)
      }

      Sys.sleep(poll)

      finished <- vapply(jobs_running, function(job) file.exists(job[["done"]]),
                         logical(1))
      for (job in jobs_running[finished]) {
        failures <- c(failures, collect(job))
      }
      jobs_running <- jobs_running[!finished]
    }

    if (length(failures) > 0) {
      stop("BayesTS failed on ", length(failures),
           if (length(failures) == 1) " path: " else " paths: ",
           paste(failures, collapse = " || "), call. = FALSE)
    }

    invisible(path)
  }
}

# One BayesTS process, started and not waited for. The command is put in a
# script of its own, which is what makes the exit status readable: a process
# started with wait = FALSE hands nothing back, so the script writes its status
# into a file the caller polls for. The script also carries the search path the
# executable needs, rather than the session changing its own PATH while several
# processes are in flight.
.bayests_start <- function(executable, library_path, command, path, args) {

  stem <- tempfile(pattern = "bayests_")
  log <- paste0(stem, ".log")
  done <- paste0(stem, ".done")

  windows <- .Platform$OS.type == "windows"
  script <- paste0(stem, if (windows) ".bat" else ".sh")

  if (windows) {
    back <- function(x) gsub("/", "\\\\", x, fixed = TRUE)
    lines <- "@echo off"
    if (length(library_path) > 0) {
      lines <- c(lines, paste0("set \"PATH=", paste(back(library_path), collapse = ";"),
                               ";%PATH%\""))
    }
    lines <- c(lines,
               paste(c(paste0("\"", back(executable), "\""), command,
                       paste0("\"", back(path), "\""), args,
                       paste0("> \"", back(log), "\" 2>&1")), collapse = " "),
               # In parentheses, because "echo 0> file" is read as a
               # redirection of stream 0 rather than as an echo of "0".
               paste0("(echo %ERRORLEVEL%)> \"", back(done), "\""))
    writeLines(lines, script)
    system2("cmd", c("/c", shQuote(script, type = "cmd")), wait = FALSE,
            stdout = NULL, stderr = NULL)
  } else {
    lines <- "#!/bin/sh"
    if (length(library_path) > 0) {
      lines <- c(lines, paste0("PATH=\"", paste(library_path, collapse = ":"),
                               ":$PATH\"; export PATH"))
    }
    lines <- c(lines,
               paste(c(shQuote(executable), command, shQuote(path), args,
                       ">", shQuote(log), "2>&1"), collapse = " "),
               paste("echo $? >", shQuote(done)))
    writeLines(lines, script)
    Sys.chmod(script, "0755")
    system2(script, wait = FALSE, stdout = NULL, stderr = NULL)
  }

  list(path = path, script = script, log = log, done = done)
}
