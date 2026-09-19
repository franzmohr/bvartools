#' Apply a Function to the Models of a Folder
#'
#' Reads every model of a folder, applies a function to it, writes it back and
#' drops it, so that a step over many models costs the memory of one model per
#' worker rather than of all of them.
#'
#' @param x an object of class 'bvarfolder', from \code{\link{open_models}}.
#' @param f a function taking a model -- a 'bvarmodel' or a 'bvecmodel' -- and
#' returning one. If it has an argument named \code{index}, the position of the
#' model in the manifest is passed to it, which is what numbers seeds and what
#' lets a caller tell the models apart.
#' @param ... further arguments passed to \code{f}.
#' @param models character vector of the models the step is applied to, as the
#' manifest names them, or \code{NULL}, the default, for all of them. The others
#' are left untouched, and so are their rows of the manifest.
#' @param write whether what \code{f} returns is written back. \code{TRUE},
#' the default, is a step of an estimation: \code{f} takes a model and returns
#' one, and the model in the file becomes what it returned. \code{FALSE} reads
#' only: \code{f} may return anything, nothing is written, and the results are
#' collected and returned.
#' @param cores the number of worker processes the models are handled on.
#' Defaults to one. With more than one the models are spread over a socket
#' cluster; every model has its own file, so the workers share nothing.
#' @param export character vector of names of objects in the calling
#' environment that the workers need, for instance the sampler that
#' \code{posterior_function} refers to. Ignored with one core.
#'
#' @details
#' The unit is one model, which is what a file of the folder holds. A model is
#' read with \code{\link{read_model_from_hdf5}}, passed to \code{f}, and written
#' back to its own file through a temporary file that replaces the old one once
#' it is complete, so an interrupted step leaves either the model as it was or
#' the model as \code{f} made it, never half of it. A file holding more than one
#' model cannot be worked on this way, since writing one of them back would have
#' to rewrite the others, and such a model is refused rather than risked.
#'
#' The manifest is rebuilt from the specifications of the models that come back,
#' so a step that changes what a model is -- its rank, its lag orders, the chain
#' it asks for -- leaves the table describing the files. It records what a model
#' is, not how many draws are left in it, so thinning the draws does not change
#' it.
#'
#' A step over a subset is a step over the whole folder applied to part of it:
#' the position a model has in the whole manifest is what is passed as
#' \code{index} and what numbers its seed, so a run taken model by model draws
#' what a run over all of them draws. That is what makes a long estimation
#' resumable, and the estimation steps for a folder pass \code{models} on.
#'
#' With more than one core, \code{f} is sent to the workers as it stands, and
#' the workers load the installed package. A function that calls something the
#' package does not export, or something only a newer version of it has, fails
#' there rather than in the session it was written in. Pass what such a function
#' needs as an argument instead.
#'
#' For the posterior itself, prefer \code{\link{bayests_files}} over a cluster:
#' the draws are the largest thing a model has, and reading them into R and
#' writing them out again costs more than the sampler. An idle R worker is not
#' cheap either -- on Windows a fresh one commits about 2 GB before it does
#' anything and about 4 GB once it has loaded this package.
#'
#' @return With \code{write = TRUE}, the object in \code{x} with its manifest
#' brought up to date, invisibly. With \code{write = FALSE}, a named list of
#' what \code{f} returned, one element per model.
#'
#' @examples
#'
#' data("e1")
#' train <- diff(log(e1)) * 100
#'
#' models <- create_bvarmodel(data = train, p = 1:2, deterministic = "const",
#'                            iterations = 20, burnin = 10)
#'
#' folder <- file.path(tempdir(), "bvartools-example-map")
#' unlink(folder, recursive = TRUE)
#' dir.create(folder, recursive = TRUE)
#' write_to_hdf5(models, folder = folder)
#'
#' stored <- open_models(folder)
#'
#' # What add_priors() for a folder does
#' stored <- map_models(stored, function(model) {
#'   add_priors(model, coef = list(v_i = 1), sigma = list(df = 3, scale = 1))
#' })
#'
#' # Reading only: how many observations each model was built on
#' map_models(stored, function(model) nrow(model$data$train$y), write = FALSE)
#'
#' @family model comparison
#' @export
map_models <- function(x, f, ..., models = NULL, write = TRUE, cores = 1,
                       export = NULL) {

  if (!inherits(x, "bvarfolder")) {
    stop("Argument 'x' must be of class 'bvarfolder'. Use open_models().")
  }
  f <- match.fun(f)

  manifest <- x[["manifest"]]
  rows <- seq_len(nrow(manifest))
  if (!is.null(models)) {
    unknown <- setdiff(models, manifest[["model"]])
    if (length(unknown) > 0) {
      stop("This folder does not hold the model(s) ",
           paste0(unknown, collapse = ", "), ".")
    }
    rows <- which(manifest[["model"]] %in% models)
  }
  if (length(rows) == 0) {
    return(if (write) invisible(x) else list())
  }

  shared <- manifest[["file"]][rows][duplicated(manifest[["file"]][rows]) |
                                       nzchar(manifest[["group"]][rows])]
  if (write && length(shared) > 0) {
    stop("The file(s) ", paste0(unique(shared), collapse = ", "),
         " hold more than one model. A step cannot write one of them back ",
         "without rewriting the others, so write such models to files of ",
         "their own first.")
  }

  # Everything a model needs travels in its own task -- the function and the
  # arguments included -- so that nothing is passed alongside the tasks and the
  # worker function takes one argument.
  arguments <- list(...)
  tasks <- lapply(rows, function(i) {
    list(file = file.path(x[["folder"]], manifest[["file"]][i]), index = i,
         group = manifest[["group"]][i], write = write, f = f,
         arguments = arguments)
  })

  cores <- max(1L, as.integer(cores))
  cores <- min(cores, length(tasks))

  specs <- if (cores < 2L) {
    lapply(tasks, .map_one_stored_model)
  } else {
    cluster <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, suppressMessages(library(bvartools)))
    if (!is.null(export)) {
      parallel::clusterExport(cluster, export, envir = parent.frame())
    }
    parallel::parLapplyLB(cluster, tasks, .map_one_stored_model)
  }

  if (!write) {
    names(specs) <- manifest[["model"]][rows]
    return(specs)
  }

  for (j in seq_along(specs)) {
    i <- rows[j]
    manifest[i, ] <- .model_row(specs[[j]], manifest[["model"]][i],
                                manifest[["file"]][i], manifest[["group"]][i])
  }
  x[["manifest"]] <- manifest

  invisible(x)
}

# One model: read, apply, write back, and return its specification for the
# manifest. The model itself never leaves the worker.
.map_one_stored_model <- function(task) {

  f <- task[["f"]]
  model <- bvartools::read_model_from_hdf5(filename = task[["file"]],
                                           group = task[["group"]])

  arguments <- c(list(model), task[["arguments"]])
  if ("index" %in% names(formals(f))) {
    arguments <- c(arguments, list(index = task[["index"]]))
  }
  model <- do.call(f, arguments)

  if (!isTRUE(task[["write"]])) {
    return(model)
  }

  if (!inherits(model, "bvarmodel") && !inherits(model, "bvecmodel")) {
    stop("The function returned an object of class ",
         paste(class(model), collapse = ", "), " for ", task[["file"]],
         " rather than a model.")
  }

  # Through a temporary file, so that an interrupted write cannot leave a model
  # that is neither the old one nor the new one. The writer refuses an existing
  # file, which is why the target is removed rather than written over.
  temporary <- paste0(task[["file"]], ".writing")
  if (file.exists(temporary)) {
    unlink(temporary)
  }
  bvartools::write_to_hdf5(model, filename = temporary)
  unlink(task[["file"]])
  if (!file.rename(temporary, task[["file"]])) {
    stop("The model of ", task[["file"]], " was written to ", temporary,
         " and could not be moved into place.")
  }

  model[["model"]]
}
