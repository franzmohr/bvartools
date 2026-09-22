#' Import Models from a Folder of HDF5 Files
#'
#' Imports every model stored below a folder.
#'
#' @param folder Path to a folder containing model data.
#' @param draws the draws to read of every model, as
#' \code{\link{read_model_from_hdf5}} takes them: \code{NULL}, the default, for
#' the whole chain, a vector of positions for those draws, or
#' \code{integer(0)} for none of them.
#'
#' @details
#'
#' The folder is walked recursively, and each HDF5 file in it is asked which of
#' its groups hold a model -- so a file holding several models contributes all
#' of them. See \code{\link{list_models_in_hdf5}} for what counts as one.
#'
#' The result is a flat named list rather than a nesting that mirrors the
#' directory tree. Each name is the path of the file relative to \code{folder}
#' without its extension, followed by the group where a model does not sit at
#' the root of its file. Names are what a caller needs to tell one model from
#' another -- which sub-model of a global model it is, say -- and a nesting
#' whose depth depended on where the caller pointed could not provide them.
#'
#' The one exception is an expanding window: the windows that
#' \code{\link{write_to_hdf5}} wrote to one directory come back as one element
#' of class 'expandingwindow', named after that directory, so that the windows of
#' different specifications are not pooled. A folder that holds a single
#' expanding window is returned as that expanding window. A model list written
#' by \code{\link{write_to_hdf5}} comes back in the order it was written in,
#' which the files record; other models are in the order of their names.
#'
#' The class of each model comes from the \code{rclass} attribute the writer
#' records, not from its file name.
#'
#' \code{draws} is what makes a folder larger than the session readable: with
#' \code{integer(0)} the models come back with their specification, their data
#' and their priors and no draws, and with a vector of positions the chain is
#' read a piece at a time. For working on the models rather than reading them,
#' \code{\link{open_models}} hands them over one at a time instead.
#'
#' @return A named list of class 'modellist', or of class 'expandingwindow' if
#' the models say they belong to one.
#'
#' @examples
#'
#' # Load data
#' data("e1")
#' train <- diff(log(e1)) * 100
#'
#' # Create and store models
#' folder <- file.path(tempdir(), "models")
#' dir.create(folder, showWarnings = FALSE)
#' for (p in 1:2) {
#'   model <- create_bvarmodel(data = train, p = p, deterministic = "const",
#'                             iterations = 10, burnin = 10)
#'   write_to_hdf5(model, filename = file.path(folder, paste0("model-", p, ".h5")))
#' }
#'
#' models <- read_models_from_folder(folder)
#' names(models)
#'
#' @export
read_models_from_folder <- function(folder, draws = NULL) {

  if (!dir.exists(folder)) {
    stop("Specified folder does not exist.")
  }

  # Relative and absolute forms of the same list, so that a model can be named
  # by where it sits below 'folder' while being read from its full path.
  relative_paths <- sort(list.files(folder, recursive = TRUE, full.names = FALSE))
  relative_paths <- relative_paths[grepl("[.]h5$", relative_paths, ignore.case = TRUE)]

  if (length(relative_paths) == 0) {
    stop("Specified folder does not contain any .h5 files.")
  }

  result <- list()
  result_names <- character()
  result_files <- character()

  for (i in relative_paths) {

    filename <- file.path(folder, i)

    # A name that does not depend on how deep the caller pointed: the file's
    # place below 'folder', and the group where the file holds more than one
    # model.
    stem <- tools::file_path_sans_ext(i)

    for (group in list_models_in_hdf5(filename)) {
      result[[length(result) + 1]] <- read_model_from_hdf5(filename = filename,
                                                           group = group,
                                                           draws = draws)
      result_names <- c(result_names,
                        if (group == "") stem else paste0(stem, ":", group))
      result_files <- c(result_files, i)
    }
  }

  if (length(result) == 0) {
    stop("Specified folder does not contain any models.")
  }

  names(result) <- result_names

  # Which models belong together is read off the models themselves: the
  # windows of an expanding window carry the class of their collection, which
  # write_to_hdf5() stores. It used to be guessed from "ExpWind" in the first
  # file's path, which made the class depend on the name of a directory; files
  # written before the attribute existed still carry the name, so it is kept as
  # a fallback. Where a model list was written, every model also carries the
  # position of its element in that list.
  directories <- dirname(result_files)
  in_window <- vapply(seq_along(result), function(i) {
    collection <- result[[i]][["model"]][["rclass_collection"]]
    if (is.null(collection)) {
      grepl("ExpWind", result_files[i], fixed = TRUE)
    } else {
      "expandingwindow" %in% collection
    }
  }, logical(1))
  index <- vapply(result, function(x) {
    i <- x[["model"]][["rindex_collection"]]
    if (is.null(i)) NA_real_ else as.numeric(i)[1]
  }, numeric(1))

  for (i in seq_along(result)) {
    result[[i]][["model"]][["rclass_collection"]] <- NULL
    result[[i]][["model"]][["rindex_collection"]] <- NULL
  }

  # Every directory of windows becomes one expanding window, in the place of
  # its first window; every other model stays an element of its own. The
  # windows of different specifications were once returned as one long
  # expanding window, which pooled them.
  key <- ifelse(in_window, paste0("window:", directories), paste0("model:", seq_along(result)))
  groups <- unique(key)
  elements <- lapply(groups, function(g) {
    members <- which(key == g)
    if (startsWith(g, "window:")) {
      windows <- result[members]
      class(windows) <- c("expandingwindow", "list")
      windows
    } else {
      result[[members]]
    }
  })
  names(elements) <- vapply(groups, function(g) {
    members <- which(key == g)
    if (startsWith(g, "window:")) directories[members[1]] else result_names[members]
  }, character(1))
  element_index <- vapply(groups, function(g) index[which(key == g)[1]], numeric(1))

  # The order of the list that was written, where the files say what it was;
  # the order of the file names otherwise. order() keeps ties, and the
  # elements without an index, in the order they were read.
  elements <- elements[order(element_index, na.last = TRUE)]

  # A folder that holds one expanding window is that expanding window, however
  # deep it sits below 'folder'.
  if (length(elements) == 1 && inherits(elements[[1]], "expandingwindow")) {
    return(elements[[1]])
  }

  class(elements) <- c("modellist", "list")
  return(elements)
}
