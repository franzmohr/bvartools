#' Import Models from a Folder of HDF5 Files
#'
#' Imports every model stored below a folder.
#'
#' @param folder Path to a folder containing model data.
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
#' The class of each model comes from the \code{rclass} attribute the writer
#' records, not from its file name.
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
read_models_from_folder <- function(folder) {

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

  for (i in relative_paths) {

    filename <- file.path(folder, i)

    # A name that does not depend on how deep the caller pointed: the file's
    # place below 'folder', and the group where the file holds more than one
    # model.
    stem <- tools::file_path_sans_ext(i)

    for (group in list_models_in_hdf5(filename)) {
      result[[length(result) + 1]] <- read_model_from_hdf5(filename = filename,
                                                           group = group)
      result_names <- c(result_names,
                        if (group == "") stem else paste0(stem, ":", group))
    }
  }

  if (length(result) == 0) {
    stop("Specified folder does not contain any models.")
  }

  names(result) <- result_names

  # The kind of collection these models form is read off the models themselves.
  # It used to be guessed by looking for "ExpWind" in the first file's path,
  # which made the class of the result depend on the name of a directory -- and
  # on the names of every directory above it, since the match was against the
  # full path. Files written before the attribute existed still carry the name,
  # so that is kept as a fallback.
  collection <- unique(unlist(lapply(result, function(x) x[["model"]][["rclass_collection"]])))
  if (is.null(collection) && any(grepl("ExpWind", relative_paths, fixed = TRUE))) {
    collection <- c("expandingwindow", "list")
  }

  if (is.null(collection)) {
    class(result) <- c("modellist", "list")
  } else {
    class(result) <- collection
  }

  return(result)
}
