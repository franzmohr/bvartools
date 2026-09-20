#' Open a Folder of Stored Models
#'
#' Opens a folder written with \code{\link{write_to_hdf5}} and returns a handle
#' to the models in it, which the estimation steps work on one at a time rather
#' than reading all of them into the session.
#'
#' @param folder path to a folder holding models in HDF5 files.
#'
#' @details
#' \code{\link{read_models_from_folder}} reads a folder into a list, draws and
#' all, which is what to do with models that fit. Many do not. A lag and rank
#' grid is dozens of models, an expanding window is one per quarter, and a model
#' with time varying coefficients keeps a coefficient path per draw, so a folder
#' can be tens of gigabytes while each model in it is a few hundred megabytes.
#'
#' A handle carries what the folder holds -- a row per model, saying where it is
#' and what it is -- and none of the draws. \code{\link{map_models}} and the
#' estimation steps for a handle then read one model, change it, write it back
#' and drop it, so a step costs one model per worker rather than the whole
#' folder. \code{\link{open_model}} is the same idea for one model whose chain is
#' too long to hold: a handle to a file, read a piece at a time.
#'
#' The manifest is built by reading the specification of every model, which is a
#' set of attributes rather than any of its data, so opening a folder is cheap.
#' A model is named by where its file sits below \code{folder}, without the
#' extension, and by its group where a file holds more than one -- the names
#' \code{\link{read_models_from_folder}} gives.
#'
#' @return A list of class 'bvarfolder' with the elements \code{folder} and
#' \code{manifest}, a data frame with the columns \code{model}, \code{file},
#' \code{group} and the specification of each model.
#'
#' @examples
#'
#' data("e1")
#' train <- diff(log(e1)) * 100
#'
#' models <- create_bvarmodel(data = train, p = 1:3, deterministic = "const",
#'                            iterations = 20, burnin = 10)
#'
#' folder <- file.path(tempdir(), "bvartools-example-folder")
#' unlink(folder, recursive = TRUE)
#' dir.create(folder, recursive = TRUE)
#' write_to_hdf5(models, folder = folder)
#'
#' stored <- open_models(folder)
#' stored
#' stored[["manifest"]][, c("model", "p", "iterations")]
#'
#' @family model comparison
#' @export
open_models <- function(folder) {

  if (!dir.exists(folder)) {
    stop("Folder ", folder, " does not exist.")
  }

  files <- sort(list.files(folder, recursive = TRUE, full.names = FALSE))
  files <- files[grepl("[.]h5$", files, ignore.case = TRUE)]
  if (length(files) == 0) {
    stop("Folder ", folder, " holds no .h5 files.")
  }

  rows <- list()
  for (file in files) {
    filename <- file.path(folder, file)
    stem <- tools::file_path_sans_ext(file)
    for (group in list_models_in_hdf5(filename)) {
      name <- if (group == "") stem else paste0(stem, ":", group)
      rows[[length(rows) + 1]] <-
        .model_row(.model_spec_in_hdf5(filename, group), name, file, group)
    }
  }

  if (length(rows) == 0) {
    stop("Folder ", folder, " holds no models.")
  }

  manifest <- do.call("rbind", rows)
  rownames(manifest) <- NULL

  structure(list(folder = folder, manifest = manifest),
            class = c("bvarfolder", "list"))
}

#' @rdname open_models
#' @param x an object of class 'bvarfolder'.
#' @param ... further arguments passed to or from other methods.
#' @export
print.bvarfolder <- function(x, ...) {

  manifest <- x[["manifest"]]
  cat(nrow(manifest), "models in", x[["folder"]], "\n")

  algorithms <- table(manifest[["algorithm"]])
  for (algorithm in names(algorithms)) {
    cat(" ", algorithms[[algorithm]], algorithm, "\n")
  }

  invisible(x)
}

#' @rdname open_models
#' @param models character vector of the models to name, as the manifest names
#' them, or \code{NULL} for all of them.
#' @details
#' \code{model_files} names the files the models lie in, which is what a program
#' that works on the files is pointed at. With
#' \code{\link{bayests_files}} a folder is estimated as
#' \preformatted{run <- bayests_files(executable = "/opt/bayests/bin/bayests")
#' run(model_files(stored), jobs = 6)}
#' and the draws never pass through R. Files are named once however many models
#' they hold, since BayesTS works through every model in a file it is given.
#'
#' @export
model_files <- function(x, models = NULL) {

  if (!inherits(x, "bvarfolder")) {
    stop("Argument 'x' must be of class 'bvarfolder'. Use open_models().")
  }

  manifest <- .select_models(x[["manifest"]], models)
  files <- unique(manifest[["file"]])
  stats::setNames(file.path(x[["folder"]], files), files)
}

# The rows of the manifest a step is applied to, in the order of the manifest.
.select_models <- function(manifest, models) {

  if (is.null(models)) {
    return(manifest)
  }
  unknown <- setdiff(models, manifest[["model"]])
  if (length(unknown) > 0) {
    stop("This folder does not hold the model(s) ",
         paste0(unknown, collapse = ", "), ".")
  }
  manifest[manifest[["model"]] %in% models, , drop = FALSE]
}

# The specification of a stored model, read from the attributes of its /model
# group. Neither its data nor its priors are touched, so this costs one open of
# the file however large the model is.
.model_spec_in_hdf5 <- function(filename, group = "") {

  group <- .normalize_hdf5_group(group)
  file <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (file$is_valid) file$close_all(), add = TRUE)
  root <- .hdf5_model_root(file, group)

  if (!"model" %in% names(root)) {
    stop("The model in ", filename, " has no specification.")
  }
  hdf5r::h5attributes(root[["model"]])
}

# One row of the manifest. The columns are the ones every model has, so that a
# folder of mixed models still gives a table; what a particular kind of model
# carries beyond them is in the model itself.
.model_row <- function(specs, name, file, group) {

  field <- function(name, mode) {
    value <- specs[[name]]
    if (is.null(value) || length(value) != 1) {
      return(as.vector(NA, mode = mode))
    }
    as.vector(value, mode = mode)
  }

  data.frame(
    model = name,
    file = file,
    group = group,
    algorithm = field("algorithm", "character"),
    type = field("type", "character"),
    k = field("k", "integer"),
    p = field("p", "integer"),
    m = field("m", "integer"),
    s = field("s", "integer"),
    n = field("n", "integer"),
    # Error correction models only. NA for a VAR, which has no rank.
    rank = field("rank", "integer"),
    varsel = field("varsel", "character"),
    structural = field("structural", "logical"),
    tvp = field("tvp", "logical"),
    error = field("error", "character"),
    iterations = field("iterations", "integer"),
    burnin = field("burnin", "integer"),
    # A model that keeps every draw carries no 'thin', so that is one rather
    # than missing.
    thin = if (is.null(specs[["thin"]])) 1L else field("thin", "integer"),
    # The discounted models only. One is the value at which the quantity a
    # discount governs does not move, so a model that carries neither is
    # reported as one rather than as missing -- and every other algorithm is a
    # model whose coefficients and error covariance are what they are, which is
    # what one says.
    delta_beta = if (is.null(specs[["delta_beta"]])) 1 else field("delta_beta", "numeric"),
    delta_sigma = if (is.null(specs[["delta_sigma"]])) 1 else field("delta_sigma", "numeric"),
    stringsAsFactors = FALSE
  )
}
