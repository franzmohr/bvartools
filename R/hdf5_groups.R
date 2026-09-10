# Addressing a model inside an HDF5 file.
#
# A model's tree -- /model, /data, /priors, /initial, /posterior -- can sit at
# the root of its file, which is where a file holding a single model puts it, or
# under a group, which is what lets one file hold several side by side.
#
# The spelling of a group name is the same one BayesTS uses, so that a group
# named here and a --group passed there mean the same thing: a leading slash, no
# trailing slash, and "" for the root of the file.

# The normalized form of a group name.
#
# Idempotent, so a name that has already been through it can be passed again.
# Throws on a name that cannot be a group -- an empty component, or a "." or
# ".." that HDF5 would resolve relative to somewhere the caller did not name.
# Rejected before a file is opened, because that is a mistake in the call rather
# than a problem with the file.
.normalize_hdf5_group <- function(group) {

  if (is.null(group)) {
    return("")
  }

  if (!is.character(group) || length(group) != 1 || is.na(group)) {
    stop("Argument 'group' must be a single character string.")
  }

  if (group == "" || group == "/") {
    return("")
  }

  if (substring(group, 1, 1) != "/") {
    group <- paste0("/", group)
  }

  while (nchar(group) > 1 && substring(group, nchar(group)) == "/") {
    group <- substring(group, 1, nchar(group) - 1)
  }

  components <- strsplit(substring(group, 2), "/", fixed = TRUE)[[1]]

  if (length(components) == 0 || any(components == "")) {
    stop("Group '", group, "' has an empty component.")
  }

  if (any(components %in% c(".", ".."))) {
    stop("Group '", group, "' contains '.' or '..', which would resolve ",
         "relative to somewhere else in the file.")
  }

  return(group)
}

# Whether `path` names something in the file.
#
# hdf5r's own exists() raises rather than returning FALSE when a level above
# `path` is missing, which is exactly the case a caller asking "is this group
# there?" has in mind. Each level is therefore checked in turn, and the first
# one that is absent answers the question.
.hdf5_exists <- function(h5_file, path) {

  if (path == "" || path == "/") {
    return(TRUE)
  }

  components <- strsplit(substring(path, 2), "/", fixed = TRUE)[[1]]

  current <- ""
  for (i in components) {
    current <- paste0(current, "/", i)
    present <- tryCatch(h5_file$exists(current), error = function(e) FALSE)
    if (!isTRUE(present)) {
      return(FALSE)
    }
  }

  return(TRUE)
}

# The object a model's paths hang off: the file itself for the root, otherwise
# the group. Everything above this layer then names "model", "data" and the rest
# against it without caring which of the two it got.
.hdf5_model_root <- function(h5_file, group) {

  if (group == "") {
    return(h5_file)
  }

  if (!.hdf5_exists(h5_file, group)) {
    stop("File ", h5_file$get_filename(), " does not contain group ", group, ".")
  }

  root <- h5_file[[group]]

  if (!inherits(root, "H5Group")) {
    stop("'", group, "' is not a group in ", h5_file$get_filename(), ".")
  }

  return(root)
}

# Creates `group` and every level above it that is not there yet, and returns it.
#
# hdf5r does not create intermediate groups on the way to a nested one, so each
# level is made in turn.
.create_hdf5_group <- function(h5_file, group) {

  if (group == "") {
    return(h5_file)
  }

  components <- strsplit(substring(group, 2), "/", fixed = TRUE)[[1]]

  path <- ""
  for (i in components) {
    path <- paste0(path, "/", i)
    if (!.hdf5_exists(h5_file, path)) {
      # Closed again straight away. The handle of a group made on the way to a
      # deeper one is of no use to anybody, and left open it is one more object
      # the file carries until it is closed.
      created <- h5_file$create_group(path)
      created$close()
    }
  }

  return(h5_file[[group]])
}

# Whether `group` holds a model: a "model" subgroup carrying an "algorithm"
# attribute.
#
# The same rule BayesTS applies, so both sides agree on what a model is. Reading
# the attribute rather than only looking for the group is what makes this answer
# the question a caller actually has -- can this be read as a model -- instead of
# the weaker "something here is called model".
.is_hdf5_model_group <- function(h5_file, group) {

  path <- paste0(group, "/model")

  if (!.hdf5_exists(h5_file, path)) {
    return(FALSE)
  }

  object <- h5_file[[path]]

  if (!inherits(object, "H5Group")) {
    return(FALSE)
  }

  return("algorithm" %in% hdf5r::h5attr_names(object))
}

# Recursive worker behind list_models_in_hdf5.
.collect_hdf5_model_groups <- function(h5_file, group) {

  if (.is_hdf5_model_group(h5_file, group)) {
    # Stop here. A model's own data, priors and posterior are its subtree, not a
    # place further models could be, and not descending into them is what keeps
    # this from walking every dataset in a file that holds a hundred.
    return(group)
  }

  handle <- if (group == "") h5_file else h5_file[[group]]

  found <- character()
  for (name in names(handle)) {
    child <- paste0(group, "/", name)
    if (inherits(h5_file[[child]], "H5Group")) {
      found <- c(found, .collect_hdf5_model_groups(h5_file, child))
    }
  }

  return(found)
}

#' Models in an HDF5 File
#'
#' Finds the groups of an HDF5 file that hold a model.
#'
#' @param filename path to an HDF5 file.
#' @param group the group to search under. Defaults to \code{""}, the root of
#' the file, which searches all of it.
#'
#' @details
#'
#' A model is a group with a \code{model} subgroup carrying an
#' \code{algorithm} attribute. The search stops at a model rather than
#' descending into its \code{data}, \code{priors} and \code{posterior}, which
#' are its own subtree and not a place further models could be.
#'
#' This is the same rule the BayesTS command line applies for its
#' \code{--all-groups} flag, so both agree on which groups of a file are models.
#'
#' @return A character vector of group names, sorted, empty if the file holds no
#' model. A model at the root of its file is reported as \code{""}, which is the
#' value \code{\link{read_model_from_hdf5}} takes for it.
#'
#' @examples
#'
#' # Load data
#' data("e1")
#' train <- diff(log(e1)) * 100
#'
#' # Create and store two models in one file
#' path_to_model <- tempfile(fileext = ".h5")
#' for (p in 1:2) {
#'   model <- create_bvarmodel(data = train, p = p, deterministic = "const",
#'                             iterations = 10, burnin = 10)
#'   write_to_hdf5(model, filename = path_to_model,
#'                 group = paste0("/models/", p))
#' }
#'
#' list_models_in_hdf5(path_to_model)
#'
#' @export
list_models_in_hdf5 <- function(filename, group = "") {

  group <- .normalize_hdf5_group(group)

  h5_file <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (h5_file$is_valid) h5_file$close_all(), add = TRUE)

  # A group that is not there is a mistake worth reporting, and is not the same
  # thing as a group that is there and holds no model.
  .hdf5_model_root(h5_file, group)

  found <- .collect_hdf5_model_groups(h5_file, group)

  # names() returns HDF5's order, which is not the caller's. Sorting makes the
  # result the same on every run and every platform.
  return(sort(found))
}
