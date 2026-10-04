#' Write and Read the Tree of a BayesTS Model File
#'
#' Writes a nested list into an HDF5 file as a model's tree, or reads such a
#' tree back as a nested list. These are the steps beneath
#' \code{\link{write_to_hdf5}} and \code{\link{read_model_from_hdf5}}, exported
#' for packages that keep models of their own in BayesTS files -- the dynamic
#' factor models of \pkg{dfmtools}, for instance.
#'
#' @param tree a named list. An element that is itself a list is a group; any
#' other element is a dataset. See 'Details'.
#' @param filename path to the HDF5 file.
#' @param group the group the tree hangs under inside its file. Defaults to
#' \code{""}, the root of the file, which is where a file holding a single
#' model puts it. The spelling is that of the BayesTS command line's
#' \code{--group} flag: a leading slash and no trailing slash.
#' @param draws the draws to read from the blocks of \code{/posterior}, as
#' their positions in the chain. Defaults to \code{NULL}, every draw; see
#' \code{\link{read_model_from_hdf5}}.
#'
#' @details
#'
#' A model package translates its object into the layout BayesTS reads --
#' \code{/model}, \code{/data}, \code{/priors}, \code{/initial} and
#' \code{/posterior}, named and shaped as the sampler expects -- and hands the
#' result to \code{write_bayests_tree()}. What is left to this function is what
#' every model shares: refusing to overwrite a model, undoing a write that fails
#' half-way, and the attributes that let a value be read back as what it was.
#' \strong{It does not check the tree against any sampler}; a file that
#' BayesTS refuses is the translation's mistake, and \code{bayests check} is
#' the way to find it.
#'
#' The tree is mapped onto the file as follows:
#' \itemize{
#'   \item An element that is a list becomes a group of the same name, and
#'   its elements go below it. \code{NULL} elements are left out.
#'   \item The element \code{.attributes} of a list is not a group but the
#'   attributes of the group it is in: a named list of values, each written
#'   as an attribute. This is where \code{/model} keeps the specification --
#'   \code{algorithm}, \code{k}, \code{p} and the rest -- and where
#'   \code{rclass}, the class \code{\link{read_model_from_hdf5}} returns,
#'   belongs, and \code{rpackage}, the package that defines that class,
#'   whose namespace the reader loads for its methods.
#'   \item Any other element becomes a dataset. HDF5 stores dimensions in the
#'   reverse of R's order, so a matrix of \code{tt} rows and \code{k} columns
#'   is a \code{(k, tt)} dataset, which is what BayesTS calls \code{(tt, k)}
#'   "on paper". A vector is a one-dimensional dataset and is marked so that
#'   it is read back as a vector. A time series carries its column names and
#'   \code{tsp}, and draws of class \code{\link[coda]{mcmc}} their start, end
#'   and thinning interval, so that both are read back as what they were.
#' }
#'
#' Without a \code{group} the tree is the whole file, so an existing file is
#' refused; with one, only that group has to be free, which is how a second
#' model is added beside the first. A write that cannot be completed raises the
#' error and undoes what the call created.
#'
#' \code{read_bayests_tree()} is the inverse, for any model file, including one
#' written by the BayesTS command line: groups become lists, the attributes of
#' each group its element \code{.attributes}, and datasets values, restored
#' through the attributes above where the file carries them and as matrices
#' where it does not.
#'
#' A file whose \code{/model} carries \code{rclass} is read by
#' \code{\link{read_model_from_hdf5}} through \code{\link{from_bayests_tree}},
#' which is the method a model package provides to turn the tree back into its
#' object.
#'
#' @return \code{write_bayests_tree()} returns the path to the written file,
#' invisibly. \code{read_bayests_tree()} returns the tree as a nested list.
#'
#' @examples
#'
#' path <- tempfile(fileext = ".h5")
#' tree <- list(
#'   "model" = list(".attributes" = list("algorithm" = "Example", "k" = 2L)),
#'   "data" = list("train" = list("y" = matrix(rnorm(20), 10, 2)))
#' )
#' write_bayests_tree(tree, filename = path)
#' str(read_bayests_tree(path))
#'
#' @seealso \code{\link{from_bayests_tree}}, \code{\link{write_to_hdf5}},
#' \code{\link{read_model_from_hdf5}}.
#'
#' @export
write_bayests_tree <- function(tree, filename, group = "") {

  if (!is.list(tree) || (length(tree) > 0 && is.null(names(tree)))) {
    stop("Argument 'tree' must be a named list, one element per group or ",
         "dataset at the top of the model's tree.")
  }

  .hdf5_write_model(filename, group, function(handles, output) {
    .hdf5_write_tree(handles, output, tree, "")
  })
}

#' @rdname write_bayests_tree
#' @export
read_bayests_tree <- function(filename, group = "", draws = NULL) {

  group <- .normalize_hdf5_group(group)
  draws <- .check_read_draws(draws)

  h5_file <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (h5_file$is_valid) h5_file$close_all(), add = TRUE)

  .hdf5_read_tree(.hdf5_model_root(h5_file, group), draws, "")
}

#' Turn the Tree of a Model File into a Model
#'
#' The step of \code{\link{read_model_from_hdf5}} that turns what
#' \code{\link{read_bayests_tree}} reads into the object of a model package.
#'
#' @param tree the tree of a model file, as \code{\link{read_bayests_tree}}
#' returns it, with the class recorded in the file's \code{/model/rclass}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#'
#' \code{\link{read_model_from_hdf5}} reads the VAR and VEC models of this
#' package itself. For a file whose \code{rclass} names any other class it
#' reads the tree, gives it that class and calls this generic, so a package
#' that writes its models with \code{\link{write_bayests_tree}} reads them
#' back by registering a method for its class. The method undoes the
#' translation its writer made: names, orderings and shapes that BayesTS
#' wants and the package's object does not.
#'
#' The package named in the file's \code{/model/rpackage} is loaded first,
#' if it is installed, so that its methods are registered whether or not it
#' is attached. A file of the BayesTS command line carries no such attribute;
#' one of a factor model is read with \pkg{dfmtools} loaded.
#'
#' Without a method the default refuses, and names the class it found no
#' method for.
#'
#' \strong{A method keeps the attributes of \code{/model} it does not know
#' in element \code{model} of what it returns.} \code{\link{write_to_hdf5}}
#' records there where a model stood in the 'modellist' it was written from,
#' and \code{\link{read_models_from_folder}} rebuilds the list from it.
#'
#' @return An object of the class of \code{tree}, as its package defines it.
#'
#' @export
from_bayests_tree <- function(tree, ...) {
  UseMethod("from_bayests_tree")
}

#' @rdname from_bayests_tree
#' @export
from_bayests_tree.default <- function(tree, ...) {
  stop("There is no from_bayests_tree() method for class '", class(tree)[1],
       "', so the model file cannot be read as one. Load the package that ",
       "defines the class, or read the raw tree with read_bayests_tree().")
}


# Opens `filename` for the model under `group`, runs `write(handles, output)`
# against it, and closes it again -- or, if `write` fails, undoes what the call
# created.
#
# Every writer of a model shares this. What must not already be there is the
# model, not the file: without a group a model is the whole file, so an
# existing file is refused; with one, a file that already holds other models is
# exactly what is being added to, and only that group has to be free.
#
# The writers used to run inside a try() that discarded its result, so a full
# disk, an unwritable path or a malformed element left a half-written file that
# looked finished. The error reaches the caller now. A file this call made is
# removed whole; a group it added to a file that was already there is unlinked
# on its own, so the models beside it survive. Either way the handle is closed,
# or the file stays locked for the rest of the session. HDF5 does not reclaim
# the space of an unlinked group, but the name is free again, which is what a
# retry needs.
.hdf5_write_model <- function(filename, group, write) {

  group <- .normalize_hdf5_group(group)

  if (dir.exists(filename)) {
    stop("Argument 'filename' is not a path to a file.")
  }

  file_existed <- file.exists(filename)
  if (group == "" && file_existed) {
    stop(paste0("File ", filename, " already exists."))
  }

  h5_file <- hdf5r::h5file(filename, mode = "a")

  completed <- FALSE
  group_existed <- FALSE
  on.exit({
    if (!completed && group != "" && !group_existed && h5_file$is_valid) {
      try(h5_file$link_delete(group), silent = TRUE)
    }
    if (h5_file$is_valid) {
      h5_file$close_all()
    }
    if (!completed && !file_existed) {
      unlink(filename)
    }
  }, add = TRUE)

  if (group != "") {
    group_existed <- .hdf5_exists(h5_file, group)
    if (group_existed) {
      stop(paste0("Group ", group, " of file ", filename, " already exists."))
    }
  }

  # Every path below is named against this rather than against the file, which
  # is all a group amounts to on the way out. What the write opens below it is
  # collected as it goes, so that it can be closed again without asking the
  # file what is open.
  handles <- .hdf5_handles()
  output <- .create_hdf5_group(h5_file, group)
  if (inherits(output, "H5Group")) {
    .hdf5_keep(handles, output)
  }

  write(handles, output)

  .hdf5_close(handles, h5_file)
  completed <- TRUE

  invisible(filename)
}

# Writes the elements of `node` below `parent`. `path` is where `parent` is in
# the tree, for the messages.
.hdf5_write_tree <- function(handles, parent, node, path) {

  # By position, since names() of a list without any is NULL and looping over
  # it would skip an unnamed element rather than refuse it.
  labels <- names(node)
  if (is.null(labels)) {
    labels <- rep("", length(node))
  }
  for (position in seq_along(node)) {
    i <- labels[position]
    value <- node[[position]]
    if (is.null(value) || identical(i, ".attributes")) {
      next
    }
    if (is.na(i) || i == "") {
      stop("An element of the tree below '", if (path == "") "/" else path,
           "' has no name, and every group and dataset needs one.")
    }
    if (is.list(value)) {
      group <- .hdf5_group(handles, parent, i)
      .hdf5_write_tree(handles, group, value, paste0(path, "/", i))
    } else {
      .hdf5_write(parent, i, value, .hdf5_tree_attrs(value))
    }
  }

  attrs <- node[[".attributes"]]
  if (!is.null(attrs)) {
    if (!is.list(attrs) || (length(attrs) > 0 && is.null(names(attrs)))) {
      stop("Element '.attributes' of '", if (path == "") "/" else path,
           "' must be a named list of values.")
    }
    for (i in names(attrs)) {
      if (!is.null(attrs[[i]])) {
        .hdf5_write_attr(parent, i, attrs[[i]])
      }
    }
  }

  invisible(NULL)
}

# The attributes a dataset of a tree is written with, decided by what it is: a
# chain carries its mcpar, a time series its names and tsp. Anything else has
# none beyond the mark .hdf5_write() puts on a vector.
.hdf5_tree_attrs <- function(value) {

  if (coda::is.mcmc(value)) {
    return(.hdf5_draws_attrs(value))
  }
  if (stats::is.ts(value)) {
    return(.hdf5_series_attrs(value))
  }

  NULL
}

# Reads the group `node` into a list. `draws` applies to the blocks below
# /posterior and nowhere else, since only they are chains.
.hdf5_read_tree <- function(node, draws, path) {

  result <- list()

  attrs <- hdf5r::h5attributes(node)
  if (length(attrs) > 0) {
    result[[".attributes"]] <- attrs
  }

  in_posterior <- path == "/posterior" || startsWith(path, "/posterior/")
  for (i in names(node)) {
    element <- node[[i]]
    child <- paste0(path, "/", i)
    if (inherits(element, "H5Group")) {
      result[[i]] <- .hdf5_read_tree(element, draws, child)
    } else if (in_posterior) {
      result[[i]] <- .read_posterior_block(element, draws)
    } else {
      result[[i]] <- .hdf5_read_tree_value(element)
    }
  }

  result
}

# A dataset outside the posterior, as what it was written as: a time series
# where it carries the attributes of one, a vector where it is marked as one,
# and a matrix otherwise.
.hdf5_read_tree_value <- function(dataset) {

  attrs <- hdf5r::h5attr_names(dataset)
  if (!all(c("tsp", "variables") %in% attrs)) {
    return(.hdf5_read_value(dataset))
  }

  series <- stats::ts(.hdf5_read_matrix(dataset))
  variables <- hdf5r::h5attr(dataset, "variables")
  dimnames(series) <- list(NULL, variables)
  stats::tsp(series) <- hdf5r::h5attr(dataset, "tsp")
  series <- .hdf5_restore_class(series, dataset)
  for (name in c("scale", "centre")) {
    if (name %in% attrs) {
      factors <- hdf5r::h5attr(dataset, name)
      names(factors) <- variables
      attr(series, name) <- factors
    }
  }

  series
}
