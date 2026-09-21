# Writing a model's tree into an HDF5 file.
#
# The tree -- /model, /data, /priors, /initial, /posterior -- is the same one
# hdf5_groups.R addresses on the way in. What is here is only about the cost of
# putting it there: a model is a few hundred kilobytes of numbers, and writing
# one used to spend more time looking things up in the file than writing them.

# The handles one write opens, so that closing them does not mean asking the
# file what it has open.
#
# H5File$close_all() answers that question by enumerating every open object of
# every kind, and the enumeration alone costs more per model than writing the
# numbers in it. What a write opened is something the write already knows.
.hdf5_handles <- function() {

  handles <- new.env(parent = emptyenv())
  handles[["open"]] <- list()

  return(handles)
}

# Records `handle` as open, and returns it.
.hdf5_keep <- function(handles, handle) {

  handles[["open"]] <- c(handles[["open"]], handle)

  return(handle)
}

# Closes what a write opened, innermost first, and then the file itself.
#
# The file is closed last and only if it is still open, because `output` is the
# file rather than a group when a model is written at the root.
.hdf5_close <- function(handles, h5_file) {

  for (handle in rev(handles[["open"]])) {
    if (handle$is_valid) {
      handle$close()
    }
  }

  if (h5_file$is_valid) {
    h5_file$close()
  }

  return(invisible(NULL))
}

# The group `name` below `parent`, created if it is not there yet.
#
# exists() answers with a single lookup, where the names() this replaces listed
# every link of the group to answer the same question.
.hdf5_group <- function(handles, parent, name) {

  if (parent$exists(name)) {
    return(.hdf5_keep(handles, parent[[name]]))
  }

  return(.hdf5_keep(handles, parent$create_group(name)))
}

# Attaches `value` to `object` as the attribute `name`, and closes it again.
#
# hdf5r's `h5attr<-` creates the attribute and leaves its handle open for the
# garbage collector. HDF5 keeps a file open for as long as anything in it is,
# closing the file itself included, so every attribute written that way kept
# the file open after the write had returned -- and on Windows locked against
# every other process, BayesTS among them -- until a gc() happened to come
# round. An attribute that is already there is replaced, as it was before.
.hdf5_write_attr <- function(object, name, value) {

  if (object$attr_exists(name)) {
    object$attr_delete(name)
  }

  attribute <- object$create_attr(name, robj = value,
                                  dtype = .hdf5_dtype(value),
                                  space = .hdf5_attr_space(value))
  attribute$close()

  return(invisible(NULL))
}

# The HDF5 type and dataspace of a value, built once rather than per write.
#
# hdf5r guesses both whenever it is not told them, and the guess is most of what
# a write costs: creating the attributes of one model -- some forty of them,
# the specification and the labels of every series -- took 84 per cent of
# write_to_hdf5(), almost all of it in guess_dtype(), the type factory it calls,
# and building a fresh dataspace. For a grid of thousands of models that was
# the whole of the run.
#
# What is cached is exactly what hdf5r would have guessed, so no file changes:
# guess_dtype() with the string length create_attr() and create_dataset() pass,
# Inf, gives a variable length C string for a character vector, H5T_LOGICAL with
# NA for a logical one, and the native int and double for the rest; and
# guess_space() for an attribute is a simple dataspace of the value's
# dimensions with maximum dimensions equal to them -- simple even for a single
# value, never scalar, which is what BayesTS has always been given. Anything
# else -- a factor, a list, a 64-bit integer, a complex number -- gets NULL, and
# hdf5r guesses as it always did.
#
# HDF5 identifiers belong to the process that made them, so the cache is
# emptied when the process id changes: a forked worker builds its own rather
# than use its parent's.
.hdf5_type_cache <- new.env(parent = emptyenv())

.hdf5_cache <- function() {
  if (!identical(.hdf5_type_cache[["pid"]], Sys.getpid())) {
    rm(list = ls(.hdf5_type_cache, all.names = TRUE), envir = .hdf5_type_cache)
    .hdf5_type_cache[["pid"]] <- Sys.getpid()
    .hdf5_type_cache[["types"]] <- new.env(parent = emptyenv())
    .hdf5_type_cache[["spaces"]] <- new.env(parent = emptyenv())
  }
  .hdf5_type_cache
}

.hdf5_dtype <- function(value) {

  if (is.factor(value) || is.list(value) || inherits(value, "integer64")) {
    return(NULL)
  }
  kind <- if (is.character(value)) {
    "character"
  } else if (is.logical(value)) {
    "logical"
  } else if (is.integer(value)) {
    "integer"
  } else if (is.double(value)) {
    "double"
  } else {
    return(NULL)
  }

  types <- .hdf5_cache()[["types"]]
  if (is.null(types[[kind]])) {
    types[[kind]] <- switch(kind,
                            character = hdf5r::H5T_STRING$new(type = "c", size = Inf),
                            logical = hdf5r::H5T_LOGICAL$new(include_NA = TRUE),
                            integer = hdf5r::h5types$H5T_NATIVE_INT,
                            double = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  }
  types[[kind]]
}

.hdf5_attr_space <- function(value) {

  # Only for the values .hdf5_dtype() has a type for: for anything else hdf5r
  # guesses the type, and the space has to be guessed with it.
  if (is.null(.hdf5_dtype(value))) {
    return(NULL)
  }

  dims <- if (is.null(dim(value))) length(value) else dim(value)
  key <- paste(dims, collapse = "x")
  spaces <- .hdf5_cache()[["spaces"]]
  if (is.null(spaces[[key]])) {
    spaces[[key]] <- hdf5r::H5S$new(type = "simple", dims = dims, maxdims = dims)
  }
  spaces[[key]]
}

# Writes the selection scheme of the covariance block where BayesTS reads it.
#
# The scheme of the coefficients is the attribute 'varsel' of /model. That of
# the covariance block of a time varying model is the attribute 'varsel' of
# /model/priors/psi, and the samplers look for it nowhere else. It used to be
# written only as the dataset /priors/psi/varsel, beside the inclusion priors,
# which the samplers never read: they took the scheme to be "none", ignored the
# inclusion priors and the starting values of the indicators, and ran without
# the selection over the covariances the model asked for. The dataset is still
# written with the rest of the prior, because that is what the reader rebuilds
# 'priors' from.
.hdf5_write_psi_varsel <- function(handles, group_model, prior_psi) {

  varsel <- prior_psi[["varsel"]]
  if (is.null(varsel)) {
    return(invisible(NULL))
  }

  group_priors <- .hdf5_group(handles, group_model, "priors")
  group_psi <- .hdf5_group(handles, group_priors, "psi")
  .hdf5_write_attr(group_psi, "varsel", varsel)

  return(invisible(NULL))
}

# Writes `value` into `group` under `name`, with `attrs` attached to it.
#
# The dataset is created once and its attributes are written on the handle that
# creation hands back. Assigning with `[[<-` and then naming the dataset again
# for each attribute opened the same dataset three times over, and an open costs
# more than the numbers written through it. Anything already under `name` is
# left as it is, which is what the writers did before.
#
# The handle is closed here rather than left to the end: a dataset is written
# once and never looked at again in the same call, and every one left open is
# one more the file has to carry.
.hdf5_write <- function(group, name, value, attrs = NULL) {

  # A model that has not been through add_priors() yet carries no value for the
  # hyperparameters its error specification calls for, and hdf5r cannot infer a
  # dataset type from NULL. Nothing is written, which is what the reader expects:
  # it rebuilds each group from the names that are actually in the file, so an
  # element that was never written simply comes back absent.
  if (is.null(value)) {
    return(invisible(NULL))
  }

  if (group$exists(name)) {
    return(invisible(NULL))
  }

  # The type is given rather than guessed, for the reason .hdf5_dtype() gives.
  # The dataspace is still hdf5r's, because a dataset is chunked and its space
  # has unlimited maximum dimensions, which .hdf5_attr_space() does not build.
  dataset <- group$create_dataset(name, value, dtype = .hdf5_dtype(value))

  # HDF5 hands every dataset back as an array, so a scalar hyperparameter, the
  # type of a prior or a vector of shapes came back as a matrix, and a round
  # trip did not return what was written. A value that had no dimensions is
  # therefore marked as such, and .hdf5_read_value() drops the dimensions again.
  # A file written elsewhere carries no mark and is read as it always was.
  if (is.null(dim(value))) {
    .hdf5_write_attr(dataset, "rshape", "vector")
  }

  for (i in names(attrs)) {
    .hdf5_write_attr(dataset, i, attrs[[i]])
  }

  dataset$close()

  return(invisible(NULL))
}

# The attributes that make a stored matrix a time series again on the way back.
#
# The error correction term of a model that has been through
# scale_error_correction() carries the factors it was divided by, and they are
# the only way back to the scale of the data, as are the means a centred term
# lost. Without them an exported model that was scaled or centred could never be
# interpreted again, so they travel with it. The names are not stored: they are
# the variable names, which are written anyway. The class is, because a series
# built by this package and one supplied by the user differ in whether it lists
# "array", and the reader cannot tell them apart.
.hdf5_series_attrs <- function(x) {

  result <- list("variables" = dimnames(x)[[2]], "tsp" = stats::tsp(x),
                 "rclass" = class(x))

  for (name in c("scale", "centre")) {
    if (!is.null(attr(x, name))) {
      result[[name]] <- unname(attr(x, name))
    }
  }

  return(result)
}

# The attributes that make stored draws a chain again on the way back.
.hdf5_draws_attrs <- function(x) {
  mcpar <- coda::mcpar(x)
  list("start" = mcpar[1], "end" = mcpar[2], "thin" = mcpar[3])
}

# A dataset of the priors or the starting values as R held it. A value that the
# writer marked as having had no dimensions comes back as a vector; anything
# else, and anything from a file written elsewhere, as a matrix.
.hdf5_read_value <- function(dataset) {

  value <- hdf5r::readDataSet(dataset)

  if ("rshape" %in% hdf5r::h5attr_names(dataset) &&
      identical(hdf5r::h5attr(dataset, "rshape"), "vector")) {
    return(as.vector(value))
  }

  return(as.matrix(value))
}

# The class a series had when it was written, where the file says so. Otherwise
# the series keeps the class the reader gave it.
.hdf5_restore_class <- function(series, dataset) {

  if ("rclass" %in% hdf5r::h5attr_names(dataset)) {
    class(series) <- hdf5r::h5attr(dataset, "rclass")
  }

  return(series)
}

# Name of a starting value inside /initial of a model file.
#
# R calls the starting error precision of every gamma model u_omega_inv. The
# model file follows BayesTS, whose constant coefficient gamma samplers
# (VarNormalGamma, VecNormalGamma) read /initial/u_sigma_inv and only the time
# varying ones /initial/u_omega_inv. The Rcpp glue maps the name for a model
# estimated inside R, but a file written under the R name left BayesTS without a
# starting precision, and it refused the file. 'direction' is "write" for the
# name in the file and "read" for the name in R.
.hdf5_initial_name <- function(model, name, direction) {
  constant_gamma <- isTRUE(model[["error"]] %in% c("gamma", "gamma+covar")) &&
    !isTRUE(model[["tvp"]])
  if (!constant_gamma) {
    return(name)
  }
  if (direction == "write" && name == "u_omega_inv") {
    return("u_sigma_inv")
  }
  if (direction == "read" && name == "u_sigma_inv") {
    return("u_omega_inv")
  }
  name
}
