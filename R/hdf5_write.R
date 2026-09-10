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

  dataset <- group$create_dataset(name, value)

  for (i in names(attrs)) {
    hdf5r::h5attr(dataset, i) <- attrs[[i]]
  }

  dataset$close()

  return(invisible(NULL))
}

# The attributes that make a stored matrix a time series again on the way back.
#
# The error correction term of a model that has been through
# scale_error_correction() carries the factors it was divided by, and they are
# the only way back to the scale of the data. Without them an exported model
# that was scaled could never be interpreted again, so they travel with it. The
# names are not stored: they are the variable names, which are written anyway.
.hdf5_series_attrs <- function(x) {

  result <- list("variables" = dimnames(x)[[2]], "tsp" = stats::tsp(x))

  if (!is.null(attr(x, "scale"))) {
    result[["scale"]] <- unname(attr(x, "scale"))
  }

  return(result)
}

# The attributes that make stored draws a chain again on the way back.
.hdf5_draws_attrs <- function(x) {
  mcpar <- coda::mcpar(x)
  list("start" = mcpar[1], "end" = mcpar[2], "thin" = mcpar[3])
}
