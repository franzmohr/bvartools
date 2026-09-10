#' Export to HDF5 File
#'
#' Exports the content of an object of class 'bvarmodel' to an HDF5 file.
#'
#' @param object list of class 'bvarmodel'.
#' @param filename path to the file, in which output should be stored.
#' @param group the group the model's tree should hang under inside its file.
#' Defaults to \code{""}, the root of the file, which is where a file holding a
#' single model puts it. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The path to the written file, invisibly.
#'
#' @details
#'
#' With a \code{group} every path is written under it instead of at the root,
#' so one file can hold several models side by side. The spelling is the one the
#' BayesTS command line uses for its \code{--group} flag: a leading slash and no
#' trailing slash, with \code{""} for the root. Intermediate groups are created
#' as needed, and \code{\link{list_models_in_hdf5}} reports which groups of a
#' file hold a model.
#'
#' What must not already be there is the model rather than the file. Without a
#' \code{group} the model is the whole file, so an existing file is refused;
#' with one, only that group has to be free, which is what lets a second model
#' be added beside the first.
#'
#' A write that cannot be completed raises the error rather than absorbing it,
#' and undoes what the call created: a file it made is removed, a group it added
#' to an existing file is unlinked on its own so the models beside it survive.
#' A retry then meets the original problem rather than the leftovers.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' train <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(data = train,
#'                           p = 1,
#'                           deterministic = "const",
#'                           iterations = 10,
#'                           burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1),
#'                     sigma = list(df = 3, scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Save model
#' path_to_model <- tempfile(fileext = ".h5")
#' write_to_hdf5(model, filename = path_to_model)
#' 
#' @export
write_to_hdf5.bvarmodel <- function(object, filename, group = "", ...) {
  
  group <- .normalize_hdf5_group(group)
  
  # Check if filename is valid
  if (dir.exists(filename)) {
    stop("Argument 'filename' is not a path to a file.")
  }
  
  # What must not already be there is the model, not the file. Without a group a
  # model is the whole file, so an existing file is still refused; with one,
  # a file that already holds other models is exactly what is being added to,
  # and only that group has to be free.
  file_existed <- file.exists(filename)
  if (group == "" && file_existed) {
    stop(paste0("File ", filename, " already exists."))
  }
  
  # Create or open an HDF5 file
  h5_file <- hdf5r::h5file(filename, mode = "a")
  
  # The body below used to sit inside a try() that discarded its result: any
  # failure -- a full disk, an unwritable path, a malformed element of 'object'
  # -- was swallowed and the function returned as though it had worked, leaving
  # a half-written file that looked finished. The error reaches the caller now.
  #
  # What has to be undone on the way out of a failure depends on what this call
  # created. A file it made is removed whole; a group it added to a file that
  # was already there is unlinked on its own, so the models beside it survive.
  # Either way the handle is closed, or the file stays locked for the rest of
  # the session. HDF5 does not reclaim the space of an unlinked group, but the
  # name is free again, which is what a retry needs.
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

  
  
  # Model information ----
  group_model <- .hdf5_group(handles, output, "model")

  # Save each available element in 'model'
  # Also useful if user adds additional elements in 'model' manually
  #
  # What the group already carries is read once instead of once per element:
  # the only thing that changes it is this loop, which keeps track itself.
  attrs_model <- hdf5r::h5attr_names(group_model)
  for (i in names(object[["model"]])) {
    if (!i %in% attrs_model) {
      hdf5r::h5attr(group_model, i) <- object[["model"]][[i]]
      attrs_model <- c(attrs_model, i)
    }
  }
  hdf5r::h5attr(group_model, "rclass") <- class(object)

  # Data ----
  group_data <- .hdf5_group(handles, output, "data")

  ## Original ----
  group_data_original <- .hdf5_group(handles, group_data, "original")
  for (i in c("endogen", "exogen", "deterministic")) {
    series <- object[["data"]][["original"]][[i]]
    if (!is.null(series)) {
      .hdf5_write(group_data_original, i, series, .hdf5_series_attrs(series))
    }
  }

  ## Train ----
  group_data_train <- .hdf5_group(handles, group_data, "train")
  for (i in c("y", "x")) {
    series <- object[["data"]][["train"]][[i]]
    if (!is.null(series)) {
      .hdf5_write(group_data_train, i, series, .hdf5_series_attrs(series))
    }
  }
  # Data without time series information
  for (i in "z") {
    if (!is.null(object[["data"]][["train"]][[i]])) {
      .hdf5_write(group_data_train, i, object[["data"]][["train"]][[i]])
    }
  }

  ## Forecast input ----
  if (!is.null(object[["data"]][["forecast"]][["z"]])) {
    group_data_forecast <- .hdf5_group(handles, group_data, "forecast")
    .hdf5_write(group_data_forecast, "z", object[["data"]][["forecast"]][["z"]])
  }

  # Priors ----
  # Nothing is written at all for a model that has not been through
  # add_priors(), because an empty group is worse than no group: the reader
  # rebuilds 'priors' from the names the file actually holds, so a group that
  # is present but empty comes back as an empty list rather than as absent.
  if (length(object[["priors"]]) > 0) {
    group_priors <- .hdf5_group(handles, output, "priors")

    ## Those kept in a group of their own ----
    for (i in c("a", "psi")) {
      if (!is.null(object[["priors"]][[i]])) {
        group_prior <- .hdf5_group(handles, group_priors, i)
        for (j in names(object[["priors"]][[i]])) {
          .hdf5_write(group_prior, j, object[["priors"]][[i]][[j]])
        }
      }
    }

    ## u_sigma_inv ----
    # Which hyperparameters there are is decided by the error specification.
    u_sigma_priors <- switch(object[["model"]][["error"]],
                             "wishart" = c("df", "scale"),
                             "gamma" = ,
                             "gamma+covar" = c("shape", "rate"),
                             "sv" = ,
                             "sv+covar" = c("mu", "v_inv", "shape", "rate", "sigma", "offset"),
                             stop("Error specification not implemented"))
    # Created only once there is something to put in it, for the same reason the
    # priors group is. Single-bracket indexing is what makes the subset safe when
    # a name is absent, where [[ would throw instead of giving NULL.
    u_sigma_values <- object[["priors"]][["u_sigma"]][u_sigma_priors]
    if (any(!vapply(u_sigma_values, is.null, logical(1)))) {
      group_priors_u_sigma <- .hdf5_group(handles, group_priors, "u_sigma")
      for (i in seq_along(u_sigma_priors)) {
        .hdf5_write(group_priors_u_sigma, u_sigma_priors[i], u_sigma_values[[i]])
      }
    }
  }

  # Initial values ----
  if (!is.null(object[["initial"]])) {
    group_initial <- .hdf5_group(handles, output, "initial")
    for (i in names(object[["initial"]])) {
      .hdf5_write(group_initial, i, object[["initial"]][[i]])
    }
  }

  # Posterior ----
  if (!is.null(object[["posterior"]])) {
    group_posterior <- .hdf5_group(handles, output, "posterior")

    ## Draws kept in a group of their own ----
    #
    # One shape for all of them. These used to be a block of the same code
    # each, and one of the blocks attached its attributes to the datasets of
    # the block above it rather than to its own.
    for (i in c("a", "psi", "u_omega_inv", "u_sigma_inv")) {
      if (i %in% names(object[["posterior"]])) {
        group_draws <- .hdf5_group(handles, group_posterior, i)
        for (j in names(object[["posterior"]][[i]])) {
          draws <- object[["posterior"]][[i]][[j]]
          .hdf5_write(group_draws, j, draws, .hdf5_draws_attrs(draws))
        }
      }
    }

    ## Draws kept on their own ----
    for (i in c("loglik", "forecast", "forecast_error")) {
      if (i %in% names(object[["posterior"]])) {
        draws <- object[["posterior"]][[i]]
        .hdf5_write(group_posterior, i, draws, .hdf5_draws_attrs(draws))
      }
    }
  }

  # Close file
  .hdf5_close(handles, h5_file)
  completed <- TRUE

  invisible(filename)
}
