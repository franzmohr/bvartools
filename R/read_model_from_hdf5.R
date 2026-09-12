#' Import Models from HDF5 Files
#' 
#' Imports model information and posterior draws from an HDF5 file.
#' 
#' @param filename Path to an HDF5 file containing model data.
#' @param group the group the model's tree hangs under inside its file.
#' Defaults to \code{""}, the root of the file, which is where a file holding a
#' single model puts it. See 'Details'.
#' 
#' @details
#' 
#' With a \code{group} every path the reader looks for is read under it
#' instead of at the root, so one file can hold several models side by side.
#' \code{\link{list_models_in_hdf5}} reports which groups of a file hold one.
#' 
#' The spelling is the one the BayesTS command line uses for its \code{--group}
#' flag: a leading slash and no trailing slash, with \code{""} for the root.
#' 
#' @export
read_model_from_hdf5 <- function(filename, group = "") {
  
  group <- .normalize_hdf5_group(group)
  
  h5_file <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (h5_file$is_valid) h5_file$close_all(), add = TRUE)
  
  # Every path below is named against this rather than against the file, which
  # is all a group amounts to on the way in.
  h5_root <- .hdf5_model_root(h5_file, group)
  
  h5_names <- names(h5_root)
  
  result <- NULL
  
  if ("model" %in% h5_names) {
    result[["model"]] <- hdf5r::h5attributes(h5_root[["model"]])

    # Written as a group of its own, because the restriction table is a table
    # and the rest of a specification is a set of values. See
    # write_to_hdf5.bvarmodel.
    if ("sign_restrictions" %in% names(h5_root[["model"]])) {
      group_sign <- h5_root[["model"]][["sign_restrictions"]]
      dataset <- group_sign[["restrictions"]]
      restrictions <- as.matrix(hdf5r::readDataSet(dataset))
      colnames(restrictions) <- hdf5r::h5attr(dataset, "columns")

      period <- NULL
      if ("period" %in% hdf5r::h5attr_names(group_sign)) {
        period <- hdf5r::h5attr(group_sign, "period")
      }

      result[["model"]][["sign_restrictions"]] <- list(
        "restrictions" = as.data.frame(restrictions),
        "max_tries" = hdf5r::h5attr(group_sign, "max_tries"),
        "period" = period
      )
    }
  } else {
    stop("File ", filename, " does not contain model specification",
         if (group != "") paste0(" in group ", group), ".")
  }
  
  if ("data" %in% h5_names) {
    
    result[["data"]] <- list()
    
    if ("original" %in% names(h5_root[["data"]])) {
      result[["data"]][["original"]] <- list()
      for (i in c("endogen", "exogen", "deterministic")) {
        if (i %in% names(h5_root[["data"]][["original"]])) {
          result[["data"]][["original"]][[i]] <- stats::ts(as.matrix(hdf5r::readDataSet(h5_root[["data"]][["original"]][[i]])), class = c("mts", "ts", "matrix"))
          dimnames(result[["data"]][["original"]][[i]]) <- list(NULL, hdf5r::h5attr(h5_root[["data"]][["original"]][[i]], "variables"))
          stats::tsp(result[["data"]][["original"]][[i]]) <- hdf5r::h5attr(h5_root[["data"]][["original"]][[i]], "tsp")
        }
      }
    }
    
    if ("train" %in% names(h5_root[["data"]])) {
      result[["data"]][["train"]] <- list()
      for (i in c("y", "w", "x")) {
        if (i %in% names(h5_root[["data"]][["train"]])) {
          dataset <- h5_root[["data"]][["train"]][[i]]
          variables <- hdf5r::h5attr(dataset, "variables")
          result[["data"]][["train"]][[i]] <- stats::ts(as.matrix(hdf5r::readDataSet(dataset)))
          dimnames(result[["data"]][["train"]][[i]]) <- list(NULL, variables)
          stats::tsp(result[["data"]][["train"]][[i]]) <- hdf5r::h5attr(dataset, "tsp")

          # A model that was exported while its error correction term was
          # scaled carries the factors it was divided by. They are named after
          # the variables, which is why the names are not stored separately.
          if ("scale" %in% hdf5r::h5attr_names(dataset)) {
            factors <- hdf5r::h5attr(dataset, "scale")
            names(factors) <- variables
            attr(result[["data"]][["train"]][[i]], "scale") <- factors
          }
        }
      }
      for (i in c("z")) {
        if (i %in% names(h5_root[["data"]][["train"]])) {
          result[["data"]][["train"]][[i]] <- as.matrix(hdf5r::readDataSet(h5_root[["data"]][["train"]][[i]]))
        }
      }
    }
    
    if ("forecast" %in% names(h5_root[["data"]])) {
      result[["data"]][["forecast"]] <- list()
      # `x` is the compact layout, `z` the SUR one written before it. Both are
      # read under their own name; the C++ side takes either.
      if ("x" %in% names(h5_root[["data"]][["forecast"]])) {
        result[["data"]][["forecast"]][["x"]] <- as.matrix(hdf5r::readDataSet(h5_root[["data"]][["forecast"]][["x"]]))
      }
      if ("z" %in% names(h5_root[["data"]][["forecast"]])) {
        result[["data"]][["forecast"]][["z"]] <- as.matrix(hdf5r::readDataSet(h5_root[["data"]][["forecast"]][["z"]]))
      }
    }
  }
  
  
  if ("priors" %in% h5_names) {
    
    result[["priors"]] <- list()
    
    for (i in names(h5_root[["priors"]])) {
      result[["priors"]][[i]] <- list()
      for (j in names(h5_root[["priors"]][[i]])) {
        result[["priors"]][[i]][[j]] <- as.matrix(hdf5r::readDataSet(h5_root[["priors"]][[i]][[j]]))
      }
    }
  }
  
  
  if ("initial" %in% h5_names) {
    
    result[["initial"]] <- list()
    
    for (i in names(h5_root[["initial"]])) {
      result[["initial"]][[i]] <- as.matrix(hdf5r::readDataSet(h5_root[["initial"]][[i]]))
    }
  }
  
  
  if ("posterior" %in% h5_names) {
    
    result[["posterior"]] <- list()
    
    for (i in names(h5_root[["posterior"]])) {
      element <- h5_root[["posterior"]][[i]]
      if (inherits(element, "H5Group")) {
        # Whatever draws the group holds, rather than the names they were
        # expected to have. Each is written with its own start, end and
        # thinning interval, so a name left out of a list here is a dataset
        # that is in the file and read back as nothing. `sigma`, which every
        # time varying model keeps beside its coefficients, was exactly that.
        for (j in names(element)) {
          dataset <- element[[j]]
          result[["posterior"]][[i]][[j]] <- coda::mcmc(as.matrix(hdf5r::readDataSet(dataset)),
                                                        start = hdf5r::h5attr(dataset, "start"),
                                                        end = hdf5r::h5attr(dataset, "end"),
                                                        thin = hdf5r::h5attr(dataset, "thin")) 
        }
      } else {
        # The draws the writer keeps on their own rather than in a group,
        # loglik and forecast among them. Told apart from a group by what they
        # are rather than by name, so that another one needs nothing here.
        result[["posterior"]][[i]] <- coda::mcmc(hdf5r::readDataSet(element))
      }
    }
  }
  
  
  h5_file$close_all()
  
  # Determine the class of the returned object based on the used algorithm
  result_class <- result[["model"]][["rclass"]]
  if (is.null(result_class)) {
    bvarmodel <- c("VarNormalAld", "VarNormalGamma", "VarNormalStochvol",
                   "VarNormalWishart", "VarTvpAld", "VarTvpGamma",
                   "VarTvpStochvol", "VarTvpWishart")
    if (result[["model"]][["algorithm"]] %in% bvarmodel) {
      result_class <- c("bvarmodel", "list")
    }
    
    bvecmodel <- c("VecNormalGamma", "VecNormalStochvol", "VecNormalWishart",
                   "VecTvpGamma", "VecTvpStochvol", "VecTvpWishart",
                   "VecKlgs2010")
    if (result[["model"]][["algorithm"]] %in% bvecmodel) {
      result_class <- c("bvecmodel", "list")
    } 
  }
  
  class(result) <- result_class
  
  return(result)
}