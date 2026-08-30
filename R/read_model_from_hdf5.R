#' Import Models from HDF5 Files
#' 
#' Imports model information and posterior draws from an HDF5 file.
#' 
#' @param filename Path to an HDF5 file containing model data.
#' 
#' @export
read_model_from_hdf5 <- function(filename) {
  
  h5_file <- hdf5r::h5file(filename, mode = "r")
  
  h5_names <- names(h5_file)
  
  result <- NULL
  
  if ("model" %in% h5_names) {
    result[["model"]] <- hdf5r::h5attributes(h5_file[["model"]])
  } else {
    stop("File ", filename, " does not contain model specification.")
  }
  
  if ("data" %in% h5_names) {
    
    result[["data"]] <- list()
    
    if ("original" %in% names(h5_file[["data"]])) {
      result[["data"]][["original"]] <- list()
      for (i in c("endogen", "exogen", "deterministic")) {
        if (i %in% names(h5_file[["data"]][["original"]])) {
          result[["data"]][["original"]][[i]] <- stats::ts(as.matrix(hdf5r::readDataSet(h5_file[["data"]][["original"]][[i]])), class = c("mts", "ts", "matrix"))
          dimnames(result[["data"]][["original"]][[i]]) <- list(NULL, hdf5r::h5attr(h5_file[["data"]][["original"]][[i]], "variables"))
          stats::tsp(result[["data"]][["original"]][[i]]) <- hdf5r::h5attr(h5_file[["data"]][["original"]][[i]], "tsp")
        }
      }
    }
    
    if ("train" %in% names(h5_file[["data"]])) {
      result[["data"]][["train"]] <- list()
      for (i in c("y", "w", "x")) {
        if (i %in% names(h5_file[["data"]][["train"]])) {
          result[["data"]][["train"]][[i]] <- stats::ts(as.matrix(hdf5r::readDataSet(h5_file[["data"]][["train"]][[i]])))
          dimnames(result[["data"]][["train"]][[i]]) <- list(NULL, hdf5r::h5attr(h5_file[["data"]][["train"]][[i]], "variables"))
          stats::tsp(result[["data"]][["train"]][[i]]) <- hdf5r::h5attr(h5_file[["data"]][["train"]][[i]], "tsp")
        }
      }
      for (i in c("z")) {
        if (i %in% names(h5_file[["data"]][["train"]])) {
          result[["data"]][["train"]][[i]] <- as.matrix(hdf5r::readDataSet(h5_file[["data"]][["train"]][[i]]))
        }
      }
    }
    
    if ("forecast" %in% names(h5_file[["data"]])) {
      result[["data"]][["forecast"]] <- list()
      if ("z" %in% names(h5_file[["data"]][["forecast"]])) {
        result[["data"]][["forecast"]][["z"]] <- as.matrix(hdf5r::readDataSet(h5_file[["data"]][["forecast"]][["z"]]))
      }
    }
  }
  
  
  if ("priors" %in% h5_names) {
    
    result[["priors"]] <- list()
    
    for (i in names(h5_file[["priors"]])) {
      result[["priors"]][[i]] <- list()
      for (j in names(h5_file[["priors"]][[i]])) {
        result[["priors"]][[i]][[j]] <- as.matrix(hdf5r::readDataSet(h5_file[["priors"]][[i]][[j]]))
      }
    }
  }
  
  
  if ("initial" %in% h5_names) {
    
    result[["initial"]] <- list()
    
    for (i in names(h5_file[["initial"]])) {
      result[["initial"]][[i]] <- as.matrix(hdf5r::readDataSet(h5_file[["initial"]][[i]]))
    }
  }
  
  
  if ("posterior" %in% h5_names) {
    
    result[["posterior"]] <- list()
    
    for (i in names(h5_file[["posterior"]])) {
      if (i %in% c("loglik", "forecast")) {
        result[["posterior"]][[i]] <- coda::mcmc(hdf5r::readDataSet(h5_file[["posterior"]][[i]]))
      } else {
        for (j in c("coeffs", "lambda")) {
          if (j %in% names(h5_file[["posterior"]][[i]])) {
            result[["posterior"]][[i]][[j]] <- coda::mcmc(as.matrix(hdf5r::readDataSet(h5_file[["posterior"]][[i]][[j]])),
                                                          start = hdf5r::h5attr(h5_file[["posterior"]][[i]][[j]], "start"),
                                                          end = hdf5r::h5attr(h5_file[["posterior"]][[i]][[j]], "end"),
                                                          thin = hdf5r::h5attr(h5_file[["posterior"]][[i]][[j]], "thin")) 
          }
        }
      }
    }
  }
  
  
  hdf5r::h5close(h5_file)
  
  # Determine the class of the returned object based on the used algorithm
  result_class <- result[["model"]][["rclass"]]
  if (is.null(result_class)) {
    bvarmodel <- c("VarNormalGamma", "VarNormalStochvol", "VarNormalWishart",
                   "VarTvpGamma", "VarTvpStochvol", "VarTvpWishart")
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