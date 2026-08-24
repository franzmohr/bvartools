#' Export to HDF5 File
#'
#' Exports the content of an object of class 'bvarmodel' to an HDF5 file.
#'
#' @param object list of class 'bvarmodel'.
#' @param filename path to the file, in which output should be stored.
#' @param ... further arguments passed to or from other methods.
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
write_to_hdf5.bvarmodel <- function(object, filename, ...) {
  
  # Check if filename is valid
  if (dir.exists(filename)) {
    stop("Argument 'filename' is not a path to a file.")
  }
  
  # Check if output file already exists
  if (file.exists(filename)) {
    stop(paste0("File ", filename, " already exists."))
  }
  
  try({
    
    # Create or open an HDF5 file
    output <- hdf5r::h5file(filename, mode = "a")
    
    # Model information ----
    if (!"model" %in% names(output)) {
      group_model <- output$create_group("model")
    } else {
      group_model <- output[["model"]]
    }
    
    # Save each available element in 'model'
    # Also useful if user adds additional elements in 'model' manually
    for (i in names(object[["model"]])) {
      if (!i %in% hdf5r::h5attr_names(group_model)) {
        hdf5r::h5attr(group_model, i) <- object[["model"]][[i]] 
      }
    }
    hdf5r::h5attr(group_model, "rclass") <- class(object)
    
    # Data ----
    if (!"data" %in% names(output)) {
      group_data <- output$create_group("data")
    } else {
      group_data <- output[["data"]]
    }
    
    ## Original ----
    if (!"original" %in% names(group_data)) {
      group_data_original <- group_data$create_group("original")
    } else {
      group_data_original <- group_data[["original"]]
    }
    for (i in c("endogen", "exogen", "deterministic")) {
      if (!is.null(object[["data"]][["original"]][[i]]) & !i %in% names(group_data_original)) {
        group_data_original[[i]] <- object[["data"]][["original"]][[i]]
        hdf5r::h5attr(group_data_original[[i]], "variables") <- dimnames(object[["data"]][["original"]][[i]])[[2]]
        hdf5r::h5attr(group_data_original[[i]], "tsp") <- stats::tsp(object[["data"]][["original"]][[i]])
      } 
    }
    
    ## Train ----
    if (!"train" %in% names(group_data)) {
      group_data_train <- group_data$create_group("train")
    } else {
      group_data_train <- group_data[["train"]]
    }
    for (i in c("y", "x")) {
      if (!is.null(object[["data"]][["train"]][[i]]) & !i %in% names(group_data_train)) {
        group_data_train[[i]] <- object[["data"]][["train"]][[i]]
        hdf5r::h5attr(group_data_train[[i]], "variables") <- dimnames(object[["data"]][["train"]][[i]])[[2]]
        hdf5r::h5attr(group_data_train[[i]], "tsp") <- stats::tsp(object[["data"]][["train"]][[i]])
      }  
    }
    # Data without time series information
    for (i in "z") {
      if (!is.null(object[["data"]][["train"]][[i]]) & !i %in% names(group_data_train)) {
        group_data_train[[i]] <- object[["data"]][["train"]][[i]]
      } 
    }
    
    ## Forecast input ----
    if (!is.null(object[["data"]][["forecast"]][["z"]])) {
      if (!"forecast" %in% names(group_data)) {
        group_data_forecast <- group_data$create_group("forecast")
      } else {
        group_data_forecast <- group_data[["forecast"]]
      }
      group_data_forecast[["z"]] <- object[["data"]][["forecast"]][["z"]]
    }
    
    
    
    # Priors ----
    if (!"priors" %in% names(output)) {
      group_priors <- output$create_group("priors")
    } else {
      group_priors <- output[["priors"]]
    }
    
    ## a ----
    if (!is.null(object[["priors"]][["a"]])) {
      if (!"a" %in% names(group_priors)) {
        group_priors_a <- group_priors$create_group("a")
      } else {
        group_priors_a <- group_priors[["a"]]
      }
      for (i in names(object[["priors"]][["a"]])) {
        if (!i %in% names(group_priors_a)) {
          group_priors_a[[i]] <- object[["priors"]][["a"]][[i]]
        }
      }
    }
    
    # psi ----
    if (!is.null(object[["priors"]][["psi"]])) {
      if (!"psi" %in% names(group_priors)) {
        group_priors_psi <- group_priors$create_group("psi")
      } else {
        group_priors_psi <- group_priors[["psi"]]
      }
      for (i in names(object[["priors"]][["psi"]])) {
        if (!i %in% names(group_priors_psi)) {
          group_priors_psi[[i]] <- object[["priors"]][["psi"]][[i]]
        }
      }
    }
    
    
    ## u_sigma_inv ----
    if (!"u_sigma" %in% names(group_priors)) {
      group_priors_u_sigma <- group_priors$create_group("u_sigma")
    } else {
      group_priors_u_sigma <- group_priors[["u_sigma"]]
    }
    if (!object[["model"]][["error"]] %in% c("wishart", "gamma", "gamma+covar", "sv", "sv+covar")) {
      stop("Error specification not implemented")
    }
    if (object[["model"]][["error"]] == "wishart") {
      for (i in c("df", "scale")) {
        if (!i %in% names(group_priors_u_sigma)) {
          group_priors_u_sigma[[i]] <- object[["priors"]][["u_sigma"]][[i]]
        }
      }
    }
    if (object[["model"]][["error"]] %in% c("gamma", "gamma+covar")) {
      for (i in c("shape", "rate")) {
        if (!i %in% names(group_priors_u_sigma)) {
          group_priors_u_sigma[[i]] <- object[["priors"]][["u_sigma"]][[i]]
        }
      }
    }
    if (object[["model"]][["error"]] %in% c("sv", "sv+covar")) {
      for (i in c("mu", "v_inv", "shape", "rate", "sigma", "offset")) {
        if (!i %in% names(group_priors_u_sigma)) {
          group_priors_u_sigma[[i]] <- object[["priors"]][["u_sigma"]][[i]]
        }
      }
    }
    
    
    # Initial values ----
    if (!is.null(object[["initial"]])) {
      if (!"initial" %in% names(output)) {
        group_initial <- output$create_group("initial")
      } else {
        group_initial <- output[["initial"]]
      }
      for (i in names(object[["initial"]])) {
        if (!i %in% names(group_initial)) {
          group_initial[[i]] <- object[["initial"]][[i]] 
        }
      } 
    }
    
    # Posterior ----
    if (!is.null(object[["posterior"]])) {
      if (!"posterior" %in% names(output)) {
        group_posterior <- output$create_group("posterior")
      } else {
        group_posterior <- output[["posterior"]]
      }
      
      if ("a" %in% names(object[["posterior"]])) {
        if (!"a" %in% names(group_posterior)) {
          group_posterior_a <- group_posterior$create_group("a")
        } else {
          group_posterior_a <- group_posterior[["a"]]
        }
        for (i in names(object[["posterior"]][["a"]])) {
          if (!i %in% names(group_posterior_a)) {
            group_posterior_a[[i]] <- object[["posterior"]][["a"]][[i]] 
            mcpar_temp <- coda::mcpar(object[["posterior"]][["a"]][[i]])
            hdf5r::h5attr(group_posterior_a[[i]], "start") <- mcpar_temp[1]
            hdf5r::h5attr(group_posterior_a[[i]], "end") <- mcpar_temp[2]
            hdf5r::h5attr(group_posterior_a[[i]], "thin") <- mcpar_temp[3]
          }
        }  
      }
      
      if ("psi" %in% names(object[["posterior"]])) {
        if (!"psi" %in% names(group_posterior)) {
          group_posterior_psi <- group_posterior$create_group("psi")
        } else {
          group_posterior_psi <- group_posterior[["psi"]]
        }
        for (i in names(object[["posterior"]][["psi"]])) {
          if (!i %in% names(group_posterior_psi)) {
            group_posterior_psi[[i]] <- object[["posterior"]][["psi"]][[i]] 
            mcpar_temp <- coda::mcpar(object[["posterior"]][["psi"]][[i]])
            hdf5r::h5attr(group_posterior_a[[i]], "start") <- mcpar_temp[1]
            hdf5r::h5attr(group_posterior_a[[i]], "end") <- mcpar_temp[2]
            hdf5r::h5attr(group_posterior_a[[i]], "thin") <- mcpar_temp[3]
          }
        }  
      }
      
      if ("u_omega_inv" %in% names(object[["posterior"]])) {
        if (!"u_omega_inv" %in% names(group_posterior)) {
          group_posterior_u_omega_inv <- group_posterior$create_group("u_omega_inv")
        } else {
          group_posterior_u_omega_inv <- group_posterior[["u_omega_inv"]]
        }
        for (i in names(object[["posterior"]][["u_omega_inv"]])) {
          if (!i %in% names(group_posterior_u_omega_inv)) {
            group_posterior_u_omega_inv[[i]] <- object[["posterior"]][["u_omega_inv"]][[i]]
            mcpar_temp <- coda::mcpar(object[["posterior"]][["u_omega_inv"]][[i]])
            hdf5r::h5attr(group_posterior_u_omega_inv[[i]], "start") <- mcpar_temp[1]
            hdf5r::h5attr(group_posterior_u_omega_inv[[i]], "end") <- mcpar_temp[2]
            hdf5r::h5attr(group_posterior_u_omega_inv[[i]], "thin") <- mcpar_temp[3]
          }
        }  
      }
      
      if ("u_sigma_inv" %in% names(object[["posterior"]])) {
        if (!"u_sigma_inv" %in% names(group_posterior)) {
          group_posterior_u_sigma_inv <- group_posterior$create_group("u_sigma_inv")
        } else {
          group_posterior_u_sigma_inv <- group_posterior[["u_sigma_inv"]]
        }
        for (i in names(object[["posterior"]][["u_sigma_inv"]])) {
          if (!i %in% names(group_posterior_u_sigma_inv)) {
            group_posterior_u_sigma_inv[[i]] <- object[["posterior"]][["u_sigma_inv"]][[i]]
            mcpar_temp <- coda::mcpar(object[["posterior"]][["u_sigma_inv"]][[i]])
            hdf5r::h5attr(group_posterior_u_sigma_inv[[i]], "start") <- mcpar_temp[1]
            hdf5r::h5attr(group_posterior_u_sigma_inv[[i]], "end") <- mcpar_temp[2]
            hdf5r::h5attr(group_posterior_u_sigma_inv[[i]], "thin") <- mcpar_temp[3]
          }
        }  
      }
      
      if ("loglik" %in% names(object[["posterior"]])) {
        if (!"loglik" %in% names(group_posterior)) {
          group_posterior[["loglik"]] <- object[["posterior"]][["loglik"]]
          mcpar_temp <- coda::mcpar(object[["posterior"]][["loglik"]])
          hdf5r::h5attr(group_posterior[["loglik"]], "start") <- mcpar_temp[1]
          hdf5r::h5attr(group_posterior[["loglik"]], "end") <- mcpar_temp[2]
          hdf5r::h5attr(group_posterior[["loglik"]], "thin") <- mcpar_temp[3]
        }
      }
      
      if ("forecast" %in% names(object[["posterior"]])) {
        if (!"forecast" %in% names(group_posterior)) {
          group_posterior[["forecast"]] <- object[["posterior"]][["forecast"]]
          mcpar_temp <- coda::mcpar(object[["posterior"]][["forecast"]])
          hdf5r::h5attr(group_posterior[["forecast"]], "start") <- mcpar_temp[1]
          hdf5r::h5attr(group_posterior[["forecast"]], "end") <- mcpar_temp[2]
          hdf5r::h5attr(group_posterior[["forecast"]], "thin") <- mcpar_temp[3]
        }
      }
      
      if ("forecast_error" %in% names(object[["posterior"]])) {
        if (!"forecast_error" %in% names(group_posterior)) {
          group_posterior[["forecast_error"]] <- object[["posterior"]][["forecast_error"]]
          mcpar_temp <- coda::mcpar(object[["posterior"]][["forecast_error"]])
          hdf5r::h5attr(group_posterior[["forecast_error"]], "start") <- mcpar_temp[1]
          hdf5r::h5attr(group_posterior[["forecast_error"]], "end") <- mcpar_temp[2]
          hdf5r::h5attr(group_posterior[["forecast_error"]], "thin") <- mcpar_temp[3]
        }
      }
    }
    
    
    # Close file
    output$close_all()
  })
}