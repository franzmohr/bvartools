#' Export to HDF5 File
#' 
#' Exports the content of an object of class 'expandingwindow' to multiple HDF5 file.
#' 
#' @param object an object of class 'expandingwindow'.
#' @param folder path to the folder where the individual models in argument
#' \code{object} should be saved.
#' @param ... further arguments passed to or from other methods.
#' 
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' train <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(data = train,
#'                           p = 0:2,
#'                           deterministic = "const",
#'                           iterations = 10,
#'                           burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Use expanding window
#' model <- use_expanding_window(model, start = 1982.5)
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
#' path_to_target_directory <- tempdir()
#' write_to_hdf5(model, folder = path_to_target_directory)
#' 
#' 
#' @export
write_to_hdf5.expandingwindow <- function(object, folder, ...){
  
  if (dir.exists(folder)) {
    
    # Get model type from the first object in the list
    model_type <- object[[1]][["model"]][["algorithm"]]
    
    # Get list of files in the folder
    dir_names <- list.dirs(folder, recursive = FALSE, full.names = FALSE)
    
    model_name <- paste0("01-ExpWind-", model_type)
    model_name <- paste0(model_name, paste0("-varsel=", object[[1]][["model"]][["varsel"]]))
    model_name <- paste0(model_name, paste0("-n=", sprintf(paste0("%0", 2, "d"), object[[1]][["model"]][["n"]])))
    model_name <- paste0(model_name, paste0("-p=", sprintf(paste0("%0", 2, "d"), object[[1]][["model"]][["p"]])))
    model_name <- paste0(model_name, paste0("-m=", sprintf(paste0("%0", 2, "d"), object[[1]][["model"]][["m"]])))
    model_name <- paste0(model_name, paste0("-s=", sprintf(paste0("%0", 2, "d"), object[[1]][["model"]][["s"]])))
    
    i <- 1
    cond <- dir.exists(file.path(folder, model_name))
    if (cond) {
      while (cond) {
        substring(model_name, 1, 2) <- sprintf(paste0("%0", 2, "d"), i)
        cond <- dir.exists(file.path(folder, model_name))
        i <- i + 1
      }
    }
    
    folder <- file.path(folder, model_name)
    dir.create(folder)
    
    # Add leading zeros to file names, so that they are listed in correct order
    # in the filesystem
    n_models <- length(object)
    num_digits <- nchar(as.character(n_models))
    model_id <- sprintf(paste0("%0", num_digits, "d"), 1:n_models)
    
    for (i in 1:n_models) {
        write_to_hdf5(object[[i]], filename = file.path(folder, paste0("Window-", model_id[i], ".h5")), ...)
    }
    
  } else {
    stop("Specified 'folder' does not exist.")
  }
  
}