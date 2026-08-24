#' Export to HDF5 File
#' 
#' Exports the content of an object of class 'modellist' to multiple HDF5 file.
#' 
#' @param object an object of class 'modellist'.
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
#'                           p = 0:4,
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
#' path_to_target_directory <- tempdir()
#' write_to_hdf5(model, folder = path_to_target_directory)
#' 
#' 
#' @export
write_to_hdf5.modellist <- function(object, folder, ...){
  
  if (dir.exists(folder)) {
    
    for (i in 1:length(object)) {
      
      if (any(class(object[[i]]) %in% c("bvarmodel", "bvecmodel"))) {
        
        # Get model type from the first object in the list
        model_type <- object[[i]][["model"]][["algorithm"]]
        
        p <- sprintf(paste0("%0", 2, "d"), object[[i]][["model"]][["p"]])
        model_details <- paste0("p=", p)
        if (object[[i]][["model"]][["m"]] > 0) {
          s <- sprintf(paste0("%0", 2, "d"), object[[i]][["model"]][["s"]])
          model_details <- paste0(model_details, "-s=", s)
        }
        model_details <- paste0(model_details, "-n=", object[[i]][["model"]][["n"]])
        if (!is.null(object[[i]][["model"]][["varsel"]])) {
          varsel_type <- paste0("-varsel=", object[[i]][["model"]][["varsel"]])
        } else {
          varsel_type <- NULL
        }
        
        # Get list of files in the folder
        file_names <- list.files(folder, recursive = FALSE, full.names = FALSE)
        file_names <- file_names[grepl(".h5", file_names)]
        
        iter <- 0
        cond <- TRUE
        while(cond) {
          iter <- iter + 1
          candidate_name <- paste0(model_type,
                                   varsel_type,
                                   "-", model_details,
                                   "-", sprintf(paste0("%0", 3, "d"), iter), ".h5")
          cond <- candidate_name %in% file_names
        }
        file_name <- file.path(folder, candidate_name)
        
        write_to_hdf5(object[[i]], filename = file_name, ...)
        
      } else {
        write_to_hdf5(object[[i]], folder = folder, ...)
      }
    }
    
  } else {
    stop("Specified 'folder' does not exist.")
  }
  
}