#' Import Models from HDF5 Files
#' 
#' Imports model information and posterior draws from an HDF5 file.
#' 
#' @param folder Path to a folder containing model data.
#' 
#' @export
read_models_from_folder <- function(folder) {
  
  if (!dir.exists(folder)) {
    stop("Specified folder does not exist.")
  }
  
  all_files <- list.files(folder, recursive = TRUE, full.names = TRUE)
  all_files <- all_files[grepl(".h5", all_files)]
  if (length(all_files) == 0) {
    stop("Specified folder does not contain any .h5 files.")
  }
  
  
  dir_paths <- list.dirs(folder, recursive = FALSE, full.names = TRUE)
  
  file_paths <- list.files(folder, recursive = FALSE, full.names = TRUE)
  file_paths <- file_paths[grepl(".h5", file_paths)]
  
  result <- NULL
  
  expanding_window <- FALSE
  
  # Read directories
  if (length(dir_paths) > 0) {
    for (i in 1:length(dir_paths)) {
      result[[i]] <- read_models_from_folder(dir_paths[i])
    } 
  }
  
  # Read h5 files
  if (length(file_paths) > 0) {
    
    pos_start <- length(result)
    
    for (i in 1:length(file_paths)) {
      result[[pos_start + i]] <- read_model_from_hdf5(filename = file_paths[i])
    }
    
    expanding_window <- grepl("ExpWind", file_paths[1])
    
  }
  
  if (expanding_window) {
    class(result) <- c("expandingwindow", "list") 
  } else {
    class(result) <- c("modellist", "list") 
  }
  
  return(result)
}