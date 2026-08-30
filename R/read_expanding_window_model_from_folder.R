#' Import Models from HDF5 Files
#' 
#' Imports model information and posterior draws from an HDF5 file.
#' 
#' @param folder Path to a folder with HDF5 files containing model data.
#' 
#' @export
read_expanding_window_model_from_folder <- function(folder) {
  
  
  dir_paths <- list.dirs(folder, recursive = FALSE)
  if (length(dir_paths) > 0) {
    stop("Specified directory contains subfolders.")
  }
  
  file_paths <- list.files(folder, full.names = TRUE)
  
  if (length(file_paths) == 0) {
    stop("No files found.")
  }

  result <- NULL
  for (i in 1:length(file_paths)) {
    result[[i]] <- import_model_from_hdf5(filename = file_paths[i])
  }
  class(result) <- c("expandwindmodellist", "list") 
  
  return(result)
}