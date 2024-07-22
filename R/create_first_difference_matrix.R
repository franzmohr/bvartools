#' Create First Difference Matrix
#' 
#' Creates a sparse first difference matrix.
#' 
#' @param k integer describing the size of an individual block.
#' @param tt integer specifying how often the blocks should be repeated.
#' 
#' @details
#' The function is used in algorithms for time varying parameters.
#' 
#' @return An object of class '\link[Matrix]{dgTMatrix}'.
#' 
#' @export
create_first_difference_matrix <- function(k, tt) {
  return(Matrix::spMatrix(nrow = tt * k,
                          ncol = tt * k,
                          i = c(1:(k * tt), (k + 1:(k * (tt - 1)))),
                          j = c(1:(k * tt), (1:(k * (tt - 1)))),
                          x = c(rep(1, k * tt), rep(-1, k * (tt - 1)))))
}
