
#' @export
get_model_specifications.expandingwindow <- function(object, ...) {
  
  result <- data.frame(type = object[[1]][["model"]][["type"]],
                       k = object[[1]][["model"]][["k"]],
                       p = object[[1]][["model"]][["p"]],
                       m = object[[1]][["model"]][["m"]],
                       s = object[[1]][["model"]][["s"]],
                       n = object[[1]][["model"]][["n"]],
                       varsel = object[[1]][["model"]][["varsel"]])

  return(result)
}