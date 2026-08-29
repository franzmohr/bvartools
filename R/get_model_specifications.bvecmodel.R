
#' @export
get_model_specifications.bvecmodel <- function(object, ...) {

  result <- data.frame(type = object[["model"]][["type"]],
                       k = object[["model"]][["k"]],
                       p = object[["model"]][["p"]],
                       m = object[["model"]][["m"]],
                       s = object[["model"]][["s"]],
                       n_unrestricted = object[["model"]][["n"]],
                       n_restricted = object[["model"]][["n_restricted"]],
                       rank = object[["model"]][["rank"]])

  # The number of observations is only available, if the object contains data
  tt <- nrow(object[["data"]][["train"]][["y"]])
  if (!is.null(tt)) {
    result[["T"]] <- tt
  }

  result[["varsel"]] <- object[["model"]][["varsel"]]

  return(result)
}
