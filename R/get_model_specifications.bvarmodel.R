
#' @export
get_model_specifications.bvarmodel <- function(object, ...) {

  result <- data.frame(type = object[["model"]][["type"]],
                       k = object[["model"]][["k"]],
                       p = object[["model"]][["p"]],
                       m = object[["model"]][["m"]],
                       s = object[["model"]][["s"]],
                       n = object[["model"]][["n"]])

  # The number of observations is only available, if the object contains data
  tt <- nrow(object[["data"]][["train"]][["y"]])
  if (!is.null(tt)) {
    result[["T"]] <- tt
  }

  result[["varsel"]] <- object[["model"]][["varsel"]]

  # Left out of a model that has none, so that the specification of every model
  # this package estimated before the restriction existed prints as it did.
  if (!is.null(object[["model"]][["n_iid"]])) {
    result[["n_iid"]] <- object[["model"]][["n_iid"]]
  }

  return(result)
}
