#' Spillover Index
#'
#' Produces the connectedness measures of Diebold and Yilmaz (2012) for every
#' model in an object of class 'modellist'.
#'
#' @param object an object of class 'modellist'.
#' @param ... arguments passed forward to \code{\link{spillover.bvarmodel}}.
#'
#' @return A list of objects of class 'bvarspillover'. Elements belonging to
#' models whose estimation failed are \code{NULL}.
#'
#' @export
spillover.modellist <- function(object, ...) {

  result <- lapply(object, function(model_i, ...) {
    if (isTRUE(model_i[["error"]])) {
      return(NULL)
    }
    out <- try(spillover(model_i, ...), silent = TRUE)
    if (inherits(out, "try-error")) {
      return(NULL)
    }
    out
  }, ...)

  return(result)
}
