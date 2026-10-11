# Applies 'f' to every model of a collection -- a 'modellist' or an
# 'expandingwindow' -- and keeps the collection's class and attributes, as
# minnesota_prior.modellist() always has. Used by the prior helpers, which
# return one prior per model in the shape of the collection.
.map_collection <- function(object, f, ...) {
  for (i in seq_along(object)) {
    object[[i]] <- f(object[[i]], ...)
  }
  object
}

# 'f' applied to every model of a collection, the results in a plain list named
# after the models, for the analysis functions whose results are not models. A
# failure stops the call and names the model it happened in, rather than
# leaving a gap that reads as a model without a result.
.analyse_collection <- function(x, f, what, ...) {
  labels <- names(x)
  if (is.null(labels)) {
    labels <- as.character(seq_along(x))
  }
  result <- lapply(seq_along(x), function(i) {
    withCallingHandlers(f(x[[i]], ...), error = function(e) {
      stop(what, "() of model '", labels[i], "': ", conditionMessage(e), call. = FALSE)
    })
  })
  names(result) <- names(x)
  result
}
