#' @include persistence_profiles.R
NULL

#' @rdname persistence_profiles.bvecmodel
#' @export
#' @method persistence_profiles bvarmodel
persistence_profiles.bvarmodel <- function(object, ...) {
  stop("A VAR model has no cointegrating relations, so it has no persistence ",
       "profiles: the statistic is the response of a cointegrating relation to ",
       "a system-wide shock, and there is none to follow. Estimate the model as ",
       "a VEC with create_bvecmodel() if the levels are cointegrated, or use ",
       "irf() for the responses of the variables themselves.", call. = FALSE)
}
