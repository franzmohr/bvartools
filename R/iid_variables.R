# Endogenous variables whose equations carry no coefficients.
#
# `iid` names them, rather than counting them, because a name says which series
# the restriction is about and a count does not -- and because the restriction
# is a statement about a variable, so the argument that carries it should look
# like the one every other function of this package takes.
#
# The sampler, on the other hand, wants a count: it restricts the equations of
# the first `n_iid` variables, and nothing else, because that is what turns
# "these equations" into "these positions of the coefficient vector" without a
# list of them. So the names have to be the leading variables of the data, and
# this is where that is checked and where a data set in the wrong order is told
# what to do about it.
.check_iid_variables <- function(iid, varnames, structural, varsel, tvp) {

  if (is.null(iid)) {
    return(0L)
  }

  if (!is.character(iid) || length(iid) == 0 || anyNA(iid)) {
    stop("Argument 'iid' must name endogenous variables of 'data'.", call. = FALSE)
  }

  unknown <- unique(iid[!iid %in% varnames])
  if (length(unknown) > 0) {
    stop("Argument 'iid' names ", if (length(unknown) > 1) "variables " else "a variable ",
         paste0("'", unknown, "'", collapse = ", "),
         ", which 'data' does not contain.", call. = FALSE)
  }

  if (anyDuplicated(iid) > 0) {
    stop("Argument 'iid' names ",
         paste0("'", unique(iid[duplicated(iid)]), "'", collapse = ", "),
         " more than once.", call. = FALSE)
  }

  if (length(iid) >= length(varnames)) {
    stop("Argument 'iid' names all ", length(varnames), " endogenous variables, which leaves ",
         "no equation with dynamics. At least one variable has to be modelled for the ",
         "restricted ones to be correlated with.", call. = FALSE)
  }

  # The leading variables of the data, in any order among themselves: which of
  # them is first makes no difference, only that none of the modelled ones comes
  # before them.
  leading <- varnames[seq_along(iid)]
  if (!setequal(iid, leading)) {
    stop("Argument 'iid' names ", paste0("'", sort(iid), "'", collapse = ", "),
         ", which are not the first ", length(iid), " columns of 'data' -- those are ",
         paste0("'", leading, "'", collapse = ", "),
         ". The restriction is carried by the position of a variable, so the variables ",
         "without dynamics have to come first. Reorder the columns of 'data'.", call. = FALSE)
  }

  if (isTRUE(structural)) {
    stop("Arguments 'iid' and 'structural' cannot be combined: the contemporaneous ",
         "coefficients of a structural model are laid out by a different rule, so which of ",
         "them belong to a restricted equation is not defined.", call. = FALSE)
  }

  # Only the constant-coefficient samplers read the restriction. Refused here
  # rather than by the sampler, so that a specification that cannot be estimated
  # is refused where it is written.
  if (isTRUE(tvp)) {
    stop("Arguments 'iid' and 'tvp' cannot be combined: an equation restricted to carry no ",
         "coefficients is only available for models whose coefficients are constant.",
         call. = FALSE)
  }

  if (!identical(varsel, "none")) {
    stop("Arguments 'iid' and 'varsel' cannot be combined: the restricted equations' ",
         "coefficients are fixed at zero, so there is nothing for variable selection to ",
         "decide about them.", call. = FALSE)
  }

  return(length(iid))
}
