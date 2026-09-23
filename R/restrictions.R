# What the two set-identification functions need of a model, and of the table
# of restrictions they are given.
#
# add_sign_restrictions() and add_sign_zero_restrictions() impose the same kind
# of restriction by different algorithms, so they ask the same things of the
# model and take the same table. Checking each in one place is what keeps their
# messages from drifting apart, which for two functions a reader is meant to
# choose between would be worse than the duplication.


# A model whose reduced form can be rotated at all.
.refuse_unrotatable <- function(object, caller) {

  if (is.null(object[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Argument 'object' must include draws of the variance-covariance matrix Sigma.")
  }

  # The rotation acts on the Choleski factor of the reduced form covariance. A
  # structural model has already spent its identification on A_0 and stores the
  # covariance of the structural errors in its place, so there is no reduced
  # form here to rotate. See irf.bvarmodel for the same restriction.
  if (object[["model"]][["structural"]]) {
    stop(caller, " are not defined for a structural model: they would rotate the ",
         "covariance of the structural errors instead of the reduced form. Estimate the model ",
         "with 'structural = FALSE' to identify it by ", tolower(caller), ".")
  }

  # Rotating a covariance that was never estimated is rotating a diagonal
  # matrix, which produces responses that are not those of any shock the data
  # speak about.
  if (object[["model"]][["error"]] %in% c("gamma", "sv", "ald")) {
    stop(caller, " need a model whose error covariances are estimated. Argument ",
         "'object' was estimated with error = \"", object[["model"]][["error"]], "\", which ",
         "leaves the off-diagonal elements at zero, so its Choleski factor is diagonal and a ",
         "rotation of it carries no information about the correlation of the errors.")
  }

  if (object[["model"]][["p"]] == 0) {
    stop(caller, " are only supported for models with p > 0.")
  }
}


# The columns every restriction table has, with the variable names translated
# into the positions the workers count in.
#
# Doing it once here means a worker never sees a name, and an error message
# never shows a number that the caller did not write. What each function makes
# of the 'sign' and 'horizon' columns is its own business and stays with it.
.check_restriction_columns <- function(restrictions, varnames) {

  if (!is.data.frame(restrictions)) {
    stop("Argument 'restrictions' must be a data frame.")
  }

  if (nrow(restrictions) == 0) {
    stop("Argument 'restrictions' does not contain any restriction.")
  }

  required <- c("impulse", "response", "sign")
  absent <- required[!required %in% names(restrictions)]
  if (length(absent) > 0) {
    stop("Argument 'restrictions' must contain the column",
         if (length(absent) > 1) "s " else " ",
         paste0("'", absent, "'", collapse = ", "), ".")
  }

  if (is.null(restrictions[["horizon"]])) {
    restrictions[["horizon"]] <- 0L
  }

  for (i in c("impulse", "response")) {
    name <- as.character(restrictions[[i]])
    unknown <- unique(name[!name %in% varnames])
    if (length(unknown) > 0) {
      stop("Column '", i, "' of argument 'restrictions' names ",
           if (length(unknown) > 1) "variables " else "a variable ",
           paste0("'", unknown, "'", collapse = ", "),
           ", which the model does not contain.")
    }
    restrictions[[i]] <- match(name, varnames)
  }

  return(restrictions)
}
