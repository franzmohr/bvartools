
# Writes the selection scheme of the covariance block where BayesTS reads it.
#
# The scheme of the coefficients is the attribute 'varsel' of /model. That of
# the covariance block of a time varying model is the attribute 'varsel' of
# /model/priors/psi, and the samplers look for it nowhere else. It used to be
# written only as the dataset /priors/psi/varsel, beside the inclusion priors,
# which the samplers never read: they took the scheme to be "none", ignored the
# inclusion priors and the starting values of the indicators, and ran without
# the selection over the covariances the model asked for. The dataset is still
# written with the rest of the prior, because that is what the reader rebuilds
# 'priors' from.
.hdf5_write_psi_varsel <- function(handles, group_model, prior_psi) {

  varsel <- prior_psi[["varsel"]]
  if (is.null(varsel)) {
    return(invisible(NULL))
  }

  group_priors <- .hdf5_group(handles, group_model, "priors")
  group_psi <- .hdf5_group(handles, group_priors, "psi")
  .hdf5_write_attr(group_psi, "varsel", varsel)

  return(invisible(NULL))
}

# Name of a starting value inside /initial of a model file.
#
# R calls the starting error precision of every gamma model u_omega_inv. The
# model file follows BayesTS, whose constant coefficient gamma samplers
# (VarNormalGamma, VecNormalGamma) read /initial/u_sigma_inv and only the time
# varying ones /initial/u_omega_inv. The Rcpp glue maps the name for a model
# estimated inside R, but a file written under the R name left BayesTS without a
# starting precision, and it refused the file. 'direction' is "write" for the
# name in the file and "read" for the name in R.
.hdf5_initial_name <- function(model, name, direction) {
  constant_gamma <- isTRUE(model[["error"]] %in% c("gamma", "gamma+covar")) &&
    !isTRUE(model[["tvp"]])
  if (!constant_gamma) {
    return(name)
  }
  if (direction == "write" && name == "u_omega_inv") {
    return("u_sigma_inv")
  }
  if (direction == "read" && name == "u_sigma_inv") {
    return("u_omega_inv")
  }
  name
}

# A constraint set, the six one-dimensional datasets of
# /data/<part>/constraints. Positions are stored as doubles, which is how
# BayesTS reads them, and one-based as R holds them.
.hdf5_write_constraints <- function(handles, group_data, part, constraints) {

  if (is.null(constraints)) {
    return(invisible(NULL))
  }
  group_part <- .hdf5_group(handles, group_data, part)
  group <- .hdf5_group(handles, group_part, "constraints")
  for (i in c("value", "group", "row", "period", "variable", "weight")) {
    .hdf5_write(group, i, as.numeric(constraints[[i]]))
  }

  invisible(NULL)
}

.hdf5_read_constraints <- function(group) {

  result <- list()
  for (i in c("value", "group", "row", "period", "variable", "weight")) {
    result[[i]] <- as.vector(hdf5r::readDataSet(group[[i]]))
  }

  result
}
