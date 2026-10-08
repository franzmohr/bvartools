#' Read a VAR or VEC Model from an Open File
#'
#' The methods of \code{\link[bayests]{read_hdf5_model}} that
#' \code{\link[bayests]{read_model_from_hdf5}} reaches for a file of a VAR or VEC model:
#' one written by \code{\link[bayests]{write_to_hdf5}}, which records the class in
#' \code{/model/rclass}, or by the BayesTS command line, whose algorithm names
#' it.
#'
#' @param x an empty object of class 'bvarmodel' or 'bvecmodel', the class to
#' read.
#' @param h5_root the group of the open file the model hangs under.
#' @param draws the positions of the draws to read, or \code{NULL} for every
#' draw.
#' @param ... not used.
#'
#' @details The groups are read one by one rather than as a tree, since this
#' layout predates the tree: a time series is restored with its \code{tsp},
#' the scale and centre of a rescaled error correction term are put back, and
#' the sign restrictions and dummy variables, which are tables, are read into
#' data frames.
#'
#' @return An object of class 'bvarmodel' or 'bvecmodel' with the elements the
#' file holds: \code{model}, \code{data}, \code{priors}, \code{initial} and,
#' where the model has been estimated, \code{posterior}. The draws of a sampler
#' are \code{\link[coda]{mcmc}} objects; the per-period posterior of a
#' discounted model is a plain matrix, since its rows are periods rather than
#' draws. \code{\link{bvartools_model}} describes the elements one by one.
#'
#' @export
read_hdf5_model.bvarmodel <- function(x, h5_root, draws = NULL, ...) {
  .read_var_layout(x, h5_root, draws)
}

#' @rdname read_hdf5_model.bvarmodel
#' @export
read_hdf5_model.bvecmodel <- function(x, h5_root, draws = NULL, ...) {
  .read_var_layout(x, h5_root, draws)
}

# The groups of a VAR or VEC model file, read one by one.
.read_var_layout <- function(x, h5_root, draws) {

  h5_names <- names(h5_root)
  
  result <- NULL
  
  if ("model" %in% h5_names) {
    result[["model"]] <- hdf5r::h5attributes(h5_root[["model"]])

    # Written as a group of its own, because the restriction table is a table
    # and the rest of a specification is a set of values. See
    # write_to_hdf5.bvarmodel.
    if ("sign_restrictions" %in% names(h5_root[["model"]])) {
      group_sign <- h5_root[["model"]][["sign_restrictions"]]
      dataset <- group_sign[["restrictions"]]
      restrictions <- .hdf5_read_matrix(dataset)
      colnames(restrictions) <- hdf5r::h5attr(dataset, "columns")

      period <- NULL
      if ("period" %in% hdf5r::h5attr_names(group_sign)) {
        period <- hdf5r::h5attr(group_sign, "period")
      }

      result[["model"]][["sign_restrictions"]] <- list(
        "restrictions" = as.data.frame(restrictions),
        "max_tries" = hdf5r::h5attr(group_sign, "max_tries"),
        "period" = period
      )
    }

    # The rules the dummy variables of add_dummy_variables() are continued by
    # in a forecast, a table like the sign restrictions.
    if ("dummy_variables" %in% names(h5_root[["model"]])) {
      result[["model"]][["dummy_variables"]] <-
        .hdf5_read_dummy_variables(h5_root[["model"]][["dummy_variables"]])
    }
  }
  
  if ("data" %in% h5_names) {
    
    result[["data"]] <- list()
    
    if ("original" %in% names(h5_root[["data"]])) {
      result[["data"]][["original"]] <- list()
      for (i in c("endogen", "exogen", "deterministic")) {
        if (i %in% names(h5_root[["data"]][["original"]])) {
          dataset <- h5_root[["data"]][["original"]][[i]]
          result[["data"]][["original"]][[i]] <- stats::ts(.hdf5_read_matrix(dataset), class = c("mts", "ts", "matrix"))
          dimnames(result[["data"]][["original"]][[i]]) <- list(NULL, hdf5r::h5attr(dataset, "variables"))
          stats::tsp(result[["data"]][["original"]][[i]]) <- hdf5r::h5attr(dataset, "tsp")
          result[["data"]][["original"]][[i]] <- .hdf5_restore_class(result[["data"]][["original"]][[i]], dataset)
        }
      }
    }
    
    if ("train" %in% names(h5_root[["data"]])) {
      result[["data"]][["train"]] <- list()
      for (i in c("y", "w", "x")) {
        if (i %in% names(h5_root[["data"]][["train"]])) {
          dataset <- h5_root[["data"]][["train"]][[i]]
          variables <- hdf5r::h5attr(dataset, "variables")
          result[["data"]][["train"]][[i]] <- stats::ts(.hdf5_read_matrix(dataset))
          dimnames(result[["data"]][["train"]][[i]]) <- list(NULL, variables)
          stats::tsp(result[["data"]][["train"]][[i]]) <- hdf5r::h5attr(dataset, "tsp")
          result[["data"]][["train"]][[i]] <- .hdf5_restore_class(result[["data"]][["train"]][[i]], dataset)

          # A model that was exported while its error correction term was
          # scaled or centred carries the factors it was divided by and the
          # means it lost. They are named after the variables, which is why the
          # names are not stored separately.
          for (name in c("scale", "centre")) {
            if (name %in% hdf5r::h5attr_names(dataset)) {
              factors <- hdf5r::h5attr(dataset, name)
              names(factors) <- variables
              attr(result[["data"]][["train"]][[i]], name) <- factors
            }
          }
        }
      }
      for (i in c("z")) {
        if (i %in% names(h5_root[["data"]][["train"]])) {
          result[["data"]][["train"]][[i]] <- .hdf5_read_matrix(h5_root[["data"]][["train"]][[i]])
        }
      }
    }
    
    if ("forecast" %in% names(h5_root[["data"]])) {
      result[["data"]][["forecast"]] <- list()
      # `x` is the compact layout, `z` the SUR one written before it. Both are
      # read under their own name; the C++ side takes either.
      if ("x" %in% names(h5_root[["data"]][["forecast"]])) {
        result[["data"]][["forecast"]][["x"]] <- .hdf5_read_matrix(h5_root[["data"]][["forecast"]][["x"]])
      }
      if ("z" %in% names(h5_root[["data"]][["forecast"]])) {
        result[["data"]][["forecast"]][["z"]] <- .hdf5_read_matrix(h5_root[["data"]][["forecast"]][["z"]])
      }
    }

    # The values the horizon realised, one row per period and one column per
    # variable. Read under /data like the rest of what a model is given, and
    # unlike the errors taken against them, which are draws and belong to the
    # posterior.
    if ("test" %in% names(h5_root[["data"]])) {
      if ("y" %in% names(h5_root[["data"]][["test"]])) {
        result[["data"]][["test"]] <- list(
          "y" = .hdf5_read_matrix(h5_root[["data"]][["test"]][["y"]]))
        # Written without names, since they are the training sample's: one
        # column per endogenous variable, in the same order.
        train_names <- colnames(result[["data"]][["train"]][["y"]])
        if (length(train_names) == ncol(result[["data"]][["test"]][["y"]])) {
          colnames(result[["data"]][["test"]][["y"]]) <- train_names
        }
      }
    }

    # What was observed of a panel not observed whole; see write_to_hdf5().
    for (i in c("train", "test", "forecast")) {
      if (i %in% names(h5_root[["data"]]) &&
          "constraints" %in% names(h5_root[["data"]][[i]])) {
        result[["data"]][[i]][["constraints"]] <-
          .hdf5_read_constraints(h5_root[["data"]][[i]][["constraints"]])
      }
    }
  }
  
  
  if ("priors" %in% h5_names) {
    
    result[["priors"]] <- list()
    
    for (i in names(h5_root[["priors"]])) {
      result[["priors"]][[i]] <- list()
      for (j in names(h5_root[["priors"]][[i]])) {
        element <- h5_root[["priors"]][[i]][[j]]
        if (inherits(element, "H5Group")) {
          # A prior of its own below the block's: a$shrinkage.
          result[["priors"]][[i]][[j]] <- list()
          for (l in names(element)) {
            result[["priors"]][[i]][[j]][[l]] <- .hdf5_read_value(element[[l]])
          }
        } else {
          result[["priors"]][[i]][[j]] <- .hdf5_read_value(element)
        }
      }
    }
  }
  
  
  if ("initial" %in% h5_names) {
    
    result[["initial"]] <- list()
    
    for (i in names(h5_root[["initial"]])) {
      result[["initial"]][[.hdf5_initial_name(result[["model"]], i, "read")]] <-
        .hdf5_read_value(h5_root[["initial"]][[i]])
    }
  }
  
  
  if ("posterior" %in% h5_names) {
    
    result[["posterior"]] <- list()
    
    for (i in names(h5_root[["posterior"]])) {
      element <- h5_root[["posterior"]][[i]]
      if (inherits(element, "H5Group")) {
        # Whatever draws the group holds, rather than the names they were
        # expected to have. Each is written with its own start, end and
        # thinning interval, so a name left out of a list here is a dataset
        # that is in the file and read back as nothing. `sigma`, which every
        # time varying model keeps beside its coefficients, was exactly that.
        for (j in names(element)) {
          dataset <- element[[j]]
          result[["posterior"]][[i]][[j]] <- .read_posterior_block(dataset, draws)
        }
      } else {
        # The draws the writer keeps on their own rather than in a group,
        # loglik among them. Told apart from a group by what they are rather
        # than by name, so that another one needs nothing here -- which is what
        # let the forecast become a group, read as posterior$forecast$forecasts
        # beside its errors, without a line of this changing.
        result[["posterior"]][[i]] <- .read_posterior_block(element, draws)
      }
    }
  }
  
  
  # Read to decide the class, and not an element of the specification.
  result[["model"]][["rclass"]] <- NULL
  result[["model"]][["rpackage"]] <- NULL
  class(result) <- class(x)
  if (inherits(result, "bvecmodel")) {
    result <- .name_beta_draws(result)
  }
  
  return(result)
}
