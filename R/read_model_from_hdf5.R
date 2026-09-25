#' Import Models from HDF5 Files
#' 
#' Imports model information and posterior draws from an HDF5 file.
#' 
#' @param filename Path to an HDF5 file containing model data.
#' @param group the group the model's tree hangs under inside its file.
#' Defaults to \code{""}, the root of the file, which is where a file holding a
#' single model puts it. See 'Details'.
#' @param draws the draws to read, as their positions in the chain. Defaults to
#' \code{NULL}, every draw. An integer vector reads those draws of every block
#' of the posterior, and \code{integer(0)} reads none of them, which gives the
#' model, its data and its priors without the draws. See 'Details'.
#' 
#' @details
#' 
#' With a \code{group} every path the reader looks for is read under it
#' instead of at the root, so one file can hold several models side by side.
#' \code{\link{list_models_in_hdf5}} reports which groups of a file hold one.
#' 
#' The spelling is the one the BayesTS command line uses for its \code{--group}
#' flag: a leading slash and no trailing slash, with \code{""} for the root.
#'
#' \code{draws} reads part of a chain. Every block of a posterior holds one row
#' per draw, and only the rows asked for are read from the file, so a caller
#' that works through a long chain in pieces -- solving a global model draw by
#' draw, for instance -- never holds more of it than the piece it is working
#' on. \code{integer(0)} reads a model without its draws, which is what a step
#' needs that only looks at what a model is.
#'
#' A partial read cannot describe the chain it came from, so the blocks it
#' returns are labelled as a chain of their own, from one to the number of draws
#' read. A full read keeps the labels the file carries.
#'
#' @return An object of class 'bvarmodel' or 'bvecmodel', depending on the
#' algorithm recorded in the file, with the elements that file holds:
#' \code{model}, \code{data}, \code{priors}, \code{initial} and, where the
#' model has been estimated, \code{posterior}. The draws of a sampler are
#' \code{\link[coda]{mcmc}} objects; the per-period posterior of a discounted
#' model is a plain matrix, since its rows are periods rather than draws.
#' \code{\link{bvartools_model}} describes the elements one by one.
#'
#' @export
read_model_from_hdf5 <- function(filename, group = "", draws = NULL) {
  
  group <- .normalize_hdf5_group(group)
  draws <- .check_read_draws(draws)
  
  h5_file <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (h5_file$is_valid) h5_file$close_all(), add = TRUE)
  
  # Every path below is named against this rather than against the file, which
  # is all a group amounts to on the way in.
  h5_root <- .hdf5_model_root(h5_file, group)
  
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
  } else {
    stop("File ", filename, " does not contain model specification",
         if (group != "") paste0(" in group ", group), ".")
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
  }
  
  
  if ("priors" %in% h5_names) {
    
    result[["priors"]] <- list()
    
    for (i in names(h5_root[["priors"]])) {
      result[["priors"]][[i]] <- list()
      for (j in names(h5_root[["priors"]][[i]])) {
        result[["priors"]][[i]][[j]] <- .hdf5_read_value(h5_root[["priors"]][[i]][[j]])
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
  
  
  h5_file$close_all()
  
  # Determine the class of the returned object based on the used algorithm
  result_class <- result[["model"]][["rclass"]]
  if (is.null(result_class)) {
    bvarmodel <- c("VarNormalAld", "VarNormalGamma", "VarNormalStochvol",
                   "VarNormalWishart", "VarTvpAld", "VarTvpGamma",
                   "VarTvpStochvol", "VarTvpWishart")
    if (result[["model"]][["algorithm"]] %in% bvarmodel) {
      result_class <- c("bvarmodel", "list")
    }
    
    bvarmodel <- c(bvarmodel, "VarTvpDiscount")

    bvecmodel <- c("VecNormalGamma", "VecNormalStochvol", "VecNormalWishart",
                   "VecTvpGamma", "VecTvpStochvol", "VecTvpWishart",
                   "VecKlgs2010", "VecTvpDiscount")
    if (result[["model"]][["algorithm"]] %in% bvecmodel) {
      result_class <- c("bvecmodel", "list")
    } 
  }
  
  # Read to decide the class, and not an element of the specification.
  result[["model"]][["rclass"]] <- NULL
  class(result) <- result_class
  if ("bvecmodel" %in% result_class) {
    result <- .name_beta_draws(result)
  }
  
  return(result)
}


# The draws to read, checked once so that a bad request fails before the file is
# opened rather than on the first block that is read.
.check_read_draws <- function(draws) {

  if (is.null(draws)) {
    return(NULL)
  }
  if (!is.numeric(draws) || anyNA(draws) || any(draws %% 1 != 0)) {
    stop("Argument 'draws' must be whole numbers, the positions of the draws to read.")
  }
  draws <- as.integer(draws)
  if (length(draws) > 0 && min(draws) < 1) {
    stop("Argument 'draws' must be positions of draws, so at least one.")
  }

  draws
}

# One block of a posterior, as a chain where it is one and as a matrix where it
# is not.
#
# Everything a sampler writes carries coda's mcpar -- the start, the end and the
# thinning interval of the chain it came from -- and is read back as an 'mcmc'
# object. The discounted models write no such attributes, and correctly: their
# posterior is one column per period of a closed form rather than one row per
# draw of a chain, so there is no start, no end and nothing thinned. A block
# without them is read as the plain matrix it is. Labelling it as a chain would
# make a row of it look like a draw, and the rows are periods.
.read_posterior_block <- function(dataset, draws) {

  block <- .read_draw_rows(dataset, draws)

  if (!"start" %in% hdf5r::h5attr_names(dataset)) {
    return(block)
  }
  if (!is.null(draws)) {
    return(coda::mcmc(block))
  }

  coda::mcmc(block,
             start = hdf5r::h5attr(dataset, "start"),
             end = hdf5r::h5attr(dataset, "end"),
             thin = hdf5r::h5attr(dataset, "thin"))
}

# One block of a posterior, as a matrix of the draws asked for. Rows are read
# from the file rather than read and then subset, which is what keeps a partial
# read small.
.read_draw_rows <- function(dataset, draws) {

  dims <- dataset$dims
  n_draws <- dims[1]
  n_columns <- if (length(dims) > 1) dims[2] else 1L

  if (is.null(draws)) {
    # Shaped from what the file says rather than from what hdf5r hands back.
    # A dataset with one row comes back as a plain vector, dimensions dropped,
    # and as.matrix() then makes a column of it -- so a block of one draw over
    # many columns was read transposed, as many draws of one column. Nothing a
    # sampler writes has one draw, which is why this went unseen until the
    # discounted models, whose /posterior/loglik is exactly that: one row, the
    # exact pointwise log marginal likelihood over the sample.
    return(matrix(hdf5r::readDataSet(dataset), nrow = n_draws, ncol = n_columns))
  }

  if (length(draws) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = n_columns))
  }
  if (max(draws) > n_draws) {
    stop("Draw ", max(draws), " was asked for, and the chain holds ", n_draws, ".")
  }

  if (length(dims) > 1) {
    matrix(dataset[draws, ], nrow = length(draws), ncol = n_columns)
  } else {
    matrix(dataset[draws], nrow = length(draws), ncol = 1L)
  }
}
