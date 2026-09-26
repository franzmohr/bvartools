#' Expected Size of a Model
#'
#' A generic function that calculates how much memory a model will take once its
#' posterior draws are complete, and so roughly how much disk space it will take
#' when written to a file.
#'
#' @param object an object of a class, for which a method should be called.
#' @param ... arguments passed forward to method.
#'
#' @details The functions that draw or write a model warn before they start when
#' the result would exceed the limit set by \code{options(bvartools.size_warning)},
#' in bytes, which is 1e9, one gigabyte, unless set otherwise, and \code{Inf}
#' turns the warnings off. \code{\link{add_posterior_coefficients}} checks the
#' size \code{expected_model_size()} predicts for the draws it is about to simulate,
#' and \code{\link{write_to_hdf5}} checks the size of the object it is about to
#' write. A package that adds a model class gets both warnings by adding a method
#' to this generic.
#'
#' @return The value returned by the method for the class of \code{object},
#' as described on the pages of the methods.
#'
#' @seealso Methods: \code{\link{expected_model_size.bvarmodel}},
#' \code{\link{expected_model_size.bvecmodel}},
#' \code{\link{expected_model_size.expandingwindow}},
#' \code{\link{expected_model_size.modellist}}.
#'
#' @export
expected_model_size <- function(object, ...) {
  UseMethod("expected_model_size")
}

#' @export
expected_model_size.default <- function(object, ...) {
  stop("expected_model_size() has no method for an object of class '",
       paste(class(object), collapse = "', '"), "'. It knows the models of class ",
       "'bvarmodel' and 'bvecmodel', which create_bvarmodel() and create_bvecmodel() ",
       "return, and the lists of them of class 'modellist' and 'expandingwindow'.")
}

#' @export
print.modelsize <- function(x, ...) {

  steps <- unique(x[["step"]])
  totals <- vapply(steps, function(step) sum(x[["bytes"]][x[["step"]] == step]), numeric(1))
  labels <- ifelse(steps == "", "data, priors and initial values", paste0(steps, "()"))

  if (is.null(x[["model"]])) {
    cat("Expected size of the model once its posterior is complete:",
        .format_bytes(sum(x[["bytes"]])), "\n\n")
  } else {
    cat("Expected size of the", length(unique(x[["model"]])),
        "models once their posteriors are complete:", .format_bytes(sum(x[["bytes"]])), "\n\n")
  }
  print(data.frame(" " = format(labels), size = .format_bytes(totals),
                   check.names = FALSE, row.names = NULL),
        row.names = FALSE, right = FALSE)

  if (is.null(x[["model"]])) {
    blocks <- as.data.frame(x)[x[["step"]] != "",c("element", "draws", "columns", "bytes")]
    blocks[["bytes"]] <- .format_bytes(blocks[["bytes"]])
    names(blocks)[4] <- "size"
    cat("\n")
    print(blocks, row.names = FALSE, right = FALSE)
  }

  invisible(x)
}

# A number of bytes as a person reads it, in powers of 1000 as a drive's
# capacity is given.
.format_bytes <- function(bytes) {
  vapply(bytes, function(b) {
    format(structure(b, class = "object_size"), units = "auto", standard = "SI")
  }, character(1))
}

# The size of every block of draws a VAR or VEC model will hold, from its
# specification alone. The numbers are the allocations of the samplers in
# src/core/models/ and of their bindings in src/*.cpp: every block is a matrix
# of doubles with one row per kept draw, and what decides its width is written
# beside it. test-expected_model_size.R holds each of them against a fitted model, so
# a sampler that stores something new fails there rather than here.
.model_size <- function(object, chains = NULL) {

  model <- object[["model"]]
  priors <- object[["priors"]]
  k <- model[["k"]]
  tt <- NROW(object[["data"]][["train"]][["y"]])
  chains <- .check_chains(object, chains)

  elements <- character(0)
  steps <- character(0)
  draws <- numeric(0)
  columns <- numeric(0)
  add <- function(element, step, n_draws, n_columns) {
    if (n_draws > 0 && n_columns > 0) {
      elements <<- c(elements, element)
      steps <<- c(steps, step)
      draws <<- c(draws, n_draws)
      columns <<- c(columns, n_columns)
    }
  }
  coef_step <- "add_posterior_coefficients"
  vec <- inherits(object, "bvecmodel")
  rank <- if (vec) model[["rank"]] else 0
  n_beta <- if (vec) model[["k_beta"]] * rank else 0

  if (.is_discount(object)) {

    # Not a chain: one row per period of a closed form, and the fixed
    # cointegration matrix of a VEC as a single row.
    n_design <- NCOL(object[["data"]][["train"]][["x"]]) + rank
    add("posterior$a$mean", coef_step, tt, n_design * k)
    add("posterior$a$scale", coef_step, tt, n_design * k)
    add("posterior$a$cov", coef_step, tt, n_design^2)
    add("posterior$u_sigma$scale", coef_step, tt, k^2)
    add("posterior$df", coef_step, tt, 1)
    add("posterior$beta$coeffs", coef_step, 1, n_beta)
    add("posterior$loglik", "add_posterior_loglik", 1, tt)

  } else {

    n <- as.numeric(model[["iterations"]]) * chains
    tvp <- isTRUE(model[["tvp"]])
    periods <- if (tvp) tt else 1
    error <- model[["error"]]
    if (is.null(error)) {
      error <- "wishart"
    }
    sv <- grepl("^sv", error)
    ald <- error == "ald"
    covar <- grepl("covar", error) && k > 1

    # The draws of the non-centred parameterisation, beside the state variances
    # of a block whose prior is 'omega_v'. The prior of the log-volatilities is
    # 'u_sigma' and their draws are kept under 'u_sigma_inv'.
    noncentred <- function(block, width, element = block) {
      if (tvp && !is.null(priors[[block]][["omega_v"]])) {
        add(paste0("posterior$", element, "$omega"), coef_step, n, width)
        add(paste0("posterior$", element, "$omega_log_zero"), coef_step, n, width)
        add(paste0("posterior$", element, "$omega_log_zero_joint"), coef_step, n, 1)
      }
    }

    # The coefficients: one per column of the SUR design, which already holds
    # the loadings of a VEC and the contemporaneous terms of a structural model.
    n_a <- NCOL(object[["data"]][["train"]][["z"]])
    add("posterior$a$coeffs", coef_step, n, n_a * periods)
    if (tvp) {
      add("posterior$a$sigma", coef_step, n, n_a)
    }
    if (!is.null(model[["varsel"]]) && model[["varsel"]] != "none") {
      add("posterior$a$lambda", coef_step, n, n_a)
    }
    if (n_a > 0) {
      noncentred("a", n_a)
    }

    add("posterior$beta$coeffs", coef_step, n, n_beta * periods)
    if (tvp && !is.null(priors[["beta"]][["rho_min"]])) {
      add("posterior$beta$rho", coef_step, n, 1)
    }

    if (covar) {
      n_psi <- k * (k - 1) / 2
      add("posterior$psi$coeffs", coef_step, n, k^2 * periods)
      if (tvp) {
        add("posterior$psi$sigma", coef_step, n, n_psi)
      }
      if (!is.null(priors[["psi"]][["varsel"]]) && priors[["psi"]][["varsel"]] != "none") {
        add("posterior$psi$lambda", coef_step, n, k^2)
      }
      noncentred("psi", n_psi)
    }

    # The error precision moves with time under stochastic volatility and the
    # asymmetric Laplace, and with a covariance block that drifts.
    add("posterior$u_sigma_inv$coeffs", coef_step, n,
        k^2 * (if (sv || ald || (tvp && covar)) tt else 1))
    if (sv) {
      add("posterior$u_sigma_inv$sigma", coef_step, n, k)
      noncentred("u_sigma", k, "u_sigma_inv")
    }
    if (error != "wishart") {
      add("posterior$u_omega_inv$coeffs", coef_step, n, k * (if (sv || ald) tt else 1))
    }
    if (ald) {
      add("posterior$u_scale$coeffs", coef_step, n, k)
    }

    add("posterior$loglik", "add_posterior_loglik", n, tt)
  }

  # Forecast draws are counted once add_forecast_input() has set a horizon.
  if (!is.null(model[["h"]])) {
    n_forecast <- as.numeric(model[["iterations"]]) * chains
    add("posterior$forecast$forecasts", "add_posterior_forecasts", n_forecast, k * model[["h"]])
  }

  result <- data.frame(element = c("data, priors and initial values", elements),
                       step = c("", steps),
                       draws = c(NA, draws),
                       columns = c(NA, columns),
                       bytes = c(as.numeric(utils::object.size(object[names(object) != "posterior"])),
                                 8 * draws * columns),
                       stringsAsFactors = FALSE)
  class(result) <- c("modelsize", "data.frame")
  result
}

# The sizes of the models in a list, one block of rows per model and a column
# naming it. A list inside the list names its models by both.
.model_list_size <- function(object, ...) {
  labels <- names(object)
  if (is.null(labels)) {
    labels <- as.character(seq_along(object))
  }
  parts <- lapply(seq_along(object), function(i) {
    part <- expected_model_size(object[[i]], ...)
    if (is.null(part[["model"]])) {
      part <- cbind(model = labels[i], as.data.frame(part), stringsAsFactors = FALSE)
    } else {
      part[["model"]] <- paste(labels[i], part[["model"]], sep = "/")
    }
    as.data.frame(part)
  })
  result <- do.call(rbind, parts)
  class(result) <- c("modelsize", "data.frame")
  result
}

# The warnings expected_model_size() exists for.
#
# Each is raised once per call from the outside: the collection methods call
# the same generic for every model in them, which would warn once per model
# about the size of one model rather than once about the size of all of them.
# The flag is only set where the size could be computed, so that a list of a
# class this package has no method for -- a 'gvarmodel' of bgvars without its
# own method -- still has its members checked one by one.
#
# A check never stops the call it guards: a size that cannot be computed is not
# a reason to refuse an estimation that may well fit.
.size_check <- new.env(parent = emptyenv())

.size_limit <- function() {
  limit <- getOption("bvartools.size_warning", 1e9)
  if (!is.numeric(limit) || length(limit) != 1 || is.na(limit)) {
    return(Inf)
  }
  limit
}

# Called by add_posterior_coefficients() before it dispatches. Returns whether
# the object was checked, which is when the caller holds the flag.
.check_posterior_size <- function(object, chains = NULL) {
  size <- tryCatch(expected_model_size(object, chains = chains), error = function(e) NULL)
  if (is.null(size)) {
    return(FALSE)
  }
  bytes <- sum(size[["bytes"]][size[["step"]] %in% c("", "add_posterior_coefficients")])
  limit <- .size_limit()
  if (bytes > limit) {
    warning("Once its posterior is drawn, this ",
            if (is.null(size[["model"]])) "model" else "list of models",
            " is expected to take about ", .format_bytes(bytes),
            " -- of memory, or of disk space where it is kept in files -- which is more ",
            "than the limit of ", .format_bytes(limit), " in options(bvartools.size_warning). ",
            "If this machine cannot spare that, lower 'iterations' or raise 'thin' where the ",
            "model was created; expected_model_size() shows which part of the model takes the space. ",
            "options(bvartools.size_warning = Inf) turns this warning off.", call. = FALSE)
  }
  TRUE
}

# Called by write_to_hdf5() before it dispatches. What is written is the object
# as it stands, so its size is measured rather than predicted.
.check_file_size <- function(object) {
  bytes <- as.numeric(utils::object.size(object))
  limit <- .size_limit()
  if (bytes > limit) {
    warning("Writing this object is expected to take up to ", .format_bytes(bytes),
            " of disk space, which is more than the limit of ", .format_bytes(limit),
            " in options(bvartools.size_warning). The files are compressed, but the draws ",
            "of a sampler compress poorly, so make sure the drive has that much room. ",
            "options(bvartools.size_warning = Inf) turns this warning off.", call. = FALSE)
  }
  invisible(NULL)
}
