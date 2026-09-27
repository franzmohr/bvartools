# Internal helpers of add_dummy_variables() and of the functions that have to
# know about the dummies it adds: the forecast input, which continues them over
# the horizon, the expanding windows, which leave out a dummy whose period a
# window has not reached, and the HDF5 files, which keep the rule each dummy was
# built by.
#
# A dummy is an ordinary deterministic term. It is appended to the end of the
# deterministic block of the regressors, so everything that finds that block by
# model$n -- the priors, the forecast, vec_to_var() -- finds the dummies without
# being told about them. What it cannot find by position is how a dummy goes on
# after the sample, which is what model$dummy_variables records: one row per
# dummy with its name, its type ("impulse", "step" or "data") and, for the
# first two, the period it refers to as a value of time().


# The periods of argument 'impulse' or 'step' as values of time(), one per
# dummy. A period is given the way stats::ts() takes 'start', c(2020, 2), or as
# a single number, 2020.25; several are a list of them.
.dummy_periods <- function(periods, argument, frequency) {

  if (is.null(periods)) {
    return(numeric(0))
  }
  if (!is.list(periods)) {
    periods <- list(periods)
  }

  vapply(periods, function(period) {
    if (!is.numeric(period) || !length(period) %in% 1:2 || anyNA(period)) {
      stop("Argument '", argument, "' must be a period such as c(2020, 2), or a list of ",
           "such periods, one per dummy variable.", call. = FALSE)
    }
    if (length(period) == 1) {
      return(as.numeric(period))
    }
    if (period[2] != round(period[2]) || period[2] < 1 || period[2] > frequency) {
      stop("The period c(", period[1], ", ", period[2], ") in argument '", argument,
           "' does not exist at the frequency of the model, which has ", frequency,
           " period(s) per year.", call. = FALSE)
    }
    period[1] + (period[2] - 1) / frequency
  }, numeric(1))
}

# The part of a dummy's name that says which period it is for: 2020Q2 for
# quarterly data, 2020M04 for monthly, 2020 for annual and 2020.2 otherwise.
.dummy_period_name <- function(time, frequency) {
  year <- floor(time + 1e-6)
  period <- round((time - year) * frequency) + 1
  switch(as.character(frequency),
         "1" = as.character(year),
         "4" = paste0(year, "Q", period),
         "12" = sprintf("%dM%02d", year, period),
         paste0(year, ".", period))
}

# The values of the impulse and step dummies in the periods 'times'.
.dummy_values <- function(type, time, times, frequency) {
  tol <- 0.1 / frequency
  if (type == "impulse") {
    return(as.numeric(abs(times - time) < tol))
  }
  as.numeric(times > time - tol)
}

# The dummies of one call as a list of the specification rows and their series
# over the periods 'times': 'data' is the user's time-series object, impulse and
# step dummies are built from their periods.
.dummy_specification <- function(impulse, step, data, frequency) {

  impulse <- .dummy_periods(impulse, "impulse", frequency)
  step <- .dummy_periods(step, "step", frequency)

  # paste0() of a prefix and nothing is the prefix, not nothing.
  label <- function(prefix, times) {
    if (length(times) == 0) character(0) else paste0(prefix, vapply(times, .dummy_period_name, "", frequency))
  }
  spec <- data.frame(name = c(label("impulse.", impulse), label("step.", step)),
                     type = c(rep("impulse", length(impulse)), rep("step", length(step))),
                     time = c(impulse, step),
                     stringsAsFactors = FALSE)

  if (!is.null(data)) {
    if (!"ts" %in% class(data)) {
      stop("Argument 'data' must be a time-series object of class 'ts'.", call. = FALSE)
    }
    if (abs(stats::frequency(data) - frequency) > 1e-8) {
      stop("Argument 'data' has frequency ", stats::frequency(data), ", but the model has ",
           frequency, ".", call. = FALSE)
    }
    data_names <- colnames(data)
    if (is.null(data_names) || any(is.na(data_names) | data_names == "")) {
      stop("The columns of argument 'data' must be named, since the names are what the ",
           "dummy variables are called in the output. Set them with colnames().",
           call. = FALSE)
    }
    spec <- rbind(spec, data.frame(name = data_names, type = "data", time = NA_real_,
                                   stringsAsFactors = FALSE))
  }

  if (nrow(spec) == 0) {
    stop("Specify at least one dummy variable in argument 'impulse', 'step' or 'data'.",
         call. = FALSE)
  }

  # The forecast input continues a term by its name, taking "const", "trend"
  # and anything with "season" in it for its own, so a dummy called like that
  # would be continued as something it is not.
  reserved <- spec[["name"]] %in% c("const", "trend") | grepl("season", spec[["name"]], fixed = TRUE)
  if (any(reserved)) {
    stop("The name(s) ", paste0("'", spec[["name"]][reserved], "'", collapse = ", "),
         " are taken by the deterministic terms of create_bvarmodel() and create_bvecmodel(). ",
         "Rename the column(s) of argument 'data'.", call. = FALSE)
  }
  if (anyDuplicated(spec[["name"]])) {
    stop("The dummy variable '", spec[["name"]][anyDuplicated(spec[["name"]])],
         "' is specified twice.", call. = FALSE)
  }

  list(spec = spec, data = data)
}

# The series of each dummy over the periods 'times', one column per row of
# 'spec'. A 'data' column that does not cover them is NA there.
.dummy_series <- function(spec, data, times, frequency) {
  result <- matrix(NA_real_, length(times), nrow(spec), dimnames = list(NULL, spec[["name"]]))
  for (i in seq_len(nrow(spec))) {
    if (spec[["type"]][i] == "data") {
      result[, i] <- .ts_values_at(data[, spec[["name"]][i]], times, frequency)
    } else {
      result[, i] <- .dummy_values(spec[["type"]][i], spec[["time"]][i], times, frequency)
    }
  }
  result
}

# The values of the series 'x' in the periods 'times', NA where it has none.
.ts_values_at <- function(x, times, frequency) {
  pos <- match(round(times * frequency), round(stats::time(x) * frequency))
  as.numeric(x)[pos]
}

# Adds the dummies of 'dummies', the output of .dummy_specification(), to one
# model. 'empty' says what happens to a dummy that is zero throughout the
# estimation sample: "stop" refuses it, "drop" leaves it out, which is what a
# window that ends before its period does.
.add_dummy_variables <- function(object, dummies, empty = "stop") {

  added <- c("priors" = "add_priors()", "initial" = "add_initial_values()",
             "posterior" = "add_posterior_coefficients()")
  present <- names(added)[names(added) %in% names(object)]
  if (length(present) > 0 || !is.null(object[["data"]][["forecast"]])) {
    stop("Dummy variables change the regressors of the model, and with them the number ",
         "of coefficients, so add_dummy_variables() must come before add_priors(), ",
         "add_initial_values(), add_posterior_coefficients() and add_forecast_input(). ",
         "Call it on the output of create_bvarmodel() or create_bvecmodel().", call. = FALSE)
  }

  spec <- dummies[["spec"]]
  k <- object[["model"]][["k"]]
  y <- object[["data"]][["train"]][["y"]]
  frequency <- stats::frequency(y)
  times <- as.numeric(stats::time(y))
  x_old <- object[["data"]][["train"]][["x"]]

  values <- .dummy_series(spec, dummies[["data"]], times, frequency)

  missing <- colnames(values)[colSums(is.na(values)) > 0]
  if (length(missing) > 0) {
    stop("The dummy variable(s) ", paste0("'", missing, "'", collapse = ", "), " in argument ",
         "'data' must have a value in every period of the estimation sample, which runs from ",
         .ts_period_label(times[1], frequency), " to ",
         .ts_period_label(times[length(times)], frequency), ".", call. = FALSE)
  }

  zero <- colSums(values != 0) == 0
  if (any(zero)) {
    if (empty == "stop") {
      stop("The dummy variable(s) ", paste0("'", spec[["name"]][zero], "'", collapse = ", "),
           " are zero in every period of the estimation sample, which runs from ",
           .ts_period_label(times[1], frequency), " to ",
           .ts_period_label(times[length(times)], frequency),
           " once the lags are taken, so the data say nothing about their coefficients.",
           call. = FALSE)
    }
    spec <- spec[!zero, , drop = FALSE]
    values <- values[, !zero, drop = FALSE]
  }
  if (nrow(spec) == 0) {
    return(object)
  }

  taken <- intersect(spec[["name"]], c(colnames(x_old), object[["model"]][["endogen"]]))
  if (length(taken) > 0) {
    stop("The model already has a regressor or variable called ",
         paste0("'", taken, "'", collapse = ", "), ".", call. = FALSE)
  }

  # A dummy the other regressors already span leaves its coefficient to the
  # prior and makes the rest arbitrary: a step dummy that is one in every period
  # is the constant again. Checked one at a time, so that the message can name
  # the first one that does it.
  current <- if (is.null(x_old)) matrix(0, length(times), 0) else as.matrix(x_old)
  for (i in seq_len(ncol(values))) {
    candidate <- cbind(current, values[, i])
    if (qr(candidate)$rank < ncol(candidate)) {
      hint <- if (spec[["type"]][i] == "step") {
        " A step dummy that is one in every period of the sample is the constant again."
      } else {
        ""
      }
      stop("The dummy variable '", spec[["name"]][i], "' is a linear combination of the ",
           "regressors the model already has, so its coefficient cannot be told apart from ",
           "theirs.", hint, call. = FALSE)
    }
    current <- candidate
  }

  # Regressors ----
  x_new <- cbind(if (is.null(x_old)) NULL else as.matrix(x_old), values)
  colnames(x_new) <- c(colnames(x_old), spec[["name"]])
  x_new <- stats::ts(x_new, class = c("mts", "ts", "matrix"))
  stats::tsp(x_new) <- stats::tsp(y)
  object <- .replace_regressors(object, x_new)

  # Specification ----
  object[["model"]][["n"]] <- as.integer(object[["model"]][["n"]] + nrow(spec))
  object[["model"]][["deterministic"]] <- c(object[["model"]][["deterministic"]], spec[["name"]])
  rownames(spec) <- NULL
  object[["model"]][["dummy_variables"]] <- rbind(object[["model"]][["dummy_variables"]], spec)

  # The series as given, over the whole span of the original data or, for a
  # 'data' dummy, the span it was given over, which may reach into the periods a
  # forecast is made for.
  original <- object[["data"]][["original"]]
  base <- if (!is.null(original[["deterministic"]])) original[["deterministic"]] else original[["endogen"]]
  base_times <- as.numeric(stats::time(base))
  full <- list()
  for (i in seq_len(nrow(spec))) {
    if (spec[["type"]][i] == "data") {
      full[[i]] <- dummies[["data"]][, spec[["name"]][i]]
    } else {
      full[[i]] <- stats::ts(.dummy_values(spec[["type"]][i], spec[["time"]][i], base_times, frequency),
                             start = stats::start(base), frequency = frequency)
    }
  }
  object[["data"]][["original"]][["deterministic"]] <-
    .merge_series(original[["deterministic"]], full, spec[["name"]])

  object
}

# 'x' with the series in 'series' appended as columns named 'names', over the
# union of their spans.
.merge_series <- function(x, series, names) {
  all <- c(if (is.null(x)) list() else list(x), series)
  result <- do.call(stats::ts.union, all)
  result <- stats::ts(as.matrix(result), class = c("mts", "ts", "matrix"))
  stats::tsp(result) <- stats::tsp(do.call(stats::ts.union, all))
  colnames(result) <- c(colnames(x), names)
  result
}

# Replaces the compact regressors of a model and rebuilds their SUR form around
# what else it holds: the columns of the cointegration term before them in a VEC
# model, and the contemporaneous terms of a structural model after them. A
# discounted model has no SUR form to rebuild.
.replace_regressors <- function(object, x_new) {

  k <- object[["model"]][["k"]]
  x_old <- object[["data"]][["train"]][["x"]]
  z <- object[["data"]][["train"]][["z"]]
  object[["data"]][["train"]][["x"]] <- x_new

  if (is.null(z)) {
    if (!.is_discount(object) && !is.null(x_new) && ncol(x_new) > 0) {
      z <- kronecker(as.matrix(x_new), diag(1, k))
      dimnames(z) <- NULL
      object[["data"]][["train"]][["z"]] <- z
    }
    return(object)
  }

  lead <- 0
  if (inherits(object, "bvecmodel") && isTRUE(object[["model"]][["rank"]] > 0)) {
    lead <- object[["model"]][["rank"]] * k
  }
  width_old <- if (is.null(x_old)) 0 else k * ncol(x_old)
  before <- z[, seq_len(lead), drop = FALSE]
  after <- z[, -seq_len(lead + width_old), drop = FALSE]
  if (lead + width_old == 0) {
    after <- z
  }

  middle <- if (is.null(x_new) || ncol(x_new) == 0) NULL else kronecker(as.matrix(x_new), diag(1, k))
  z <- cbind(before, middle, after)
  dimnames(z) <- NULL
  object[["data"]][["train"]][["z"]] <- z

  object
}

# Leaves out the dummies of a model that are zero in every period of its
# estimation sample, which is what a window that ends before the period of an
# impulse or step dummy has: a forecast made at its end could not have known
# about the event. The coefficients are what gets dropped, so priors already
# added cannot be kept.
.drop_empty_dummy_variables <- function(object) {

  spec <- object[["model"]][["dummy_variables"]]
  x <- object[["data"]][["train"]][["x"]]
  if (is.null(spec) || is.null(x)) {
    return(object)
  }

  empty <- spec[["name"]][colSums(as.matrix(x[, spec[["name"]], drop = FALSE]) != 0) == 0]
  if (length(empty) == 0) {
    return(object)
  }

  if (!is.null(object[["priors"]])) {
    stop("The dummy variable(s) ", paste0("'", empty, "'", collapse = ", "), " are zero in ",
         "every period of the first windows, which end before the period they refer to, and ",
         "are left out of those windows. Their priors were already added for the whole ",
         "sample, so call use_expanding_window() before add_priors().", call. = FALSE)
  }

  keep <- !colnames(x) %in% empty
  x_new <- x[, keep, drop = FALSE]
  if (ncol(x_new) == 0) {
    x_new <- NULL
  }
  object <- .replace_regressors(object, x_new)
  object[["model"]][["n"]] <- as.integer(object[["model"]][["n"]] - length(empty))
  object[["model"]][["deterministic"]] <- setdiff(object[["model"]][["deterministic"]], empty)
  if (length(object[["model"]][["deterministic"]]) == 0) {
    object[["model"]][["deterministic"]] <- NULL
  }
  spec <- spec[!spec[["name"]] %in% empty, , drop = FALSE]
  rownames(spec) <- NULL
  object[["model"]][["dummy_variables"]] <- if (nrow(spec) > 0) spec else NULL

  original <- object[["data"]][["original"]][["deterministic"]]
  if (!is.null(original)) {
    original <- original[, !colnames(original) %in% empty, drop = FALSE]
    object[["data"]][["original"]][["deterministic"]] <- if (ncol(original) > 0) original else NULL
  }

  object
}

# Warns about dummies that window() left without a single non-zero value.
.warn_empty_dummy_variables <- function(object) {
  spec <- object[["model"]][["dummy_variables"]]
  x <- object[["data"]][["train"]][["x"]]
  if (is.null(spec) || is.null(x)) {
    return(invisible(NULL))
  }
  empty <- spec[["name"]][colSums(as.matrix(x[, spec[["name"]], drop = FALSE]) != 0) == 0]
  if (length(empty) > 0) {
    warning("The dummy variable(s) ", paste0("'", empty, "'", collapse = ", "), " are zero in ",
            "every period of the window, so the data say nothing about their coefficients.",
            call. = FALSE)
  }
  invisible(NULL)
}

# The values of the dummies of a model in the periods 'times' of the forecast.
# An impulse dummy is one in its period only and a step dummy from its period
# on; a 'data' dummy takes the values it was given, and a forecast period it
# was not given a value for is refused with the periods it has to cover.
.continue_dummy_variables <- function(object, times) {

  spec <- object[["model"]][["dummy_variables"]]
  frequency <- stats::frequency(object[["data"]][["train"]][["y"]])
  original <- object[["data"]][["original"]][["deterministic"]]

  values <- .dummy_series(spec, original, times, frequency)

  missing <- colnames(values)[colSums(is.na(values)) > 0]
  if (length(missing) > 0) {
    stop("The dummy variable(s) ", paste0("'", missing, "'", collapse = ", "), " were given ",
         "as argument 'data' of add_dummy_variables() without values for the forecast periods ",
         .ts_period_label(times[2], frequency), " to ",
         .ts_period_label(times[length(times)], frequency), ". Give them there, or pass all ",
         "deterministic terms of the forecast in argument 'deterministic'.", call. = FALSE)
  }

  values
}

# Writes the specification of the dummies to group /model/dummy_variables, one
# attribute per column, since the table does not fit in an attribute of /model.
.hdf5_write_dummy_variables <- function(handles, group_model, spec) {
  if (is.null(spec)) {
    return(invisible(NULL))
  }
  group <- .hdf5_group(handles, group_model, "dummy_variables")
  .hdf5_write_attr(group, "name", spec[["name"]])
  .hdf5_write_attr(group, "type", spec[["type"]])
  .hdf5_write_attr(group, "time", as.numeric(spec[["time"]]))
  invisible(NULL)
}

.hdf5_read_dummy_variables <- function(group) {
  time <- as.numeric(hdf5r::h5attr(group, "time"))
  time[is.nan(time)] <- NA_real_
  data.frame(name = as.character(hdf5r::h5attr(group, "name")),
             type = as.character(hdf5r::h5attr(group, "type")),
             time = time,
             stringsAsFactors = FALSE)
}
