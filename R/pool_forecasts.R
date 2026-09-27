#' Pool Forecasts
#'
#' Combines the forecasts of several models into one equal-weight pool, which
#' is compared with the models like any of them.
#'
#' @param ... two or more models with forecasts: objects of class
#' \code{"expandingwindow"}, or single models of class \code{"bvarmodel"} or
#' \code{"bvecmodel"}, or \code{"modellist"} objects of either. Arguments may be
#' named; the names label the members of the pool.
#'
#' @details A pool of forecasts is a mixture of predictive distributions: the
#' predictive distribution of the pool is the average of those of its members.
#' Combining forecasts in this way is one of the most reliable ways of improving
#' them, because the errors of different models are not perfectly correlated
#' and no single specification is best in every period (Bates and Granger,
#' 1969; Hall and Mitchell, 2007; Geweke and Amisano, 2011).
#'
#' The pool is built from draws the members already hold. For every forecast,
#' the same number of draws is taken at random from each member -- as many as
#' the member with the fewest draws has -- and stacked, so that every member
#' has the same weight. The draws are selected with R's random number
#' generator, so \code{set.seed()} makes a pool reproducible.
#'
#' \strong{The members must forecast the same thing.} The pool contains the
#' endogenous variables all members share, and the forecast horizons all of them
#' reach. A variable forecast by only some of the members is left out, and a
#' pool of models that share no variable is refused. A VEC model forecasts the
#' levels of its variables, so it can be pooled with VAR models of the same
#' levels, but not with VAR models of their growth rates. Forecasts aggregated
#' to annual figures can only be pooled with forecasts aggregated in the same
#' way.
#'
#' The windows of expanding window exercises are matched by the end of their
#' estimation samples, so that each pooled forecast combines forecasts made with
#' the same information. Windows that not every member has are left out, with a
#' message. Single models and expanding windows cannot be pooled with each
#' other, and external forecasts, which are point forecasts rather than
#' predictive distributions, cannot be pooled at all: they are compared beside a
#' pool, by combining both with \code{\link{combine_models}}.
#'
#' The log predictive densities that \code{\link{add_predictive_loglik}} and
#' the forecasts of BayesTS provide are pooled as the forecasts are, which gives
#' the log of the average predictive density of the members -- the score of the
#' mixture, not of its draws. They are kept only if every member carries them
#' and all members have the same endogenous variables, because a predictive
#' density is a joint density of all the variables of a model.
#'
#' The weights are equal. Weights that favour members which forecast better in
#' the past must only use the periods before each forecast, or they flatter the
#' pool with information it did not have.
#'
#' \strong{A pool is formed after the members' forecasts are drawn}, and the
#' functions that estimate a model leave it unchanged. It can be passed to
#' \code{\link{add_forecast_errors}} and \code{\link{selection_criteria}}, and
#' included in a list of models with \code{\link{combine_models}}. It cannot be
#' written to a file: write its members and form the pool again after reading
#' them back.
#'
#' @return An object of class \code{"forecastpool"}, which behaves like an
#' expanding window exercise, or, if single models were pooled, a single pooled
#' forecast of class \code{"poolwindow"}. Each pooled forecast holds the draws
#' of the pool in \code{posterior$forecast$forecasts} and, in \code{model},
#' the names of the members in \code{members} and the number of draws taken from
#' each in \code{draws}.
#'
#' @references
#' Bates, J. M., & Granger, C. W. J. (1969). The combination of forecasts.
#' \emph{Journal of the Operational Research Society, 20}(4), 451--468.
#' \doi{10.1057/jors.1969.103}
#'
#' Geweke, J., & Amisano, G. (2011). Optimal prediction pools. \emph{Journal of
#' Econometrics, 164}(1), 130--141. \doi{10.1016/j.jeconom.2011.02.017}
#'
#' Hall, S. G., & Mitchell, J. (2007). Combining density forecasts.
#' \emph{International Journal of Forecasting, 23}(1), 1--13.
#' \doi{10.1016/j.ijforecast.2006.08.001}
#'
#' @examples
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' train <- window(e1, end = c(1978, 4))
#'
#' # Two specifications, each estimated over expanding windows
#' model <- create_bvarmodel(train, p = 1:2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#' model <- add_priors(model,
#'                     coef = list(v_i = 1 / 10, v_i_det = 1 / 100),
#'                     sigma = list(df = "k", scale = 1))
#' model <- use_expanding_window(model, start = c(1978, 1))
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model)
#' model <- add_forecast_input(model, n_ahead = 2)
#' model <- add_posterior_forecasts(model)
#'
#' # The equal-weight pool of both
#' set.seed(123)
#' pool <- pool_forecasts(p1 = model[[1]], p2 = model[[2]])
#'
#' # Scored and compared like the models
#' all_models <- combine_models(model, pool)
#' all_models <- add_forecast_errors(all_models, test_sample = e1)
#' selection_criteria(all_models)
#'
#' @family model comparison
#' @export
pool_forecasts <- function(...) {

  members <- .pool_members(list(...))
  n_members <- length(members)

  kinds <- vapply(members, function(x) {
    if (inherits(x, "expandingwindow")) "windows" else "single"
  }, character(1))
  if (length(unique(kinds)) > 1) {
    stop("A pool combines either expanding window exercises or single models, not both. ",
         "Pool the windows with the windows and the single models with the single models.")
  }
  windows <- kinds[1] == "windows"

  # Every member as a list of forecasts, named by the end of their estimation
  # sample, which is what matches the forecasts of different members
  forecasts <- lapply(seq_len(n_members), function(i) {
    x <- if (windows) unclass(members[[i]]) else list(members[[i]])
    ends <- vapply(x, .pool_sample_end, numeric(1))
    stats::setNames(x, format(ends, digits = 10))
  })
  names(forecasts) <- names(members)

  freq <- unique(unlist(lapply(forecasts, function(x) {
    lapply(x, function(y) stats::frequency(y[["data"]][["train"]][["y"]]))
  })))
  if (length(freq) > 1) {
    stop("The members of the pool have data of different frequencies. ",
         "Only forecasts of the same periods can be pooled.")
  }

  endogen <- lapply(forecasts, function(x) x[[1]][["model"]][["endogen"]])
  common <- Reduce(intersect, endogen)
  if (length(common) == 0) {
    stop("The members of the pool share no endogenous variable: ",
         paste0(names(members), " (", vapply(endogen, paste, character(1), collapse = ", "), ")",
                collapse = "; "),
         ". A pool mixes forecasts of the same variables. A VEC model forecasts levels, ",
         "which cannot be pooled with forecasts of growth rates.")
  }
  same_variables <- all(vapply(endogen, function(x) setequal(x, common), logical(1)))

  aggregation <- unique(lapply(forecasts, function(x) x[[1]][["model"]][["aggregation"]][["frequency"]]))
  if (length(aggregation) > 1) {
    stop("Only some members of the pool forecast annual figures aggregated from ",
         "higher frequency data. Aggregate all of them or none.")
  }

  ends <- Reduce(intersect, lapply(forecasts, names))
  if (length(ends) == 0) {
    stop("The members of the pool have no forecast made at the end of the same estimation sample. ",
         "Pooled forecasts must use the same information.")
  }
  n_dropped <- sum(vapply(forecasts, function(x) length(setdiff(names(x), ends)), integer(1)))
  if (n_dropped > 0) {
    message(n_dropped, " forecast", ifelse(n_dropped > 1, "s were", " was"),
            " left out, because not every member of the pool has a forecast from the same ",
            "estimation sample.")
  }
  ends <- ends[order(as.numeric(ends))]

  result <- lapply(ends, function(end) {
    .pool_one(lapply(forecasts, function(x) x[[end]]), names(members), common, same_variables)
  })

  if (!windows) {
    return(result[[1]])
  }

  class(result) <- c("forecastpool", "expandingwindow", "list")
  return(result)
}


# The members of a pool from the arguments of pool_forecasts(): list elements
# of a model list become members of their own, and every member gets a name.
.pool_members <- function(input) {

  arg_names <- names(input)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(input))
  }

  members <- list()
  member_names <- character(0)
  for (i in seq_along(input)) {
    x <- input[[i]]
    if (inherits(x, "externalforecast") || inherits(x, "externalwindow")) {
      stop("External forecasts are point forecasts and cannot be pooled with predictive ",
           "distributions. Compare them beside the pool instead, by combining both with ",
           "'combine_models'.")
    }
    if (inherits(x, "forecastpool") || inherits(x, "poolwindow")) {
      stop("A pool cannot be a member of another pool, which would give its members ",
           "unequal weights. Pass the members of both pools to 'pool_forecasts' instead.")
    }
    if (inherits(x, "modellist")) {
      inner <- names(x)
      if (is.null(inner)) {
        inner <- rep("", length(x))
      }
      for (j in seq_along(x)) {
        members <- c(members, list(x[[j]]))
        member_names <- c(member_names, if (nzchar(inner[j])) inner[j] else "")
      }
    } else if (inherits(x, "expandingwindow") || inherits(x, "bvarmodel") ||
               inherits(x, "bvecmodel")) {
      members <- c(members, list(x))
      member_names <- c(member_names, arg_names[i])
    } else {
      stop("Argument ", i, " of 'pool_forecasts' is not a model. The members of a pool are ",
           "objects of class 'expandingwindow', 'bvarmodel' or 'bvecmodel', or a 'modellist' of them.")
    }
  }

  if (length(members) < 2) {
    stop("A pool needs at least two members.")
  }
  empty <- !nzchar(member_names)
  member_names[empty] <- paste0("Model ", which(empty))
  names(members) <- member_names

  return(members)
}


# The end of the estimation sample of a model, which identifies the information
# its forecast was made with.
.pool_sample_end <- function(object) {
  y <- object[["data"]][["train"]][["y"]]
  if (is.null(stats::tsp(y))) {
    stop("A member of the pool has training data without a time index, so its forecasts ",
         "cannot be matched with those of the other members.")
  }
  stats::tsp(y)[2]
}


# One pooled forecast from the forecasts of the members made at the end of the
# same estimation sample.
.pool_one <- function(parts, member_names, common, same_variables) {

  draws <- lapply(seq_along(parts), function(i) {
    x <- .forecast_draws(parts[[i]])
    if (is.null(x)) {
      stop("Member '", member_names[i], "' of the pool has no forecasts for the sample ending in ",
           format(.pool_sample_end(parts[[i]])), ". Use 'add_posterior_forecasts' on it first.")
    }
    as.matrix(x)
  })

  k_m <- vapply(parts, function(x) x[["model"]][["k"]], numeric(1))
  h_m <- vapply(seq_along(parts), function(i) ncol(draws[[i]]) / k_m[i], numeric(1))
  h <- min(h_m)
  k <- length(common)

  # The same number of draws from every member gives each the same weight
  n_draws <- min(vapply(draws, nrow, integer(1)))
  picks <- lapply(draws, function(x) sort(sample.int(nrow(x), n_draws)))

  pooled <- do.call(rbind, lapply(seq_along(parts), function(i) {
    endogen_i <- parts[[i]][["model"]][["endogen"]]
    cols <- as.vector(outer(match(common, endogen_i), (seq_len(h) - 1) * k_m[i], "+"))
    draws[[i]][picks[[i]], cols, drop = FALSE]
  }))
  dimnames(pooled) <- list(NULL, paste0(rep(common, h), "_", rep(seq_len(h), each = k)))

  first <- parts[[1]]
  original <- first[["data"]][["original"]][["endogen"]]
  end <- .pool_sample_end(first)
  original <- original[, common, drop = FALSE]

  model <- list(type = "Pool",
                algorithm = "pool",
                k = k,
                p = 0L,
                m = 0L,
                s = 0L,
                n = 0L,
                varsel = "none",
                endogen = common,
                structural = FALSE,
                tvp = FALSE,
                h = h,
                members = member_names,
                draws = n_draws)
  if (!is.null(first[["model"]][["aggregation"]])) {
    model[["aggregation"]] <- first[["model"]][["aggregation"]]
  }

  data <- list(original = list(endogen = original),
               train = list(y = stats::window(original, end = end)))
  test <- first[["data"]][["test"]][["y"]]
  if (!is.null(test)) {
    test <- as.matrix(test)
    if (all(common %in% colnames(test))) {
      data[["test"]] <- list(y = test[seq_len(min(h, nrow(test))), common, drop = FALSE])
    }
  }

  result <- list(model = model, data = data,
                 posterior = list(forecast = list(forecasts = coda::mcmc(pooled))))

  # Log predictive densities, pooled like the forecasts: the log of the mean of
  # exp() of the stacked draws is the log of the average density of the members
  if (same_variables) {
    scores <- lapply(parts, function(x) x[["posterior"]][["forecast"]][["loglik"]])
    if (all(!vapply(scores, is.null, logical(1)))) {
      scores <- lapply(scores, as.matrix)
      n_periods <- min(vapply(scores, ncol, integer(1)))
      n_scores <- min(vapply(scores, nrow, integer(1)))
      result[["posterior"]][["forecast"]][["loglik"]] <- coda::mcmc(do.call(rbind, lapply(scores, function(s) {
        s[sort(sample.int(nrow(s), n_scores)), seq_len(n_periods), drop = FALSE]
      })))
    }

    predictive <- lapply(parts, function(x) x[["predictive"]])
    if (all(!vapply(predictive, is.null, logical(1)))) {
      periods <- unique(vapply(predictive, function(p) as.numeric(p[["period"]]), numeric(1)))
      if (length(periods) == 1) {
        n_pd <- min(vapply(predictive, function(p) length(p[["loglik"]]), integer(1)))
        result[["predictive"]] <- list(
          loglik = unlist(lapply(predictive, function(p) p[["loglik"]][sort(sample.int(length(p[["loglik"]]), n_pd))])),
          period = predictive[[1]][["period"]])
      }
    }
  }

  class(result) <- c("poolwindow", "bvarmodel", "list")
  return(result)
}
