#' @include selection_criteria.R
#'
#' @param x an object used to select a method. Usually, the result of a call to
#' \code{\link{selection_criteria}}.
#' @param digits the minimum number of significant digits to be printed in values.
#' @param relative an integer specifying the model that is used as the reference
#' of relative forecast performance. Default is `0`, which indicates that results
#' are not displayed in relation to each other.
#'
#'
#' @export
#' @rdname selection_criteria
print.selcritlist <- function(x, digits = max(3L, getOption("digits") - 3L), relative = 0, ...){

  n_models <- length(x)
  # Not every model has to contain both types of criteria. External forecasts, for
  # example, are not estimated, so they do not have in-sample criteria
  use_ll <- any(unlist(lapply(x, function(y) {!is.null(y[["LL"]])})))
  use_fe <- any(unlist(lapply(x, function(y) {!is.null(y[["FE"]])})))
  ci <- attr(x[[1]], "ci")

  if (use_ll) {

    cat("\n------------------------------------------\n")
    cat("In-sample")
    cat("\n------------------------------------------\n")

    template <- as.data.frame(matrix(NA, n_models, 5))
    names(template) <- c("", "Mean", "Median",
                         paste0("Quantile (", ci[1], ")"),
                         paste0("Quantile (", ci[2], ")"))
    template[, 1] <- paste0("Model ", 1:n_models)

    criteria <- c("LL", "AIC", "BIC", "HQ")
    headers <- c("Log-likelihood",
                 "Akaike Information Criterion (AIC)",
                 "Bayesian Information Criterion (BIC)",
                 "Hannan-Quinn Criterion (HQ)")

    for (j in 1:length(criteria)) {

      cat(ifelse(j == 1, "\n", "\n\n"), headers[j], "\n\n", sep = "")

      temp <- template
      for (i in 1:n_models) {
        # Models, which do not contain the criterion, remain missing
        if (!is.null(x[[i]][[criteria[j]]])) {
          temp[i, 2:5] <- as.matrix(x[[i]][[criteria[j]]])[1, ]
        }
      }
      print(temp, digits = digits, row.names = FALSE, ...)

    }

  }

  if (use_fe) {

    cat("\n\n------------------------------------------\n")
    cat("Out-of-sample")
    cat("\n------------------------------------------\n")

    # The models do not have to contain the same variables and forecast horizons,
    # so the union of both is used
    pairs <- do.call("rbind", lapply(x, function(y) {
      if (is.null(y[["FE"]])) {
        return(NULL)
      }
      y[["FE"]][, c("variable", "h")]
    }))
    pairs <- unique(pairs)

    template <- as.data.frame(matrix(NA, nrow(pairs), 2 + n_models))
    names(template) <- c("Variable", "h", paste0("Model ", 1:n_models))
    template[, "Variable"] <- pairs[, "variable"]
    template[, "h"] <- pairs[, "h"]

    pair_main <- paste0(template[, "Variable"], template[, "h"])

    cat("\nMean absolute forecast errors (MAFE)\n\n")

    mafe <- template
    for (i in 1:n_models) {
      if (is.null(x[[i]][["AFE"]])) {
        next
      }
      pair_right <- paste0(x[[i]][["AFE"]][, "variable"], x[[i]][["AFE"]][, "h"])
      pos <- match(pair_right, pair_main)
      mafe[pos, 2 + i] <- x[[i]][["AFE"]][, "mean"]
    }
    if (relative > 0) {
      mafe[, 3:ncol(mafe)] <- mafe[, 3:ncol(mafe)] / mafe[, 2 + relative]
    }
    print(mafe, digits = digits, row.names = FALSE, ...)


    cat("\n\nRoot mean squared forecast errors (RMSFE)\n\n")
    rsfe <- template
    for (i in 1:n_models) {
      if (is.null(x[[i]][["RSFE"]])) {
        next
      }
      pair_right <- paste0(x[[i]][["RSFE"]][, "variable"], x[[i]][["RSFE"]][, "h"])
      pos <- match(pair_right, pair_main)
      rsfe[pos, 2 + i] <- x[[i]][["RSFE"]][, "mean"]
    }
    if (relative > 0) {
      rsfe[, 3:ncol(rsfe)] <- rsfe[, 3:ncol(rsfe)] / rsfe[, 2 + relative]
    }
    print(rsfe, digits = digits, row.names = FALSE, ...)

  }

}
