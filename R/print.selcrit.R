#' @include selection_criteria.R
#'
#' @export
#' @rdname selection_criteria
print.selcrit <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  use_ll <- !is.null(x[["LL"]]) || !is.null(x[["LML"]])
  use_fe <- !is.null(x[["FE"]])
  ci <- attr(x, "ci")
  
  
  if (use_ll) {
    
    cat("\n------------------------------------------\n")
    cat("In-sample")
    cat("\n------------------------------------------\n\n")
    
    # Only the criteria the object actually carries. WAIC needs more than one
    # draw to estimate the variance it penalises with, so it can be absent.
    # LML is the discounted models' one in-sample criterion and stands where LL
    # stands for the samplers, without the deviance-based criteria beside it:
    # there is no chain to estimate an effective number of parameters from, and
    # the marginal likelihood has already paid for the complexity a count of
    # parameters would charge for.
    criterion <- c("LL", "LML", "AIC", "HQ", "BIC", "WAIC", "LOOIC")
    criterion <- criterion[!vapply(x[criterion], is.null, logical(1))]

    result <- as.data.frame(matrix(NA, length(criterion), 5))
    names(result) <- c("Criterion", "Mean", "Median",
                       paste0("Quantile (", ci[1], ")"),
                       paste0("Quantile (", ci[2], ")"))

    for (i in 1:length(criterion)) {
      result[i, 1] <- criterion[i]
      result[i, 2:5] <- x[[criterion[i]]][, c("mean", "median", "qlower", "qupper")]
    }
    
    .print_criteria_table(result, digits = digits, ...)

    .print_waic_diagnostics(x[["WAIC"]])
    .print_loo_diagnostics(x[["LOOIC"]])

  }

  if (!is.null(x[["LPL"]])) {

    cat("\n\n------------------------------------------\n")
    cat("Predictive")
    cat("\n------------------------------------------\n\n")

    result <- data.frame(Criterion = "LPL", Mean = x[["LPL"]][["mean"]],
                         qlower = x[["LPL"]][["qlower"]], qupper = x[["LPL"]][["qupper"]])
    names(result)[3:4] <- paste0("Quantile (", ci, ")")
    .print_criteria_table(result, digits = digits, ...)

    terms <- attr(x[["LPL"]], "terms")
    cat("\nOne-step-ahead log predictive likelihood over ", nrow(terms), " periods; ",
        "numerical standard error ", format(attr(x[["LPL"]], "nse"), digits = digits), ".\n",
        sep = "")

  }

  if (use_fe) {
    
    cat("\n\n------------------------------------------\n")
    cat("Out-of-sample")
    cat("\n------------------------------------------\n")
    
    result <- as.data.frame(matrix(NA, nrow(x[["FE"]]), 4))
    result[, 1:2] <- x[["FE"]][, 1:2]
    names(result) <- c("Variable", "h", "MAFE", "RMSFE")
    result[, "MAFE"] <- x[["AFE"]][, "mean"]
    result[, "RMSFE"] <- x[["RSFE"]][, "mean"]
    result <- result[order(result[, "h"]),]
    result <- result[order(result[, "Variable"]),]
    print(result, digits = digits, row.names = FALSE, ...)
    
  }
  
} 
# Criteria that are point estimates have no credible band, and 'NA' in that
# column reads as a value that could not be computed rather than as one that
# does not exist. The columns are formatted first so that the blanks do not
# cost the remaining numbers their alignment.
.print_criteria_table <- function(x, digits, ...) {

  out <- x
  for (j in seq_along(out)) {
    if (is.numeric(out[[j]])) {
      column <- format(out[[j]], digits = digits)
      column[is.na(out[[j]])] <- ""
      out[[j]] <- format(column, justify = "right")
    }
  }

  print(out, row.names = FALSE, ...)
}
