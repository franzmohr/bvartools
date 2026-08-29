#' @include selection_criteria.R
#'
#' @export
#' @rdname selection_criteria
print.selcrit <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  use_ll <- !is.null(x[["LL"]])
  use_fe <- !is.null(x[["FE"]])
  ci <- attr(x, "ci")
  
  
  if (use_ll) {
    
    cat("\n------------------------------------------\n")
    cat("In-sample")
    cat("\n------------------------------------------\n\n")
    
    result <- as.data.frame(matrix(NA, 4, 5))
    names(result) <- c("Criterion", "Mean", "Median",
                       paste0("Quantile (", ci[1], ")"),
                       paste0("Quantile (", ci[2], ")"))
    
    
    criterion <- c("LL", "AIC", "HQ", "BIC")
    for (i in 1:length(criterion)) {
      result[i, 1] <- criterion[i]
      result[i, 2:5] <- x[[criterion[i]]][, c("mean", "median", "qlower", "qupper")]
    }
    
    print(result, digits = digits, row.names = FALSE, ...)
    
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