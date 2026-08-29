#' @include selection_criteria.R
#' 
#' @param x an object used to select a method. Usually, the result of a call to
#' \code{\link{selection_criteria}}.
#' @param relative an integer specifying the model that is used as the reference
#' of relative forecast performance. Default is `0`, which indicates that results
#' are not displayed in relation to each other.
#' 
#'
#' @export
#' @rdname selection_criteria
print.selcritlist <- function(x, digits = max(3L, getOption("digits") - 3L), relative = 0, ...){
  
  n_models <- length(x)
  use_ll <- !is.null(x[[1]][["LL"]])
  use_fe <- !is.null(x[[1]][["FE"]])
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
    
    cat("\nLog-likelihood\n\n")
    ll <- template
    for (i in 1:n_models) {
      ll[i, 2:5] <- as.matrix(x[[i]][["LL"]])[1, ]
    }
    print(ll, digits = digits, row.names = FALSE, ...)
    
    
    cat("\n\nAkaike Information Criterion (AIC)\n\n")
    aic <- template
    for (i in 1:n_models) {
      aic[i, 2:5] <- as.matrix(x[[i]][["AIC"]])[1, ]
    }
    print(aic, digits = digits, row.names = FALSE, ...)
    
    
    cat("\n\nBayesian Information Criterion (BIC)\n\n")
    bic <- template
    for (i in 1:n_models) {
      bic[i, 2:5] <- as.matrix(x[[i]][["BIC"]])[1, ]
    }
    print(bic, digits = digits, row.names = FALSE, ...)
    
    
    cat("\n\nHannan-Quinn Criterion (HQ)\n\n")
    hq <- template
    for (i in 1:n_models) {
      hq[i, 2:5] <- as.matrix(x[[i]][["HQ"]])[1, ]
    }
    print(hq, digits = digits, row.names = FALSE, ...)
    
    
  }
  
  if (use_fe) {
    
    cat("\n\n------------------------------------------\n")
    cat("Out-of-sample")
    cat("\n------------------------------------------\n")
    
    n_rows <- nrow(x[[1]][["FE"]])
    varnames <- unique(unlist(lapply(x, function(y) {unique(y[["FE"]][, "variable"])})))
    
    template <- as.data.frame(matrix(NA, n_rows, 2 + n_models))
    names(template) <- c("Variable", "h", paste0("Model ", 1:n_models))
    template[, "Variable"] <- varnames
    template[, "h"] <- x[[1]][["FE"]][, "h"]
    
    pair_main <- paste0(template[, "Variable"], template[, "h"])
    
    cat("\nMean absolute forecast errors (MAFE)\n\n")
    
    mafe <- template
    for (i in 1:n_models) {
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
      pair_right <- paste0(x[[i]][["AFE"]][, "variable"], x[[i]][["AFE"]][, "h"])
      pos <- match(pair_right, pair_main)
      rsfe[pos, 2 + i] <- x[[i]][["RSFE"]][, "mean"]
    }
    if (relative > 0) {
      rsfe[, 3:ncol(rsfe)] <- rsfe[, 3:ncol(rsfe)] / rsfe[, 2 + relative]
    }
    print(rsfe, digits = digits, row.names = FALSE, ...)
    
  }
  
} 
