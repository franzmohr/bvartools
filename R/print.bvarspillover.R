#' Print a Spillover Index
#'
#' Prints the connectedness table of Diebold and Yilmaz (2012).
#'
#' @param x an object of class 'bvarspillover', usually the result of a call to
#' \code{\link{spillover}}.
#' @param digits the number of decimal places. Defaults to 1, which is what this
#' literature reports.
#' @param ... not used.
#'
#' @details The body of the table is the posterior mean share of the forecast
#' error variance of the variable in the row that is attributed to a shock to the
#' variable in the column, in percent. The \code{from} column and the \code{to}
#' row are the directional spillovers, and the corner is the total index.
#'
#' Only the means appear here. The credible intervals are in the \code{total},
#' \code{from}, \code{to} and \code{net} elements of the object.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.bvarspillover <- function(x, digits = 1, ...) {

  spec <- x[["specification"]]

  cat("Spillover index (Diebold and Yilmaz, 2012)\n\n")
  cat("Horizon:      ", spec[["n_ahead"]], "\n")
  cat("Decomposition:",
      if (spec[["type"]] == "gir") "generalised" else "orthogonalised", "\n")
  if (!is.null(spec[["period"]])) {
    cat("Period:       ", spec[["period"]], "\n")
  }
  cat("Draws:        ", spec[["draws"]], "\n\n")

  # The table is already in percent and already carries its margins and its
  # corner, so there is nothing to do here but round it.
  print(round(x[["table"]], digits))

  cat("\nRows are responses, columns are shocks. The corner is the total index.\n")
  cat("Means only; see $total, $from, $to and $net for credible intervals.\n")

  invisible(x)
}
