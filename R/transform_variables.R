#' Apply Transformations
#'
#' Transforms the columns of a time-series object using the seven
#' transformation codes of FRED-MD/FRED-QD.
#'
#' @param x a time-series object.
#' @param code a named integer vector of transformation codes, one of
#' \code{1:7} -- see 'Details'. Names must match columns of \code{x}; a column
#' with no entry is left untransformed (code \code{1}).
#'
#' @details The seven codes, applied column by column, are
#' \describe{
#'   \item{1}{no transformation, \eqn{x_t}}
#'   \item{2}{first difference, \eqn{\Delta x_t}}
#'   \item{3}{second difference, \eqn{\Delta^2 x_t}}
#'   \item{4}{logarithm, \eqn{\log(x_t)}}
#'   \item{5}{first difference of the logarithm, \eqn{\Delta \log(x_t)}}
#'   \item{6}{second difference of the logarithm, \eqn{\Delta^2 \log(x_t)}}
#'   \item{7}{first difference of the growth rate,
#'   \eqn{\Delta(x_t / x_{t - 1} - 1)}}
#' }
#' Codes \code{4:6} need a strictly positive series.
#'
#' A difference drops as many leading observations as its order -- one for
#' codes \code{2}, \code{5} and \code{7}, two for \code{3} and \code{6} --
#' rather than shortening the series: those leading periods become \code{NA},
#' so every column keeps the time index of \code{x} and can still be combined
#' with \code{\link[stats]{ts.intersect}} or passed to
#' \code{\link{create_bvarmodel}}.
#'
#' @return A time-series object of the same shape as \code{x}.
#'
#' @examples
#'
#' data("us_macrodata")
#'
#' code <- c(Dp = 2, u = 1, r = 2)
#' x <- transform_variables(us_macrodata, code)
#'
#' @references
#'
#' McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for
#' macroeconomic research. \emph{Journal of Business & Economic Statistics,
#' 34}(4), 574--589.
#'
#' McCracken, M. W., & Ng, S. (2021). FRED-QD: A quarterly database for
#' macroeconomic research. \emph{Federal Reserve Bank of St. Louis Review,
#' 103}(1), 1--44.
#'
#' @export
transform_variables <- function(x, code) {

  # Input checks ----
  if (!"ts" %in% class(x)) {
    stop("Argument 'x' must be an object of class 'ts'.")
  }

  mat <- as.matrix(x)
  names_x <- dimnames(mat)[[2]]
  if (is.null(names_x)) {
    if (NCOL(mat) > 1) {
      stop("Argument 'x' must have column names when it has more than one ",
           "series.")
    }
    names_x <- NA_character_
  }

  if (!is.numeric(code)) {
    stop("Argument 'code' must be a named numeric vector of transformation ",
         "codes 1 to 7.")
  }
  if (any(is.na(code)) || any(code != as.integer(code)) ||
      any(!code %in% 1:7)) {
    stop("Argument 'code' must only contain the integers 1 to 7. See ",
         "'Details'.")
  }
  code_names <- names(code)
  if (NCOL(mat) > 1) {
    if (is.null(code_names) || any(code_names == "")) {
      stop("Argument 'code' must be named after the columns of 'x'.")
    }
    if (anyDuplicated(code_names) > 0) {
      stop("Argument 'code' has duplicate names: ",
           paste0("'", unique(code_names[duplicated(code_names)]), "'",
                   collapse = ", "), ".")
    }
    unknown <- setdiff(code_names, names_x)
    if (length(unknown) > 0) {
      stop("Argument 'code' names not found among the columns of 'x': ",
           paste0("'", unknown, "'", collapse = ", "), ".")
    }
  }

  # Transformation ----
  tsp_x <- stats::tsp(x)
  result <- mat
  for (j in seq_len(NCOL(mat))) {
    code_j <- if (NCOL(mat) == 1) {
      if (is.null(code_names)) code[[1]] else code[[names_x[j]]]
    } else {
      if (names_x[j] %in% code_names) code[[names_x[j]]] else 1L
    }
    result[, j] <- .transform_variables_code(result[, j], code_j)
  }

  result <- stats::ts(result, start = tsp_x[1], frequency = tsp_x[3])
  dimnames(result)[[2]] <- dimnames(mat)[[2]]
  result
}

# One column, transformed by one of the seven FRED-MD/FRED-QD codes. Leading
# observations lost to differencing come back as NA rather than shortening the
# series -- see 'Details' of transform_variables.
.transform_variables_code <- function(x, code) {
  n <- length(x)
  out <- rep(NA_real_, n)
  if (code == 1) {
    out <- x
  } else if (code == 2) {
    out[-1] <- diff(x, differences = 1)
  } else if (code == 3) {
    out[-(1:2)] <- diff(x, differences = 2)
  } else if (code == 4) {
    out <- log(x)
  } else if (code == 5) {
    out[-1] <- diff(log(x), differences = 1)
  } else if (code == 6) {
    out[-(1:2)] <- diff(log(x), differences = 2)
  } else if (code == 7) {
    growth <- x[-1] / x[-n] - 1
    out[-(1:2)] <- diff(growth, differences = 1)
  }
  out
}
