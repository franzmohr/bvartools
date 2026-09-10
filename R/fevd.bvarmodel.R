#' Forecast Error Variance Decomposition
#' 
#' Produces the forecast error variance decomposition for an object of class 'bvarmodel'.
#' 
#' @param x an object of class 'bvarmodel'.
#' @param response name of the response variable.
#' @param n_ahead number of steps ahead.
#' @param type type of the impulse responses used to calculate forecast error variable decompositions.
#' Possible choices are orthogonalised \code{oir} (default) and generalised \code{gir} impulse responses.
#' @param normalise_gir logical. Should the GIR-based FEVD be normalised?
#' @param period integer. Index of the period, for which the variance decomposition should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param max_groups integer. Maximum number of variables the decomposition should contain.
#' The \code{max_groups - 1} variables with the largest contributions across the whole horizon
#' are kept and the contributions of the remaining variables are added up in a further column
#' named \code{"Other"}. This keeps the legend of the corresponding plot readable for models
#' with many variables. Default is \code{NULL}, so that a column is returned for every variable.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The function produces forecast error variance decompositions (FEVD) for the VAR model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}. For non-structural models matrix \eqn{A_0} is set to the identiy matrix
#' and can therefore be omitted, where not relevant.
#' 
#' If the FEVD is based on the orthogonalised impulse resonse (OIR), the FEVD will be calculated as
#' \deqn{\omega^{OIR}_{jk, h} = \frac{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i P e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma \Phi_i^{\prime} e_j )},}
#' where \eqn{\Phi_i} is the forecast error impulse response for the \eqn{i}th period,
#' \eqn{P} is the lower triangular Choleski decomposition of the variance-covariance
#' matrix \eqn{\Sigma}, \eqn{e_j} is a selection vector for the response variable and
#' \eqn{e_k} a selection vector for the impulse variable.
#'
#' If \code{type = "sir"}, the structural FEVD will be
#' calculated as \deqn{\omega^{SIR}_{jk, h} = \frac{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} A_0^{-1\prime} \Phi_i^{\prime} e_j )},}
#' where \eqn{\sigma_{jj}} is the diagonal element of the \eqn{j}th variable of the variance covariance matrix.
#'
#' If \code{type = "gir"}, the generalised FEVD will be
#' calculated as \deqn{\omega^{GIR}_{jk, h} = \frac{\sigma^{-1}_{jj} \sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i \Sigma \Phi_i^{\prime} e_j )},}
#' where \eqn{\sigma_{jj}} is the diagonal element of the \eqn{j}th variable of the variance covariance matrix.
#' 
#' If \code{type = "sgir"}, the structural generalised FEVD will be
#' calculated as \deqn{\omega^{SGIR}_{jk, h} = \frac{\sigma^{-1}_{jj} \sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} \Sigma e_k )^2}{\sum_{i = 0}^{h-1} (e_j^{\prime} \Phi_i A_0^{-1} \Sigma A_0^{-1\prime} \Phi_i^{\prime} e_j )}}.
#' 
#' Since GIR-based FEVDs do not add up to unity, they can be normalised by setting \code{normalise_gir = TRUE}.
#' 
#' @return A time-series object of class 'bvarfevd'.
#' 
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' e1 <- window(e1, end = c(1978, 4))
#' 
#' # Generate model data
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 100, burnin = 10)
#' # Chosen number of iterations and burnin should be much higher.
#' 
#' # Add prior specifications
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#' 
#' # Obtain posterior draws
#' object <- add_posterior_coefficients(model)
#' 
#' # Obtain FEVD
#' vd <- fevd(object, response = "cons")
#' 
#' # Plot FEVD
#' plot(vd)
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Pesaran, H. H., & Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' @export
fevd.bvarmodel <- function(x, response = NULL, n_ahead = 5, type = "oir", normalise_gir = FALSE, period = NULL,
                           max_groups = NULL, ...) {
  
  
  if (is.null(x[["posterior"]][["u_sigma_inv"]][["coeffs"]])) {
    stop("Argument 'object' must include draws of the variance-covariance matrix Sigma.")
  }
  
  if (!type %in% c("oir", "sir", "gir", "sgir")) {
    stop("The specified type of the used impulse response is not known.")
  }
  
  if(is.null(response)) {
    stop("Please provide a valid response variable.")
  }
  
  if (x[["model"]][["p"]] == 0 & !x[["model"]][["structural"]]) {
    stop("Variance decompositions only supported for models with p > 0 or structural models.")
  }

  need_A0 <- FALSE
  if (type %in% c("sgir", "sir")) {
    if (!x[["model"]][["structural"]]) {
      stop("Structural FEVD requires a structural model as input.")
    }
    need_A0 <- TRUE
  }

  max_groups <- .check_max_groups(max_groups)

  varnames <- x[["model"]][["endogen"]]
  response <- which(varnames == response)
  if (length(response) == 0){stop("Response variable not available.")}
  
  k <- x[["model"]][["k"]]

  # The draws in the shape .vardecomp wants them, shared with spillover() so
  # that the two cannot disagree about which slice of a row is `period`.
  A <- .collect_draws(x, period = period, need_A0 = need_A0)

  phi <- lapply(A, .vardecomp, h = n_ahead, type = type, response = response)
  
  result <- matrix(rowMeans(matrix(unlist(phi), (n_ahead + 1) * k)), n_ahead + 1)
  
  if (type %in% c("gir", "sgir")) {
    if (normalise_gir) {
      result <- t(apply(result, 1, function(x) {x / sum(x)}))
    }
  }
  
  colnames(result) <- varnames # Name columns

  if (!is.null(max_groups) && max_groups < k) {
    result <- .limit_fevd_groups(result, max_groups)
  }

  result <- stats::ts(result, start = 0, frequency = 1)
  
  class(result) <- append("bvarfevd", class(result))
  return(result)
}

# Reduce a variance decomposition to at most `max_groups` columns, so that the
# legend of the plot stays readable for a model with many variables.
#
# The `max_groups - 1` columns with the largest contribution across the whole
# horizon are kept in the order of the endogenous variables and the remaining
# ones are added up in a column "Other", which occupies the last of the slots.
# Row sums are left untouched, so a decomposition that added up to one before
# still does afterwards.
.limit_fevd_groups <- function(x, max_groups) {

  keep <- sort(order(colSums(x), decreasing = TRUE)[seq_len(max_groups - 1)])
  rest <- setdiff(seq_len(ncol(x)), keep)

  label <- "Other"
  while (label %in% colnames(x)[keep]) {
    label <- paste0(label, "_")
  }

  result <- cbind(x[, keep, drop = FALSE], rowSums(x[, rest, drop = FALSE]))
  colnames(result) <- c(colnames(x)[keep], label)

  return(result)
}

# Validate the `max_groups` argument of fevd() and plot.bvarfevd(), which offer
# it with the same meaning and should refuse the same input. Returns NULL, which
# stands for "show every variable", or the value as an integer.
.check_max_groups <- function(max_groups) {

  if (is.null(max_groups)) {
    return(NULL)
  }

  if (length(max_groups) != 1 || !is.numeric(max_groups) || is.na(max_groups) || max_groups < 1) {
    stop("Argument 'max_groups' must be a single positive integer.")
  }

  return(as.integer(max_groups))
}
