#' Impulse Response Function
#' 
#' Computes the impulse response coefficients for an object of class 'bvarmodel'.
#' 
#' @param x an object of class 'bvarmodel'.
#' @param impulse name of the impulse variable.
#' @param response name of the response variable.
#' @param n_ahead number of steps ahead. Zero is allowed and returns the impact response
#' alone.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param shock size of the shock. For \code{"oir"}, \code{"gir"} and \code{"sgir"} it is counted
#' in standard deviations of the shock, so the default of 1 is a shock of one standard deviation.
#' For \code{"feir"} and \code{"sir"} it is counted in units of the reduced form or structural error
#' of the impulse variable, and for \code{"sign"} and \code{"custom"} it rescales the columns of the
#' impact matrix. \code{"sd"} and \code{"nsd"} are a positive and a negative shock of one standard
#' deviation and are not available for \code{"sign"} and \code{"custom"}.
#' @param type type of the impulse response. Possible choices are forecast error \code{"feir"}
#' (default), orthogonalised \code{"oir"}, structural \code{"sir"}, generalised \code{"gir"},
#' structural generalised \code{"sgir"}, sign restricted \code{"sign"} and \code{"custom"}
#' impulse responses. For a structural model only \code{"sir"} and \code{"sgir"} are available;
#' the other five require a non-structural model. See 'Details'.
#' @param cumulative logical specifying whether a cumulative IRF should be calculated.
#' @param keep_draws logical specifying whether the function should return all draws of
#' the posterior impulse response function. Defaults to \code{FALSE} so that
#' the median and the credible intervals of the posterior draws are returned.
#' @param period integer. Index of the period, for which the IR should be generated.
#' Only used for TVP or SV models. Default is \code{NULL}, so that the posterior draws of the last time period
#' are used.
#' @param impact the impact matrix of a \code{"custom"} impulse response, either a single
#' \eqn{K \times K} matrix that identifies every posterior draw the same way, or a list of
#' such matrices with one entry per draw. Ignored for every other value of \code{type}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details The function produces different types of impulse responses for the VAR model
#' \deqn{A_0 y_t = \sum_{i = 1}^{p} A_{i} y_{t-i} + u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}.
#' 
#' Forecast error impulse responses \eqn{\Phi_i} are obtained by recursions
#' \deqn{\Phi_i = \sum_{j = 1}^{i} \Phi_{i-j} A_j,   i = 1, 2,...,h}
#' with \eqn{\Phi_0 = I_K}.
#' 
#' Orthogonalised impulse responses \eqn{\Theta^o_i} are calculated as \eqn{\Theta^o_i = \Phi_i P},
#' where P is the lower triangular Choleski decomposition of \eqn{\Sigma}, so that they are the
#' responses to orthogonalised shocks of one standard deviation.
#' 
#' Structural impulse responses \eqn{\Theta^s_i} are calculated as \eqn{\Theta^s_i = \Phi_i A_0^{-1}}.
#' 
#' For a structural model the posterior draws describe the structural form: the coefficients are
#' the \eqn{A_i} of the equation above and \eqn{\Sigma} is the covariance of the structural errors.
#' The recursion for \eqn{\Phi_i} needs the reduced form, i.e. \eqn{A_0^{-1} A_i} and
#' \eqn{A_0^{-1} \Sigma A_0^{-1\prime}}, which only \code{"sir"} and \code{"sgir"} form. The other
#' types are therefore not available for a structural model, and the two structural types not for
#' any other.
#' 
#' Sign restricted impulse responses \eqn{\Theta^r_i} are calculated as
#' \eqn{\Theta^r_i = \Phi_i P Q}, where \eqn{P} is the lower triangular Choleski decomposition of
#' \eqn{\Sigma} and \eqn{Q} the rotation that \code{\link{add_sign_restrictions}} accepted for that
#' draw. Draws for which no admissible rotation was found are left out, so the responses cover
#' fewer draws than the posterior holds and the credible interval is one over the set of models the
#' restrictions admit rather than over a single identified model.
#'
#' Custom impulse responses \eqn{\Theta^c_i} are calculated as \eqn{\Theta^c_i = \Phi_i P}, where
#' \eqn{P} is the matrix supplied in argument \code{impact}. This is the route by which an
#' identification that the package does not derive itself reaches the recursion; \code{shock}
#' rescales the result but nothing normalises the columns of \eqn{P}, so an impact matrix that
#' means to deliver unit shocks has to arrive that way.
#' 
#' (Structural) Generalised impulse responses to a shock to variable \eqn{j} are calculated as
#' \eqn{\Theta^g_{i} = \sigma_{jj}^{-1/2} \Phi_i A_0^{-1} \Sigma e_j}, where \eqn{\sigma_{jj}} is the
#' \eqn{j}th diagonal element of \eqn{\Sigma}, the variance of the error of the impulse variable, and
#' \eqn{e_j} is a selection vector containing one in its \eqn{j}th element and zero otherwise. They are
#' therefore the responses to a shock of one standard deviation (Pesaran and Shin, 1998). If the
#' \code{"bvarmodel"} object does not contain draws of \eqn{A_0}, it is assumed to be an identity matrix.
#' 
#' @return A time-series object of class 'bvarirf' running from period 0 to \code{n_ahead},
#' with the lower bound, the median and the upper bound of the credible band of the
#' response in three columns named after their quantiles, e.g. \code{"2.5\%"},
#' \code{"50\%"} and \code{"97.5\%"} for \code{ci = .95}. If \code{keep_draws = TRUE}, a
#' matrix of class 'bvarirf' with one row per draw and one column per period instead.
#'
#' @examples
#' 
#' # Load data
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#' 
#' # Create model
#' model <- create_bvarmodel(e1, p = 2, deterministic = "const",
#'                           iterations = 20, burnin = 10)
#' # Number of iterations and burnin should be much higher.
#' 
#' # Add priors
#' model <- add_priors(model,
#'                     coef = list(v_i = 1, v_i_det = 1 / 10),
#'                     sigma = list(df = "k", scale = 1))
#' 
#' # Add initial values
#' model <- add_initial_values(model)
#'
#' # Obtain posterior draws 
#' model <- add_posterior_coefficients(model)
#' 
#' # Obtain IR
#' ir <- irf(model, impulse = "invest", response = "cons")
#' 
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' Pesaran, H. H., Shin, Y. (1998). Generalized impulse response analysis in linear multivariate models. \emph{Economics Letters, 58}, 17-29.
#' 
#' @family post-estimation analysis
#' @export
#' @method irf bvarmodel
irf.bvarmodel <- function(x, impulse = NULL, response = NULL, n_ahead = 5, ci = .95, shock = 1,
                          type = "feir", cumulative = FALSE, keep_draws = FALSE, period = NULL,
                          impact = NULL, ...) {
  
  if (!type %in% c("feir", "oir", "gir", "sir", "sgir", "sign", "custom")) {
    stop("Argument 'type' not known.")
  }

  if (type == "custom" && is.null(impact)) {
    stop("Impulse responses of type \"custom\" need an impact matrix in argument 'impact'.")
  }

  # A sign restricted identification is a custom one whose impact matrices the
  # model is already carrying, so it is assembled here and travels the same
  # path. `type` keeps its own name until then, so that the checks below and
  # any message they produce speak of what was asked for.
  if (type == "sign") {
    impact <- .sign_impact(x, "Impulse responses")
  }
  
  # A horizon of zero is the impact period on its own, which is well defined:
  # the response is Phi_0 P = P, with no recursion behind it. A negative
  # horizon is not, and reaches the C++ worker as a matrix of no rows.
  if (length(n_ahead) != 1 || !is.numeric(n_ahead) || is.na(n_ahead) || n_ahead < 0) {
    stop("Argument 'n_ahead' must be a single integer of at least 0.")
  }
  
  if (x[["model"]][["p"]] == 0 & !x[["model"]][["structural"]]) {
    stop("Impulse responses only supported for models with p > 0 or structural models.")
  }
  
  need_A0 <- FALSE
  if (type %in% c("sgir", "sir")) {
    if (!x[["model"]][["structural"]]) {
      stop("Structural IR requires a structural model as input.")
    }
    need_A0 <- TRUE
  } else {
    # A structural model stores the contemporaneous block separately, so its
    # coefficient draws are the structural A_i and its covariance draws the
    # covariance of the structural errors. The recursion behind these types
    # wants the reduced form -- A_0^{-1} A_i and A_0^{-1} Sigma A_0^{-1}' --
    # and reading the structural quantities in their place silently produces
    # responses that belong to no model at all.
    if (x[["model"]][["structural"]]) {
      stop("Impulse responses of type \"", type, "\" are not defined for a structural model: ",
           "they would be calculated from the structural coefficients and the covariance of the ",
           "structural errors instead of the reduced form the recursion needs. Use type \"sir\" ",
           "or \"sgir\" for a structural model, or estimate the model with 'structural = FALSE'.")
    }
  }
  
  if (!(is.numeric(shock) | shock %in% c("sd", "nsd"))) {
    stop("Invalid specification of argument 'shock'.")
  }

  # The standard deviation based sizes read a scale off Sigma, which is a
  # statement about the reduced form errors and not about the shocks a custom
  # impact matrix defines. Whatever scale those shocks have is already in the
  # columns of that matrix, so there is nothing here left to infer.
  if (type %in% c("custom", "sign") && !is.numeric(shock)) {
    stop("Argument 'shock' must be numeric for an impulse response of type \"", type, "\": ",
         "the size of a shock is carried by the impact matrix.")
  }
  
  # "sign" needs it to build the Choleski factor its rotation acts on.
  if (type %in% c("oir", "gir", "sgir", "sign") | shock %in% c("sd", "nsd")) {
    need_Sigma <- TRUE
  } else {
    need_Sigma <- FALSE
  }
  
  varnames <- x[["model"]][["endogen"]]
  impulse <- which(varnames == impulse)
  if (length(impulse) == 0){stop("Impulse variable not available.")}
  response <- which(varnames == response)
  if (length(response) == 0){stop("Response variable not available.")}
  
  # Shared with fevd() and spillover(), so that the three cannot disagree about
  # which slice of a row belongs to `period`, nor about folding a structural
  # model's contemporaneous block into the coefficients the recursion uses.
  A <- .collect_draws(x, period = period, need_A0 = need_A0,
                      need_Sigma = need_Sigma, impact = impact)

  # Size of the shock, one value per draw. A numeric shock is the same for every
  # draw. "oir", "gir" and "sgir" count a shock in standard deviations already,
  # so a one standard deviation shock is a shock of size one there. "feir" and
  # "sir" count it in units of the error, whose standard deviation is read off
  # that draw's Sigma.
  for (i in seq_along(A)) {
    if (is.numeric(shock)) {
      A[[i]][["shock"]] <- shock
    } else {
      if (type %in% c("oir", "gir", "sgir")) {
        A[[i]][["shock"]] <- 1
      } else {
        A[[i]][["shock"]] <- sqrt(diag(A[[i]][["Sigma"]])[impulse])
      }

      if (shock == "nsd") {
        A[[i]][["shock"]] <- -A[[i]][["shock"]]
      }
    }
  }

  result <- lapply(A, .ir, h = n_ahead, type = if (type == "sign") "custom" else type,
                   impulse = impulse, response = response)
  
  result <- t(matrix(unlist(result), n_ahead + 1))
  
  if (cumulative) {
    # apply() returns a plain vector when every row holds a single value, which
    # t() would turn into one row of draws instead of one column of horizons --
    # the shape at n_ahead = 0. Rebuild the draws x horizon matrix explicitly.
    result <- matrix(apply(result, 1, cumsum), nrow = nrow(result), byrow = TRUE)
  }
  
  if (!keep_draws) {
    result <- .summarise_irf_draws(result, ci)
  }
  
  class(result) <- append("bvarirf", class(result))
  return(result)
}

# The quantiles an impulse response reports, over the draws of the responses.
# Shared with the method for a stored model, which takes them over the draws of
# every piece of the chain together, so that the two cannot summarise the same
# responses differently.
.summarise_irf_draws <- function(result, ci) {

  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  pr <- c(ci_low, .5, ci_high)

  stats::ts(t(apply(result, 2, stats::quantile, probs = pr)), start = 0, frequency = 1)
}
