#' Single Quantiles of a Quantile Grid
#'
#' Splits a structural quantile VAR, a model created with \code{quantile_grid = TRUE} in
#' \code{\link{create_bvarmodel}}, into one model per quantile.
#'
#' The draws of a grid are stacked level by level, and each level's draws are the chain the
#' single-quantile model would have drawn. The models returned hold those draws with the
#' specification of that model, so that everything written for a single quantile --
#' \code{\link[=summary.bvarmodel]{summary}}, \code{\link[=plot.bvarmodel]{plot}},
#' \code{\link{add_posterior_loglik}} -- reads them as such. What describes the grid as a
#' whole, its forecasts and its log likelihood, is not carried over. The impulse responses of
#' a structural grid are differences of forecasts at a fixed level with and without a
#' scenario; see \code{forecast_quantile} in \code{\link{add_posterior_forecasts}}.
#'
#' @param object an object of class 'bvarmodel' holding a grid of quantiles, usually the
#' result of a call to \code{\link{add_posterior_coefficients}}.
#'
#' @return A list of class 'modellist' with one object of class 'bvarmodel' per quantile, in
#' increasing order.
#'
#' @examples
#'
#' data("us_macrodata")
#'
#' object <- create_bvarmodel(data = us_macrodata, p = 1, deterministic = "const",
#'                            structural = TRUE, error = "ald",
#'                            quantile = c(0.1, 0.5, 0.9), quantile_grid = TRUE,
#'                            iterations = 20, burnin = 10)
#' object <- add_priors(object, coef = list(v_i = 1),
#'                      sigma = list(shape = 3, rate = .01))
#' object <- add_initial_values(object)
#' object <- add_posterior_coefficients(object)
#'
#' levels <- split_quantile_grid(object)
#' summary(levels[[3]])
#'
#' @seealso \code{\link{create_bvarmodel}}
#' @family post-estimation analysis
#' @export
split_quantile_grid <- function(object) {
  tau <- object[["model"]][["quantiles"]]
  if (is.null(tau)) {
    stop("Argument 'object' does not hold a grid of quantiles; see 'quantile_grid' in ",
         "create_bvarmodel().")
  }
  n_levels <- length(tau)

  # Columns j of a stacked block, relabelled as the chain they are.
  level_block <- function(draws, j) {
    if (is.null(draws)) {
      return(NULL)
    }
    n <- ncol(draws) / n_levels
    tsp_draws <- coda::mcpar(draws)
    coda::mcmc(as.matrix(draws)[, (j - 1) * n + seq_len(n), drop = FALSE],
               start = tsp_draws[1], end = tsp_draws[2], thin = tsp_draws[3])
  }

  result <- lapply(seq_len(n_levels), function(j) {
    level <- object
    level[["model"]][["quantiles"]] <- NULL
    level[["model"]][["forecast_quantile"]] <- NULL
    level[["model"]][["quantile"]] <- tau[j]
    posterior <- object[["posterior"]]
    if (!is.null(posterior)) {
      a <- .drop_null(list(coeffs = level_block(posterior[["a"]][["coeffs"]], j),
                           lambda = level_block(posterior[["a"]][["lambda"]], j)))
      posterior <- list(a = if (length(a) > 0) a else NULL,
                        u_scale = list(coeffs = level_block(posterior[["u_scale"]][["coeffs"]], j)))
      posterior <- .drop_null(posterior)
    }
    level[["posterior"]] <- posterior
    level
  })
  names(result) <- paste0("q", tau)
  class(result) <- c("modellist", "list")
  result
}

# The chain labels of a model's draws. A grid of quantiles keeps no error
# precisions, so its scales carry them instead.
.draws_mcpar <- function(object) {
  draws <- object[["posterior"]][["u_sigma_inv"]][["coeffs"]]
  if (is.null(draws)) {
    draws <- object[["posterior"]][["u_scale"]][["coeffs"]]
  }
  coda::mcpar(draws)
}

# The draws of a grid of quantiles are stacked level by level, which nothing
# written for a single set of coefficients can read.
.refuse_quantile_grid <- function(x, caller) {
  if (!is.null(x[["model"]][["quantiles"]])) {
    stop(caller, " of a grid of quantiles work on its single levels, which ",
         "split_quantile_grid() returns.", call. = FALSE)
  }
}
