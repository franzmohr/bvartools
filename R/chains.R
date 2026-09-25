#' Convergence Diagnostics of Several Chains
#'
#' Compares the chains of a model estimated with \code{chains} above one in
#' \code{\link{add_posterior_coefficients}} and reports, for every parameter, the
#' potential scale reduction factor and two effective sample sizes.
#'
#' @param object an object of class 'bvarmodel' or 'bvecmodel', whose posterior was
#' simulated with \code{chains} of at least two.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' A single chain can look converged and not be: its autocorrelations and its
#' effective sample size describe how well it explores the region it is in, not
#' whether that region is the posterior. Chains started from the same values but
#' drawing different random numbers settle in different places if the posterior has
#' several modes or the sampler has not left the neighbourhood of its start, and
#' comparing them is the check a single chain cannot provide.
#'
#' The statistic is the rank-normalised split \eqn{\hat{R}} of Vehtari et al.
#' (2021). Every chain is cut into two halves, the pooled draws are replaced by
#' their ranks and put back on a normal scale, and the variance between the means
#' of the halves is compared with the variance within them,
#' \deqn{\hat{R} = \sqrt{\frac{\frac{n - 1}{n} W + \frac{1}{n} B}{W}},}
#' where \eqn{n} is the length of a half, \eqn{W} the average variance within the
#' halves and \eqn{B} \eqn{n} times the variance of their means. Splitting the chains
#' catches a chain that is still drifting.
#'
#' \strong{The rank normalisation is what makes the number mean the same thing
#' whatever scale the parameter is on.} The plain version of this statistic is
#' built on variances, so it is undefined for a posterior heavy-tailed enough to
#' have none and it moves when the same draws are transformed monotonically --
#' which is awkward for a package whose parameters are freely rescaled and whose
#' variance blocks are reported both as \code{omega} and as its square. What is
#' reported is the larger of the rank-normalised \eqn{\hat{R}} and the one computed
#' on the draws folded around their median, so that chains agreeing about the
#' centre but not about the spread are caught as well.
#'
#' Values close to one indicate that the chains describe the same distribution;
#' above 1.01, Vehtari et al. (2021) recommend running the chains longer or
#' reconsidering the model. A parameter that does not vary in any chain, such as a
#' coefficient that variable selection excluded throughout or an equation that
#' \code{iid} left without coefficients, has no \eqn{\hat{R}} and is reported as
#' \code{NA}.
#'
#' The signed standard deviations \code{omega} of a block estimated under the
#' non-centred prior \code{omega_v} are symmetric around zero by construction --
#' the sampler switches their sign at random -- so their \eqn{\hat{R}} is close to
#' one whether or not the chains agree. Their squares, \code{sigma}, which are in
#' the posterior beside them, are the ones to read.
#'
#' \strong{Two effective sample sizes are reported, and for this package the
#' second is usually the one that matters.} \code{ess_bulk} is computed on the
#' rank-normalised draws and says how much independent information the sample
#' carries about the centre of the posterior. \code{ess_tail} is the smaller of
#' the effective sample sizes at the 5th and 95th percentiles and says the same
#' about the extremes. Almost nothing this package reports is a point:
#' \code{\link{irf}}, \code{\link{fevd}} and the forecasts all come back as
#' quantiles, and it is \code{ess_tail} that says whether those quantiles have
#' settled. A sample can carry a comfortable \code{ess_bulk} and still have a
#' band that moves from one run to the next.
#'
#' Neither is a convergence diagnostic: chains stuck in different modes can each
#' have a large one.
#'
#' The chains start from the same initial values, those of
#' \code{\link{add_initial_values}}, and differ in their random numbers. Dispersed
#' starting values would make the comparison more demanding.
#'
#' @return A data frame with one row per parameter of every block of posterior draws,
#' other than the forecasts and the log-likelihood, and the columns \code{block},
#' \code{parameter}, \code{rhat}, \code{ess_bulk} and \code{ess_tail}.
#'
#' @references
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B.
#' (2013). \emph{Bayesian data analysis} (3rd ed.). Boca Raton: CRC Press.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021).
#' Rank-normalization, folding, and localization: An improved \eqn{\hat{R}} for
#' assessing convergence of MCMC. \emph{Bayesian Analysis, 16}(2), 667--718.
#' \doi{10.1214/20-BA1221}
#'
#' @examples
#' data("e1")
#' e1 <- diff(log(e1)) * 100
#'
#' model <- create_bvarmodel(e1, p = 1, deterministic = "const",
#'                           iterations = 200, burnin = 100)
#' model <- add_priors(model, coef = list(v_i = 0, v_i_det = 0),
#'                     sigma = list(df = "k", scale = 1))
#' model <- add_initial_values(model)
#' model <- add_posterior_coefficients(model, chains = 2)
#'
#' diag <- chain_diagnostics(model)
#' max(diag$rhat, na.rm = TRUE)
#'
#' @family posterior simulation
#' @export
chain_diagnostics <- function(object, ...) {
  if (!any(c("bvarmodel", "bvecmodel") %in% class(object))) {
    stop("Argument 'object' must be an object of class 'bvarmodel' or 'bvecmodel'.")
  }
  chains <- object[["model"]][["chains"]]
  if (is.null(chains) || as.integer(chains) < 2) {
    stop("The model was estimated with a single chain. Use add_posterior_coefficients() ",
         "with 'chains' of at least 2 to compare several.")
  }
  chains <- as.integer(chains)

  blocks <- .chain_blocks(object[["posterior"]])
  if (length(blocks) == 0) {
    stop("Argument 'object' does not contain posterior draws.")
  }

  result <- lapply(names(blocks), function(name) {
    draws <- .draws_matrix(blocks[[name]])
    n_total <- nrow(draws)
    if (n_total %% chains != 0) {
      stop("The ", n_total, " draws of '", name, "' do not divide into ", chains,
           " chains of equal length. Were they thinned by a factor that the length ",
           "of a chain is not a multiple of?")
    }
    n <- n_total %/% chains
    if (n < 4) {
      stop("Each chain has ", n, " draws, too few to split into halves.")
    }
    parameter <- colnames(draws)
    if (is.null(parameter)) {
      parameter <- as.character(seq_len(ncol(draws)))
    }
    stats <- vapply(seq_len(ncol(draws)), function(j) {
      # The pooled block stacks the chains one after another, so filling a
      # matrix column by column puts each chain in a column of its own, which
      # is the shape the diagnostics of the posterior package take.
      x <- matrix(draws[, j], nrow = n)
      # posterior warns when it caps an effective sample size that came out
      # above the number of draws, which an antithetic chain legitimately
      # produces. The cap is the right answer and there is nothing for the user
      # to do about it, so the warning is not passed on -- once per parameter it
      # would bury the diagnostics it is attached to.
      withCallingHandlers(
        c(posterior::rhat(x), posterior::ess_bulk(x), posterior::ess_tail(x)),
        warning = function(w) {
          if (grepl("capped", conditionMessage(w), fixed = TRUE)) {
            invokeRestart("muffleWarning")
          }
        }
      )
    }, numeric(3))

    data.frame(block = name, parameter = parameter,
               rhat = stats[1, ], ess_bulk = stats[2, ], ess_tail = stats[3, ],
               stringsAsFactors = FALSE, row.names = NULL)
  })
  do.call(rbind, result)
}

# The seed of each chain. The first is the model's own, so that one chain draws
# what it always has; the others are spread far from the seed + 1, seed + 2, ...
# that add_seed() gives the neighbouring models of a list.
.chain_seeds <- function(seed, chains) {
  offsets <- (seq_len(chains) - 1) * 1000003
  as.integer((as.numeric(seed) + offsets) %% .Machine$integer.max)
}

# The number of chains a call asks for: its argument, or what the model says.
.check_chains <- function(object, chains) {
  if (is.null(chains)) {
    chains <- object[["model"]][["chains"]]
  }
  if (is.null(chains)) {
    return(1L)
  }
  if (length(chains) != 1 || !is.numeric(chains) || !is.finite(chains) ||
      chains < 1 || abs(chains - round(chains)) > 1e-8) {
    stop("Argument 'chains' must be a single whole number of at least 1.")
  }
  chains <- as.integer(round(chains))
  if (chains > 1 && .is_discount(object)) {
    stop("A discounted model has a closed-form posterior rather than a chain, so there ",
         "are no chains to compare. Use 'chains = 1'.")
  }
  chains
}

# Draws the posterior of 'object' once per chain with 'draw' and pools the chains,
# one after the other, into the blocks every later step reads.
.run_chains <- function(object, chains, draw) {
  seed <- object[["model"]][["seed"]]
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }
  seeds <- .chain_seeds(seed, chains)
  fits <- lapply(seeds, function(s) {
    single <- object
    single[["model"]][["seed"]] <- s
    single[["model"]][["chains"]] <- NULL
    draw(single)
  })
  result <- fits[[1]]
  result[["model"]][["seed"]] <- as.integer(seed)
  result[["model"]][["chains"]] <- chains
  result[["posterior"]] <- .pool_chains(lapply(fits, function(x) x[["posterior"]]))
  result
}

# Stacks the draws of the chains block by block. The thinning interval is that of
# a chain, and the draws are labelled on through the chains.
.pool_chains <- function(posteriors) {
  first <- posteriors[[1]]
  if (coda::is.mcmc(first) || is.matrix(first)) {
    draws <- do.call(rbind, lapply(posteriors, .draws_matrix))
    mc <- if (coda::is.mcmc(first)) coda::mcpar(first) else c(1, NROW(first), 1)
    return(coda::mcmc(draws, start = mc[1], end = mc[1] + (nrow(draws) - 1) * mc[3], thin = mc[3]))
  }
  if (is.list(first)) {
    for (name in names(first)) {
      first[[name]] <- .pool_chains(lapply(posteriors, function(x) x[[name]]))
    }
  }
  first
}

# Every block of draws of a posterior other than the forecasts and the
# log-likelihood, flattened to "block$element" names.
.chain_blocks <- function(posterior, prefix = NULL) {
  blocks <- list()
  for (name in setdiff(names(posterior), c("forecast", "loglik"))) {
    element <- posterior[[name]]
    label <- paste(c(prefix, name), collapse = "$")
    if (coda::is.mcmc(element) || is.matrix(element)) {
      blocks[[label]] <- element
    } else if (is.list(element)) {
      blocks <- c(blocks, .chain_blocks(element, label))
    }
  }
  blocks
}

# What summary() reports about the chains: nothing for one chain.
.chains_summary <- function(object) {
  chains <- object[["model"]][["chains"]]
  if (is.null(chains) || as.integer(chains) < 2) {
    return(NULL)
  }
  diag <- tryCatch(chain_diagnostics(object), error = function(e) NULL)
  if (is.null(diag) || all(is.na(diag[["rhat"]]))) {
    return(list(chains = as.integer(chains)))
  }
  worst <- which.max(diag[["rhat"]])
  thinnest <- which.min(diag[["ess_tail"]])
  list(chains = as.integer(chains), n = sum(!is.na(diag[["rhat"]])),
       above = sum(diag[["rhat"]] > 1.01, na.rm = TRUE),
       max_rhat = diag[["rhat"]][worst],
       worst = paste0(diag[["block"]][worst], "[", diag[["parameter"]][worst], "]"),
       # The bands irf() and the forecasts report are quantiles, so the smallest
       # tail effective sample size is what says whether they have settled --
       # which the largest R-hat, about agreement rather than precision, does not.
       min_ess_tail = if (length(thinnest) == 1) diag[["ess_tail"]][thinnest] else NA_real_,
       thinnest = if (length(thinnest) == 1) {
         paste0(diag[["block"]][thinnest], "[", diag[["parameter"]][thinnest], "]")
       } else NA_character_)
}

.print_chains_summary <- function(x) {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  cat("\nChains:", x[["chains"]])
  if (!is.null(x[["max_rhat"]])) {
    cat(sprintf(". Largest R-hat %.3f (%s); %d of %d parameters above 1.01.",
                x[["max_rhat"]], x[["worst"]], x[["above"]], x[["n"]]))
    if (!is.null(x[["min_ess_tail"]]) && is.finite(x[["min_ess_tail"]])) {
      cat(sprintf("\nSmallest tail effective sample size %.0f (%s).",
                  x[["min_ess_tail"]], x[["thinnest"]]))
      if (x[["min_ess_tail"]] < 400) {
        cat("\nBelow about 400 the outer quantiles of a credible band are not settled,",
            "so the bands of irf(), fevd() and the forecasts will move between runs.")
      }
    }
    if (x[["above"]] > 0) {
      cat("\nThe chains disagree: run them longer or reconsider the model before using",
          "the draws. See ?chain_diagnostics.")
    }
  }
  cat("\n")
  invisible(NULL)
}
