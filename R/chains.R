#' Convergence Diagnostics of Several Chains
#'
#' Compares the chains of a model estimated with \code{chains} above one in
#' \code{\link{add_posterior_coefficients}} and reports, for every parameter, the
#' split potential scale reduction factor and the effective sample size.
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
#' The statistic is the split \eqn{\hat{R}} of Gelman et al. (2013, section 11.4):
#' every chain is cut into two halves, and the variance between the means of the
#' halves is compared with the variance within them,
#' \deqn{\hat{R} = \sqrt{\frac{\frac{n - 1}{n} W + \frac{1}{n} B}{W}},}
#' where \eqn{n} is the length of a half, \eqn{W} the average variance within the
#' halves and \eqn{B} \eqn{n} times the variance of their means. Splitting the chains
#' also catches a chain that is still drifting. Values close to one indicate that
#' the chains describe the same distribution; above 1.01, Vehtari et al. (2021)
#' recommend running the chains longer or reconsidering the model. A parameter that
#' does not vary in any chain, such as a coefficient that variable selection
#' excluded throughout, has no \eqn{\hat{R}} and is reported as \code{NA}.
#'
#' The signed standard deviations \code{omega} of a block estimated under the
#' non-centred prior \code{omega_v} are symmetric around zero by construction --
#' the sampler switches their sign at random -- so their \eqn{\hat{R}} is close to
#' one whether or not the chains agree. Their squares, \code{sigma}, which are in
#' the posterior beside them, are the ones to read.
#'
#' The effective sample size is the sum over the chains of
#' \code{\link[coda]{effectiveSize}} of each chain. It is not a convergence
#' diagnostic: chains stuck in different modes can each have a large one.
#'
#' The chains start from the same initial values, those of
#' \code{\link{add_initial_values}}, and differ in their random numbers. Dispersed
#' starting values would make the comparison more demanding.
#'
#' @return A data frame with one row per parameter of every block of posterior draws,
#' other than the forecasts and the log-likelihood, and the columns \code{block},
#' \code{parameter}, \code{rhat} and \code{ess}.
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
    data.frame(block = name, parameter = parameter,
               rhat = .split_rhat(draws, chains),
               ess = .chain_ess(draws, chains),
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

# Split R-hat, column by column: each chain cut into two halves of equal length.
.split_rhat <- function(draws, chains) {
  n <- nrow(draws) %/% chains
  half <- n %/% 2
  rows <- unlist(lapply(seq_len(chains), function(c) {
    start <- (c - 1) * n + (n - 2 * half)
    list(start + seq_len(half), start + half + seq_len(half))
  }), recursive = FALSE)
  means <- sapply(rows, function(r) colMeans(draws[r, , drop = FALSE]))
  vars <- sapply(rows, function(r) apply(draws[r, , drop = FALSE], 2, stats::var))
  if (is.null(dim(means))) {
    means <- matrix(means, nrow = 1)
    vars <- matrix(vars, nrow = 1)
  }
  w <- rowMeans(vars)
  b <- half * apply(means, 1, stats::var)
  rhat <- sqrt(((half - 1) / half * w + b / half) / w)
  rhat[!is.finite(rhat) | w <= 0] <- NA_real_
  rhat
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
  list(chains = as.integer(chains), n = sum(!is.na(diag[["rhat"]])),
       above = sum(diag[["rhat"]] > 1.01, na.rm = TRUE),
       max_rhat = diag[["rhat"]][worst],
       worst = paste0(diag[["block"]][worst], "[", diag[["parameter"]][worst], "]"))
}

.print_chains_summary <- function(x) {
  if (is.null(x)) {
    return(invisible(NULL))
  }
  cat("\nChains:", x[["chains"]])
  if (!is.null(x[["max_rhat"]])) {
    cat(sprintf(". Largest split R-hat %.3f (%s); %d of %d parameters above 1.01.",
                x[["max_rhat"]], x[["worst"]], x[["above"]], x[["n"]]))
    if (x[["above"]] > 0) {
      cat("\nThe chains disagree: run them longer or reconsider the model before using",
          "the draws. See ?chain_diagnostics.")
    }
  }
  cat("\n")
  invisible(NULL)
}

.chain_ess <- function(draws, chains) {
  n <- nrow(draws) %/% chains
  ess <- 0
  for (c in seq_len(chains)) {
    chain <- draws[(c - 1) * n + seq_len(n), , drop = FALSE]
    ess <- ess + coda::effectiveSize(coda::mcmc(chain))
  }
  unname(ess)
}
