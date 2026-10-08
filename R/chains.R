
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
