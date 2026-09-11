

# Names of the requested columns of the regressor matrix of a model, or NULL if
# the matrix does not carry names for them. Used wherever a model specification
# does not name a block of coefficients itself.

.regressor_columns <- function(object, pos) {

  x_names <- dimnames(object[["data"]][["train"]][["x"]])[[2]]
  if (is.null(x_names) || any(pos > length(x_names))) {
    return(NULL)
  }

  x_names[pos]
}


# Flattens the blocks of a model into one vector of regressor names, which is
# what the summary methods label their columns with.
# add_block prefixes every name with the symbol of the block it belongs to.

.flatten_regressor_blocks <- function(blocks, add_block = FALSE) {

  if (length(blocks) == 0) {
    return(NULL)
  }

  unlist(lapply(names(blocks), function(i) {
    if (add_block) {
      paste0(i, "\n", blocks[[i]][["labels"]])
    } else {
      blocks[[i]][["labels"]]
    }
  }), use.names = FALSE)
}


# Blocks of coefficients of a 'bvarmodel' object, in the order in which they
# appear among the regressors: the lags of the endogenous variables, the
# exogenous variables, the deterministic terms and the contemporaneous
# endogenous variables of a structural model.
#
# Each element carries the symbol the block is known by as its name, a title for
# a reader, and one label per regressor of the block.

.get_regressor_blocks_bvarmodel <- function(object) {

  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n"]]
  y_names <- object[["model"]][["endogen"]]
  blocks <- list()

  if (p > 0) {
    temp_names <- NULL
    for (i in 1:p) {
        temp_names <- c(temp_names, paste(y_names, ".l", i, sep = ""))
    }
    blocks[["A"]] <- list(title = "Lagged endogenous variables",
                          labels = temp_names)
  }

  if (m > 0) {
    exogen_names <- object[["model"]][["exogen"]]
    if (length(exogen_names) == m) {
      temp_names <- paste0(exogen_names, ".l0")
      if (s > 0) {
        temp_names <- c(temp_names, paste0(exogen_names, ".l", rep(1:s, each = m)))
      }
    } else {
      # A model whose specification does not name one variable per exogenous
      # block -- a GVAR sub-model, whose exogenous block already contains the
      # lags of the weakly exogenous and the global variables -- still carries
      # the names on the columns of its regressor matrix.
      temp_names <- .regressor_columns(object, k * p + 1:(m * (s + 1)))
      if (is.null(temp_names)) {
        exogen_names <- paste0("x", 1:m)
        temp_names <- paste0(exogen_names, ".l0")
        if (s > 0) {
          temp_names <- c(temp_names, paste0(exogen_names, ".l", rep(1:s, each = m)))
        }
      }
    }
    blocks[["B"]] <- list(title = "Exogenous variables", labels = temp_names)
  }

  if (n > 0) {
    temp_names <- object[["model"]][["deterministic"]]
    # A model that does not carry the names still has the terms, so take them
    # from the regressor matrix and fall back to placeholders only after that,
    # rather than return a vector too short for the coefficients it labels.
    if (length(temp_names) != n) {
      temp_names <- .regressor_columns(object, k * p + m * (s + 1) + 1:n)
    }
    if (length(temp_names) != n) {
      temp_names <- paste0("det.", 1:n)
    }
    blocks[["C"]] <- list(title = "Deterministic terms", labels = temp_names)
  }

  if (object[["model"]][["structural"]]) {
    blocks[["A0"]] <- list(title = "Contemporaneous endogenous variables",
                           labels = y_names)
  }

  return(blocks)
}


# Extracts the names of the regressors from a 'bvarmodel' object
# add_block adds the letter of the block of endogenous, exogensous, deterministic, structural and sigma coefficients

.get_regressor_names_bvarmodel <- function(object, add_block = FALSE) {

  .flatten_regressor_blocks(.get_regressor_blocks_bvarmodel(object),
                            add_block = add_block)
}


# Blocks of coefficients of a 'bvecmodel' object, in the order in which they
# appear among the regressors.
#
# The regressors of the error correction term are those of the cointegration
# matrix Pi and not those of the loading matrix alpha, since the draws of the
# former are what is reported for a VEC model.

.get_regressor_blocks_bvecmodel <- function(object) {

  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n"]]
  rank <- object[["model"]][["rank"]]
  blocks <- list()

  named <- function(names, count, fallback) {
    if (length(names) != count) {
      names <- fallback
    }
    return(names)
  }

  if (rank > 0) {
    ect_names <- dimnames(object[["data"]][["train"]][["w"]])[[2]]
    k_beta <- k + m + object[["model"]][["n_restricted"]]
    fallback <- c(paste0("l.", object[["model"]][["endogen"]]),
                  if (m > 0) paste0("l.", object[["model"]][["exogen"]]),
                  if (object[["model"]][["n_restricted"]] > 0) paste0("l.d", 1:object[["model"]][["n_restricted"]]))
    blocks[["Pi"]] <- list(title = "Cointegration matrix",
                           labels = named(ect_names, k_beta, fallback))
  }

  reg_names <- dimnames(object[["data"]][["train"]][["x"]])[[2]]
  pos <- 0

  n_gamma <- k * (p - 1)
  if (n_gamma > 0) {
    fallback <- paste0("d.", rep(object[["model"]][["endogen"]], times = p - 1),
                       ".l", rep(.lag_label(1:(p - 1), p - 1), each = k))
    blocks[["Gamma"]] <- list(title = "Lagged differenced endogenous variables",
                              labels = named(reg_names[pos + 1:n_gamma], n_gamma, fallback))
    pos <- pos + n_gamma
  }

  n_upsilon <- m * s
  if (n_upsilon > 0) {
    fallback <- paste0("d.", rep(object[["model"]][["exogen"]], times = s),
                       ".l", rep(.lag_label(0:(s - 1), s - 1), each = m))
    blocks[["Upsilon"]] <- list(title = "Differenced exogenous variables",
                                labels = named(reg_names[pos + 1:n_upsilon], n_upsilon, fallback))
    pos <- pos + n_upsilon
  }

  if (n > 0) {
    blocks[["C"]] <- list(title = "Unrestricted deterministic terms",
                          labels = named(reg_names[pos + 1:n], n, paste0("det.", 1:n)))
  }

  if (object[["model"]][["structural"]]) {
    blocks[["A0"]] <- list(title = "Contemporaneous endogenous variables",
                           labels = named(object[["model"]][["endogen"]], k,
                                          paste0("y", 1:k)))
  }

  return(blocks)
}


# Extracts the names of the regressors from a 'bvecmodel' object
# add_block adds the name of the block of cointegration, endogenous, exogenous,
# deterministic and structural coefficients

.get_regressor_names_bvecmodel <- function(object, add_block = FALSE) {

  .flatten_regressor_blocks(.get_regressor_blocks_bvecmodel(object),
                            add_block = add_block)
}
