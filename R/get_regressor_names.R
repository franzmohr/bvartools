

# Extracts the names of the regressors from a 'bvarmodel' object
# add_block adds the letter of the block of endogenous, exogensous, deterministic, structural and sigma coefficients

.get_regressor_names_bvarmodel <- function(object, add_block = FALSE) {
  
  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n"]]
  tvp <- object[["model"]][["tvp"]]
  y_names <- dimnames(object[["data"]][["original"]][["endogen"]])[[2]]
  x_names <- NULL
  
  if (p > 0) {
    temp_names <- NULL
    for (i in 1:p) {
        temp_names <- c(temp_names, paste(y_names, ".l", i, sep = ""))
    } 
    if (add_block) {
      temp_names <- paste0("A\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  if (m > 0) {
    temp_names <- NULL
    exogen_names <- dimnames(object[["data"]][["original"]][["exogen"]])[[2]]
    temp_names <- paste0(exogen_names, ".l0")
    if (s > 0) {
      temp_names <- c(temp_names, paste0(exogen_names, ".l", rep(1:s, each = m))) 
    }
    if (add_block) {
      temp_names <- paste0("B\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  if (n > 0) {
    temp_names <- dimnames(object[["data"]][["original"]][["deterministic"]])[[2]]
    if (add_block) {
      temp_names <- paste0("C\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  if (object[["model"]][["structural"]]) {
    temp_names <- y_names
    if (add_block) {
      temp_names <- paste0("A0\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  return(x_names)
}


# Extracts the names of the regressors from a 'bvecmodel' object
# add_block adds the name of the block of cointegration, endogenous, exogenous,
# deterministic and structural coefficients
# The regressors of the error correction term are those of the cointegration
# matrix Pi and not those of the loading matrix alpha, since the draws of the
# former are what is reported for a VEC model

.get_regressor_names_bvecmodel <- function(object, add_block = FALSE) {

  k <- object[["model"]][["k"]]
  p <- object[["model"]][["p"]]
  m <- object[["model"]][["m"]]
  s <- object[["model"]][["s"]]
  n <- object[["model"]][["n"]]
  rank <- object[["model"]][["rank"]]

  x_names <- NULL

  add <- function(names, count, fallback, block) {
    if (count == 0) {
      return(NULL)
    }
    if (length(names) != count) {
      names <- fallback
    }
    if (add_block) {
      names <- paste0(block, "\n", names)
    }
    return(names)
  }

  if (rank > 0) {
    ect_names <- dimnames(object[["data"]][["train"]][["w"]])[[2]]
    k_beta <- k + m + object[["model"]][["n_restricted"]]
    fallback <- c(paste0("l.", object[["model"]][["endogen"]]),
                  if (m > 0) paste0("l.", object[["model"]][["exogen"]]),
                  if (object[["model"]][["n_restricted"]] > 0) paste0("l.d", 1:object[["model"]][["n_restricted"]]))
    x_names <- c(x_names, add(ect_names, k_beta, fallback, "Pi"))
  }

  reg_names <- dimnames(object[["data"]][["train"]][["x"]])[[2]]
  pos <- 0

  n_gamma <- k * (p - 1)
  if (n_gamma > 0) {
    fallback <- paste0("d.", rep(object[["model"]][["endogen"]], times = p - 1),
                       ".l", rep(.lag_label(1:(p - 1), p - 1), each = k))
    x_names <- c(x_names, add(reg_names[pos + 1:n_gamma], n_gamma, fallback, "Gamma"))
    pos <- pos + n_gamma
  }

  n_upsilon <- m * s
  if (n_upsilon > 0) {
    fallback <- paste0("d.", rep(object[["model"]][["exogen"]], times = s),
                       ".l", rep(.lag_label(0:(s - 1), s - 1), each = m))
    x_names <- c(x_names, add(reg_names[pos + 1:n_upsilon], n_upsilon, fallback, "Upsilon"))
    pos <- pos + n_upsilon
  }

  if (n > 0) {
    x_names <- c(x_names, add(reg_names[pos + 1:n], n, paste0("det.", 1:n), "C"))
  }

  if (object[["model"]][["structural"]]) {
    x_names <- c(x_names, add(object[["model"]][["endogen"]], k,
                              paste0("y", 1:k), "A0"))
  }

  return(x_names)
}
