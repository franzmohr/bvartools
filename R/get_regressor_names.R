
# Extracts the names of the regressors from a 'bvar' object
# add_block adds the letter of the block of endogenous, exogensous, deterministic, structural and sigma coefficients

.get_regressor_names_var <- function(object, add_block = FALSE) {
  
  k <- object[["specifications"]][["dims"]][["K"]]
  tvp <- object[["specifications"]][["tvp"]]
  y_names <- dimnames(object[["y"]])[[2]]
  x_names <- NULL
  
  m <- n <- o <- 0
  if (!is.null(object[["A"]])) {
    if (tvp[["A"]]) {
      m <- NCOL(object[["A"]][[1]]) / k
    } else {
      m <- NCOL(object[["A"]]) / k 
    }
    p <- m / k
    temp_names <- NULL
    if (!is.null(object[["x"]])) {
      temp_names <- c(temp_names, dimnames(object[["x"]])[[2]][1:m])
    } else {
      for (i in 1:p) {
        temp_names <- c(temp_names, paste(y_names, ".l", i, sep = ""))
      } 
    }
    if (add_block) {
      temp_names <- paste0("A\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  if (!is.null(object[["B"]])) {
    if (tvp[["B"]]) {
      n <- NCOL(object[["B"]][[1]]) / k
    } else {
      n <- NCOL(object[["B"]]) / k
    }
    temp_names <- NULL
    if (!is.null(object[["x"]])) {
      temp_names <- c(temp_names, dimnames(object[["x"]])[[2]][m + 1:n])
    } else {
      temp_names <- c(temp_names, paste("x", 1:n, sep = ""))
    }
    if (add_block) {
      temp_names <- paste0("B\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  if (!is.null(object[["C"]])) {
    if (tvp[["C"]]) {
      o <- NCOL(object[["C"]][[1]]) / k
    } else {
      o <- NCOL(object[["C"]]) / k 
    }
    temp_names <- NULL
    if (!is.null(object[["x"]])) {
      temp_names <- c(temp_names, dimnames(object[["x"]])[[2]][m + n + 1:o])
    } else {
      temp_names <- c(temp_names, paste("det", 1:o, sep = ""))
    }
    if (add_block) {
      temp_names <- paste0("C\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  if (!is.null(object[["A0"]])) {
    temp_names <- NULL
    temp_names <- dimnames(object[["y"]])[[2]]
    if (add_block) {
      temp_names <- paste0("A0\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  return(x_names)
}

.get_regressor_names_vec <- function(object, add_block = FALSE) {
  
  k <- object[["specifications"]][["dims"]][["K"]]
  tvp <- object[["specifications"]][["tvp"]]
  y_names <- dimnames(object[["y"]])[[2]]
  ect_names <- NULL
  x_names <- NULL
  
  n_pi <- n_pi_x <- n_pi_d <- m <- n <- o <- 0
  
  # Cointegration coefficients ----
  
  if (!is.null(object[["Pi"]])) {
    if (tvp[["Pi"]]) {
      n_pi <- NCOL(object[["Pi"]][[1]]) / k
    } else {
      n_pi <- NCOL(object[["Pi"]]) / k 
    }
    temp_names <- NULL
    if (!is.null(object[["w"]])) {
      temp_names <- c(temp_names, dimnames(object[["w"]])[[2]][1:n_pi])
    } else {
      temp_names <- c(temp_names, paste0("l.", gsub("d.", "", y_names)))
    }
    if (add_block) {
      temp_names <- paste0("Pi\n", temp_names)
    }
    ect_names <- c(ect_names, temp_names)
  }
  
  if (!is.null(object[["Pi_x"]])) {
    if (tvp[["Pi_x"]]) {
      n_pi_x <- NCOL(object[["Pi_x"]][[1]]) / k
    } else {
      n_pi_x <- NCOL(object[["Pi_x"]]) / k 
    }
    temp_names <- NULL
    if (!is.null(object[["w_x"]])) {
      temp_names <- c(temp_names, dimnames(object[["w_x"]])[[2]][1:n_pi_x])
    } else {
      temp_names <- c(temp_names, paste0("l.x", 1:n_pi_x))
    }
    if (add_block) {
      temp_names <- paste0("Pi_x\n", temp_names)
    }
    ect_names <- c(ect_names, temp_names)
  }
  
  if (!is.null(object[["Pi_d"]])) {
    if (tvp[["Pi_d"]]) {
      n_pi_d <- NCOL(object[["Pi_d"]][[1]]) / k
    } else {
      n_pi_d <- NCOL(object[["Pi_d"]]) / k 
    }
    temp_names <- NULL
    if (!is.null(object[["w_d"]])) {
      temp_names <- c(temp_names, dimnames(object[["w_d"]])[[2]][1:n_pi_d])
    } else {
      temp_names <- c(temp_names, paste0("l.d", 1:n_pi_d))
    }
    
    if (add_block) {
      temp_names <- paste0("Pi_d\n", temp_names)
    }
    ect_names <- c(ect_names, temp_names)
  }
  
  # Non-cointegration coefficients ----
  
  if (!is.null(object[["Gamma"]])) {
    if (tvp[["Gamma"]]) {
      m <- NCOL(object[["Gamma"]][[1]]) / k
    } else {
      m <- NCOL(object[["Gamma"]]) / k 
    }
    p <- m / k
    temp_names <- NULL
    if (!is.null(object[["x"]])) {
      temp_names <- c(temp_names, dimnames(object[["x"]])[[2]][1:m])
    } else {
      for (i in 1:p) {
        temp_names <- c(temp_names, paste(y_names, ".l", i, sep = ""))
      } 
    }
    if (add_block) {
      temp_names <- paste0("Gamma\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  if (!is.null(object[["Upsilon"]])) {
    if (tvp[["Upsilon"]]) {
      n <- NCOL(object[["Upsilon"]][[1]]) / k
    } else {
      n <- NCOL(object[["Upsilon"]]) / k
    }
    temp_names <- NULL
    if (!is.null(object[["x_x"]])) {
      temp_names <- c(temp_names, dimnames(object[["x_x"]])[[2]])
    } else {
      temp_names <- c(temp_names, paste0("x", 1:n))
    }
    if (add_block) {
      temp_names <- paste0("Upsilon\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  if (!is.null(object[["C"]])) {
    if (tvp[["C"]]) {
      o <- NCOL(object[["C"]][[1]]) / k
    } else {
      o <- NCOL(object[["C"]]) / k 
    }
    temp_names <- NULL
    if (!is.null(object[["x_d"]])) {
      temp_names <- c(temp_names, dimnames(object[["x_d"]])[[2]])
    } else {
      temp_names <- c(temp_names, paste("det.", 1:o, sep = ""))
    }
    if (add_block) {
      temp_names <- paste0("C\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  if (!is.null(object[["A0"]])) {
    temp_names <- NULL
    temp_names <- dimnames(object[["y"]])[[2]]
    if (add_block) {
      temp_names <- paste0("A0\n", temp_names)
    }
    x_names <- c(x_names, temp_names)
  }
  
  return(list("pi" = ect_names, "x" = x_names))
}