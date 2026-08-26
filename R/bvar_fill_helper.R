# Helper functions used by bvar() and bvec() to turn user provided posterior
# draws into the layout of a 'bvarmodel' or 'bvecmodel' object.


# Ensures that a time-series object is a matrix with column names, which is the
# form in which create_bvarmodel() and create_bvecmodel() store their data
.name_ts <- function(x, prefix) {

  if (is.null(x)) {
    return(NULL)
  }

  if (!is.null(dimnames(x)[[2]])) {
    return(x)
  }

  ts_info <- stats::tsp(x)
  result <- stats::ts(as.matrix(x), class = c("mts", "ts", "matrix"))
  stats::tsp(result) <- ts_info
  if (NCOL(result) == 1) {
    dimnames(result) <- list(NULL, prefix)
  } else {
    dimnames(result) <- list(NULL, paste0(prefix, 1:NCOL(result)))
  }

  return(result)
}


# Turns a subset of the columns of a data matrix into a time-series object with
# the time-series attributes of that matrix
.subset_ts <- function(x, pos, names = NULL) {

  ts_info <- stats::tsp(x)
  if (is.null(names)) {
    names <- dimnames(x)[[2]][pos]
  }

  result <- stats::ts(as.matrix(x[, pos, drop = FALSE]), class = c("mts", "ts", "matrix"))
  stats::tsp(result) <- ts_info
  dimnames(result) <- list(NULL, names)

  return(result)
}


# Splits a coefficient argument of bvar() or bvec() into the matrices of
# coefficient, inclusion and state variance draws
.split_draws <- function(draws, name) {

  if (is.null(draws)) {
    return(NULL)
  }

  if (is.list(draws)) {
    if (!"coeffs" %in% names(draws)) {
      stop("If argument '", name, "' is a list, it must contain an element 'coeffs'.")
    }
    result <- list("coeffs" = as.matrix(draws[["coeffs"]]),
                   "lambda" = if (is.null(draws[["lambda"]])) NULL else as.matrix(draws[["lambda"]]),
                   "sigma" = if (is.null(draws[["sigma"]])) NULL else as.matrix(draws[["sigma"]]))
  } else {
    result <- list("coeffs" = as.matrix(draws), "lambda" = NULL, "sigma" = NULL)
  }

  return(result)
}


# Determines the number of coefficients per period of a block of draws and
# whether those draws are time varying. Argument 'unit' is the number of rows a
# single period of the block is a multiple of, i.e. K for coefficients of one
# equation block and K^2 for square matrices.
.block_dims <- function(block, tt, unit, name) {

  if (is.null(block)) {
    return(NULL)
  }

  n <- nrow(block[["coeffs"]])

  # Draws are considered time varying if they contain a plausible number of
  # coefficients for every period of the sample
  tvp <- n %% tt == 0 && n / tt >= unit && (n / tt) %% unit == 0
  if (tvp) {
    n <- n / tt
  }

  if (n %% unit != 0) {
    stop("Row number of the draws of '", name, "' is not a multiple of ", unit, ".")
  }

  block[["n"]] <- n
  block[["tvp"]] <- tvp

  return(block)
}


# Reduces draws of the K x K matrix of structural coefficients to the free
# elements below its diagonal, which is how they enter the data matrix in SUR
# form and, hence, the vector of estimated coefficients
.a0_free_elements <- function(block, k) {

  pos <- which(lower.tri(matrix(NA_real_, k, k)))
  n_periods <- nrow(block[["coeffs"]]) / (k * k)
  pos <- as.vector(outer(pos, (1:n_periods - 1) * k * k, "+"))

  for (i in c("coeffs", "lambda", "sigma")) {
    if (!is.null(block[[i]])) {
      block[[i]] <- block[[i]][pos, , drop = FALSE]
    }
  }

  return(block)
}


# Binds the blocks of coefficient draws of a model in the order in which their
# regressors appear in the data matrix and transposes them into the S x N
# layout of an 'mcmc' object. For models with time varying parameters the
# coefficients of every period are appended to each other, where blocks of
# constant coefficients repeat their draws. Argument 'fill' is the value used
# for blocks that do not provide the requested element.
.bind_blocks <- function(blocks, tt, tvp, element = "coeffs", fill = NA_real_) {

  blocks <- blocks[!vapply(blocks, is.null, logical(1))]
  if (length(blocks) == 0) {
    return(NULL)
  }

  if (element != "coeffs") {
    if (!any(vapply(blocks, function(x) {!is.null(x[[element]])}, logical(1)))) {
      return(NULL)
    }
  }

  draws <- ncol(blocks[[1]][["coeffs"]])
  n_period <- sum(vapply(blocks, function(x) {x[["n"]]}, numeric(1)))

  # Inclusion parameters and state variances are not time varying, since they
  # describe the state equation of a coefficient and not the coefficient itself
  if (element != "coeffs" | !tvp) {

    result <- matrix(NA_real_, draws, n_period)
    pos <- 0
    for (i in blocks) {
      temp <- i[[element]]
      if (is.null(temp)) {
        temp <- matrix(fill, i[["n"]], draws)
      }
      result[, pos + 1:i[["n"]]] <- t(temp)
      pos <- pos + i[["n"]]
    }

    return(result)
  }

  result <- matrix(NA_real_, draws, tt * n_period)
  pos <- 0
  for (i in blocks) {
    for (j in 1:tt) {
      rows <- if (i[["tvp"]]) (j - 1) * i[["n"]] + 1:i[["n"]] else 1:i[["n"]]
      result[, (j - 1) * n_period + pos + 1:i[["n"]]] <- t(i[["coeffs"]][rows, , drop = FALSE])
    }
    pos <- pos + i[["n"]]
  }

  return(result)
}


# Inverts draws of the covariance matrix of the error term, since 'bvarmodel'
# and 'bvecmodel' objects contain the draws of its inverse
.invert_sigma_draws <- function(sigma, k) {

  draws <- ncol(sigma)
  n_periods <- nrow(sigma) / (k * k)
  result <- matrix(NA_real_, draws, nrow(sigma))

  for (i in 1:n_periods) {
    pos <- (i - 1) * k * k + 1:(k * k)
    temp <- sigma[pos, , drop = FALSE]
    temp <- vapply(1:draws,
                   function(j) {as.numeric(solve(matrix(temp[, j], k)))},
                   numeric(k * k))
    result[, pos] <- t(temp)
  }

  return(result)
}


# Binds time-series objects into a single one, dropping those that are not
# specified and keeping the column names of the remaining ones
.cbind_ts <- function(objects) {

  objects <- objects[!vapply(objects, is.null, logical(1))]
  if (length(objects) == 0) {
    return(NULL)
  }

  ts_info <- stats::tsp(objects[[1]])
  names <- unlist(lapply(objects, function(x) {dimnames(x)[[2]]}))

  result <- do.call(cbind, lapply(objects, as.matrix))
  result <- stats::ts(result, class = c("mts", "ts", "matrix"))
  stats::tsp(result) <- ts_info
  dimnames(result) <- list(NULL, names)

  return(result)
}


# Binds the draws of the cointegration matrix of the endogenous, unmodelled and
# restricted deterministic regressors into the draws of the matrix of all
# cointegration relations, which is stored as its vectorised form
.bind_beta <- function(blocks, sizes, rank, k_beta) {

  names <- c("beta", "beta_x", "beta_d")

  pos <- which(!vapply(blocks, is.null, logical(1)))
  if (length(pos) == 0) {
    stop("No draws of the cointegration matrix provided.")
  }
  draws <- ncol(blocks[[pos[1]]][["coeffs"]])

  result <- matrix(NA_real_, k_beta * rank, draws)
  for (j in 1:rank) {
    pos <- (j - 1) * k_beta
    for (i in seq_along(blocks)) {
      if (sizes[i] == 0) {
        next
      }
      if (is.null(blocks[[i]])) {
        stop("Argument '", names[i], "' must be specified for models with a cointegration rank larger than zero.")
      }
      result[pos + 1:sizes[i], ] <- blocks[[i]][["coeffs"]][(j - 1) * sizes[i] + 1:sizes[i], , drop = FALSE]
      pos <- pos + sizes[i]
    }
  }

  return(result)
}


# Checks that all blocks of draws of a model contain the same number of draws
# and returns that number
.count_draws <- function(blocks) {

  blocks <- blocks[!vapply(blocks, is.null, logical(1))]
  if (length(blocks) == 0) {
    stop("No posterior draws provided.")
  }

  draws <- unique(vapply(blocks, function(x) {ncol(x[["coeffs"]])}, numeric(1)))
  if (length(draws) > 1) {
    stop("The provided coefficient draws do not correspond to the same number of MCMC iterations.")
  }

  return(as.integer(draws))
}


# Validates the specification of the covariance matrix of the error term of a
# model that is constructed from posterior draws. Since those draws contain the
# full covariance matrix, they only reveal whether it is time varying, so the
# more specific specifications have to be provided by the user.
.check_error_spec <- function(error, sv, k) {

  if (is.null(error)) {
    if (sv) {
      return(ifelse(k == 1, "sv", "sv+covar"))
    }
    return("wishart")
  }

  if (!error %in% c("wishart", "gamma", "gamma+covar", "sv", "sv+covar")) {
    stop("Invalid specification of argument 'error'.")
  }
  if (sv & !error %in% c("sv", "sv+covar")) {
    stop("Argument 'error' must be 'sv' or 'sv+covar' if the draws of the covariance matrix are time varying.")
  }
  if (!sv & error %in% c("sv", "sv+covar")) {
    stop("Argument 'error' can only be 'sv' or 'sv+covar' if the draws of the covariance matrix are time varying.")
  }

  return(error)
}


# Name of the posterior simulation algorithm of a VEC model with the given
# specification of the error term, as determined by create_bvecmodel()
.vec_algorithm <- function(error, tvp) {
  return(sub("^Var", "Vec", .var_algorithm(error, tvp)))
}
