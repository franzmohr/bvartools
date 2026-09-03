# A small estimated VAR to run the measures on.
make_bvar <- function(p = 2, iterations = 60, burnin = 10, seed = 1) {
  set.seed(seed)
  data("e1", envir = environment())
  y <- diff(log(get("e1", envir = environment()))) * 100
  y <- stats::window(y, end = c(1978, 4))

  model <- create_bvarmodel(y, p = p, deterministic = "const",
                            iterations = iterations, burnin = burnin)
  model <- add_priors(model,
                      coef = list(v_i = 1, v_i_det = 1 / 10),
                      sigma = list(df = "k", scale = 1))
  model <- add_initial_values(model)
  add_posterior_coefficients(model)
}

# One draw in the shape .collect_draws produces, so that the worker can be
# driven directly with coefficients chosen by hand.
one_draw <- function(A, Sigma) {
  list(list(A = A, Sigma = Sigma))
}

test_that("every row of the table sums to one", {
  object <- make_bvar()
  A <- bvartools:::.collect_draws(object)

  for (type in c("gir", "oir")) {
    tables <- lapply(A, bvartools:::.spillover_table, h = 10, type = type)
    row_totals <- unlist(lapply(tables, rowSums))
    expect_equal(row_totals, rep(1, length(row_totals)), info = type)
  }
})

test_that("the orthogonalised shares add up before normalisation", {
  # Under a Choleski factorisation the decomposition is exhaustive by
  # construction, so the normalisation step must be a no-op there.
  set.seed(2)
  k <- 3
  A <- matrix(stats::rnorm(k * k, sd = 0.2), k)
  Sigma <- crossprod(matrix(stats::rnorm(k * k), k)) + diag(k)

  table <- bvartools:::.spillover_table(one_draw(A, Sigma)[[1]], h = 8,
                                        type = "oir")
  expect_equal(rowSums(table), rep(1, k))
})

test_that("from, to and the total obey the accounting identity", {
  object <- make_bvar()
  sp <- spillover(object, n_ahead = 10, keep_draws = TRUE)

  from <- as.matrix(sp$from)
  to <- as.matrix(sp$to)
  total <- as.numeric(sp$total)

  expect_equal(rowSums(from), total)
  expect_equal(rowSums(to), total)
  expect_equal(as.matrix(sp$net), to - from, ignore_attr = TRUE)
})

test_that("a system with no transmission has a zero index", {
  # Diagonal coefficients and a diagonal covariance: nothing moves between
  # variables, so the table is the identity and the index is zero.
  k <- 3
  A <- cbind(diag(0.5, k), diag(0.2, k))
  Sigma <- diag(c(1, 4, 0.25))

  table <- bvartools:::.spillover_table(one_draw(A, Sigma)[[1]], h = 10,
                                        type = "gir")
  expect_equal(table, diag(k))

  # And the same through the whole path, on a model whose draws are replaced by
  # that system.
  object <- make_bvar(p = 1)
  store <- nrow(object$posterior$u_sigma_inv$coeffs)
  kk <- object$model$k^2
  a_diag <- as.numeric(diag(0.5, object$model$k))
  object$posterior$a$coeffs[, 1:kk] <- matrix(a_diag, store, kk, byrow = TRUE)
  object$posterior$u_sigma_inv$coeffs <-
    matrix(as.numeric(diag(1, object$model$k)), store, kk, byrow = TRUE)

  sp <- spillover(object, n_ahead = 5)
  expect_equal(unname(sp$total["50%", 1]), 0)
})

test_that("a two variable system matches the index computed by hand", {
  k <- 2
  A <- matrix(c(0.5, 0.3, 0.2, 0.4), k)  # column major: a11 a21 a12 a22
  Sigma <- matrix(c(1, 0.3, 0.3, 2), k)
  h <- 4

  table <- bvartools:::.spillover_table(one_draw(A, Sigma)[[1]], h = h,
                                        type = "gir")

  # The Pesaran and Shin decomposition, written out independently of the worker.
  phi <- list(diag(k))
  for (i in 1:h) {
    phi[[i + 1]] <- phi[[i]] %*% A
  }
  num <- matrix(0, k, k)
  den <- numeric(k)
  for (i in seq_along(phi)) {
    num <- num + (phi[[i]] %*% Sigma)^2
    den <- den + diag(phi[[i]] %*% Sigma %*% t(phi[[i]]))
  }
  theta <- num / den                      # row j by its own mse
  theta <- t(t(theta) / diag(Sigma))      # column i by the shock variance
  theta <- theta / rowSums(theta)

  expect_equal(table, theta)

  expected_total <- 100 * (sum(theta) - sum(diag(theta))) / k
  expect_equal(100 * (sum(table) - sum(diag(table))) / k, expected_total)
})

test_that("the scaling is Pesaran-Shin and not the one .vardecomp uses", {
  # With unit variances the impulse scaling is one, so the new worker and the
  # row-normalised .vardecomp must agree exactly.
  k <- 3
  set.seed(4)
  A <- matrix(stats::rnorm(k * k, sd = 0.2), k)
  corr <- diag(k)
  corr[1, 2] <- corr[2, 1] <- 0.4
  corr[1, 3] <- corr[3, 1] <- 0.2

  table <- bvartools:::.spillover_table(one_draw(A, corr)[[1]], h = 6,
                                        type = "gir")
  for (j in 1:k) {
    vd <- bvartools:::.vardecomp(list(A = A, Sigma = corr), h = 6,
                                 type = "gir", response = j)
    row <- vd[nrow(vd), ]
    expect_equal(as.numeric(table[j, ]), as.numeric(row / sum(row)))
  }

  # With unequal variances they must differ, by exactly the impulse variance.
  Sigma <- diag(sqrt(c(1, 9, 0.25))) %*% corr %*% diag(sqrt(c(1, 9, 0.25)))
  table <- bvartools:::.spillover_table(one_draw(A, Sigma)[[1]], h = 6,
                                        type = "gir")
  for (j in 1:k) {
    vd <- bvartools:::.vardecomp(list(A = A, Sigma = Sigma), h = 6,
                                 type = "gir", response = j)
    row <- vd[nrow(vd), ] / diag(Sigma)
    expect_equal(as.numeric(table[j, ]), as.numeric(row / sum(row)))
  }
})

test_that("the index is computed per draw and not from the mean table", {
  object <- make_bvar()
  draws <- spillover(object, keep_draws = TRUE)
  summarised <- spillover(object)

  expect_equal(unname(summarised$total["50%", 1]),
               unname(stats::median(as.numeric(draws$total))))
  # A ratio of sums, so the spread across draws is not degenerate.
  expect_gt(stats::sd(as.numeric(draws$total)), 0)
})

test_that("the returned object has the promised shape", {
  object <- make_bvar()
  sp <- spillover(object, n_ahead = 8, ci = 0.9)
  k <- object$model$k
  varnames <- object$model$endogen

  expect_s3_class(sp, "bvarspillover")
  expect_equal(dim(sp$total), c(3, 1))
  expect_equal(dim(sp$from), c(3, k))
  expect_equal(colnames(sp$from), varnames)
  expect_equal(rownames(sp$from), c("5%", "50%", "95%"))
  expect_equal(dim(sp$table), c(k + 1, k + 1))
  expect_equal(dim(sp$pairwise), c(k, k))
  expect_equal(sp$specification$n_ahead, 8)
  expect_equal(sp$specification$k, k)

  # Net pairwise spillovers are antisymmetric by construction.
  expect_equal(sp$pairwise, -t(sp$pairwise))
})

test_that("the index sits in a plausible range on real data", {
  object <- make_bvar()
  sp <- spillover(object)
  total <- unname(sp$total["50%", 1])

  # Pinned at either end means a margin or the normalisation is wrong, which no
  # identity check would catch.
  expect_gt(total, 1)
  expect_lt(total, 99)
})

test_that("print and plot work", {
  object <- make_bvar()
  sp <- spillover(object)

  expect_output(print(sp), "Spillover index")
  expect_output(print(sp), "total index")

  png(tempfile())
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(sp))
  expect_invisible(plot(sp, which = "to"))
  expect_invisible(plot(sp, which = "from"))
})

test_that("impossible arguments are refused", {
  object <- make_bvar()

  expect_error(spillover(object, type = "sir"), "'gir' or 'oir'")
  expect_error(spillover(object, n_ahead = 0), "at least 1")
  expect_error(spillover(object, ci = 1), "between 0 and 1")

  no_sigma <- object
  no_sigma$posterior$u_sigma_inv <- NULL
  expect_error(spillover(no_sigma), "variance-covariance")
})

test_that("fevd still agrees with itself after the refactor", {
  # .collect_draws was lifted out of fevd.bvarmodel, so the decomposition it
  # produces has to be unchanged.
  object <- make_bvar()
  vd <- fevd(object, response = object$model$endogen[1], n_ahead = 5)

  expect_s3_class(vd, "bvarfevd")
  expect_equal(dim(vd), c(6, object$model$k))
  expect_equal(as.numeric(rowSums(vd)), rep(1, 6))
})

test_that("the sample length is read from the regressors, not from y's shape", {
  # train$y is one row per period and one column per variable, so dividing its
  # row count by k -- which is what this code did while gen_var stacked the
  # series -- understates the sample by a factor of k. The regressor matrix is
  # unambiguous, so the count comes from there.
  object <- make_bvar(p = 2)
  k <- object$model$k

  expect_equal(bvartools:::.train_periods(object, k),
               nrow(object$data$train$z) / k)
  expect_equal(bvartools:::.train_periods(object, k),
               nrow(object$data$train$y))

  # A stacked single column still reads as the same sample.
  stacked <- object
  stacked$data$train$z <- NULL
  stacked$data$train$y <- matrix(as.numeric(t(object$data$train$y)), ncol = 1)
  expect_equal(bvartools:::.train_periods(stacked, k),
               nrow(object$data$train$y))
})

test_that("period selects a block of a time varying model", {
  set.seed(1)
  data("e1", envir = environment())
  y <- diff(log(get("e1", envir = environment()))) * 100
  y <- stats::window(y, end = c(1978, 4))

  model <- create_bvarmodel(y, p = 1, deterministic = "const", tvp = TRUE,
                            iterations = 60, burnin = 10)
  model <- add_priors(model,
                      coef = list(v_i = 1, v_i_det = 1 / 10,
                                  shape = 3, rate = 0.0001),
                      sigma = list(df = "k", scale = 1))
  object <- add_posterior_coefficients(add_initial_values(model))

  # One coefficient vector per period is stored end to end, so the width of the
  # posterior divided by the width of one period is the sample length.
  tt <- bvartools:::.train_periods(object, object$model$k)
  expect_equal(tt, ncol(object$posterior$a$coeffs) / ncol(object$data$train$z))

  # The last period is the default, and every period in the sample is reachable.
  expect_equal(spillover(object, n_ahead = 5, period = tt)$total,
               spillover(object, n_ahead = 5)$total)
  expect_silent(spillover(object, n_ahead = 5, period = 1))
  expect_error(spillover(object, n_ahead = 5, period = tt + 1), "Implausible")

  # Different periods give different coefficients and so different measures.
  expect_false(identical(spillover(object, n_ahead = 5, period = 2)$total,
                         spillover(object, n_ahead = 5, period = tt)$total))
})

test_that("the printed table is in one set of units and adds up", {
  object <- make_bvar()
  sp <- spillover(object)
  k <- object$model$k
  tab <- sp$table

  # Body, margins and corner are all percentages: a body row plus nothing else
  # is 100, and each margin sums to the corner.
  expect_equal(unname(rowSums(tab[seq_len(k), seq_len(k)])), rep(100, k))
  expect_equal(sum(tab[seq_len(k), k + 1]), tab[k + 1, k + 1])
  expect_equal(sum(tab[k + 1, seq_len(k)]), tab[k + 1, k + 1])
  expect_false(anyNA(tab))

  # And the corner is the mean of the index, so it sits inside its interval.
  expect_gt(tab[k + 1, k + 1], sp$total["2.5%", 1])
  expect_lt(tab[k + 1, k + 1], sp$total["97.5%", 1])
})
