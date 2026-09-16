# Structural models with four endogenous variables.
#
# The free elements of A0 are stored column by column: (2,1), (3,1), (4,1),
# (3,2), (4,2), (4,3). With three variables that is also the row by row order,
# which is how three places read them the wrong way without a test noticing.

us_four_variables <- function() {
  us <- bvartools::us_macrodata
  data <- stats::ts.intersect(us, dr = diff(us[, "r"]))
  colnames(data) <- c("Dp", "u", "r", "dr")
  data
}

structural_var_k4 <- function(varsel = "none") {
  model <- create_bvarmodel(us_four_variables(), p = 1, deterministic = "const",
                            structural = TRUE, error = "gamma", varsel = varsel,
                            iterations = 40, burnin = 10)
  model <- suppressWarnings(
    add_priors(model, coef = list(v_i = 1, v_i_det = 0.1), sigma = list(shape = 3, rate = 0.01),
               varsel = if (varsel == "bvs") list(inprior = 0.5) else NULL))
  set.seed(20260914)
  model <- add_initial_values(model)
  add_posterior_coefficients(model)
}

# A0 of draw i from the layout of the draws.
a0_of_draw <- function(model, i) {
  k <- model[["model"]][["k"]]
  a <- as.matrix(model[["posterior"]][["a"]][["coeffs"]])
  A0 <- diag(k)
  A0[lower.tri(A0)] <- a[i, ncol(a) - k * (k - 1) / 2 + seq_len(k * (k - 1) / 2)]
  A0
}

test_that("the contemporaneous regressors are the free elements of A0 column by column", {
  model <- structural_var_k4()
  y <- as.matrix(model[["data"]][["train"]][["y"]])
  z <- model[["data"]][["train"]][["z"]]

  expect_equal(z[, ncol(z) - 5:0], kronecker(-y, diag(4))[, which(lower.tri(diag(4)))])
})

test_that("structural responses and decompositions read A0 in the layout of the draws", {
  model <- structural_var_k4()
  k <- 4
  a <- as.matrix(model[["posterior"]][["a"]][["coeffs"]])
  draws <- bvartools:::.collect_draws(model, need_A0 = TRUE)

  for (i in c(1, nrow(a))) {
    expect_equal(draws[[i]][["A0"]], a0_of_draw(model, i))
  }

  # The response of dr to a unit structural shock to Dp, written out from that A0.
  i <- 2
  a0_inv <- solve(a0_of_draw(model, i))
  B <- a0_inv %*% matrix(a[i, 1:(k * k)], k)
  expected <- numeric(4)
  phi <- diag(k)
  for (h in 0:3) {
    expected[h + 1] <- (phi %*% a0_inv)[4, 1]
    phi <- phi %*% B
  }
  response <- irf(model, impulse = "Dp", response = "dr", n_ahead = 3, type = "sir",
                  keep_draws = TRUE)
  expect_equal(as.numeric(response[i, ]), expected)

  # And the impact decomposition of the variance of dr, averaged over the draws.
  omega <- function(i) diag(1 / diag(matrix(model[["posterior"]][["u_sigma_inv"]][["coeffs"]][i, ], k)))
  shares <- Reduce(`+`, lapply(seq_len(nrow(a)), function(i) {
    impact <- solve(a0_of_draw(model, i)) %*% sqrt(omega(i))
    impact[4, ]^2 / sum(impact[4, ]^2)
  })) / nrow(a)
  decomposition <- fevd(model, response = "dr", n_ahead = 0, type = "sir")
  expect_equal(as.numeric(decomposition[1, ]), shares)
})

test_that("summary puts each element of A0 and its inclusion probability in its own cell", {
  model <- structural_var_k4()
  a <- as.matrix(model[["posterior"]][["a"]][["coeffs"]])
  a0_mean <- Reduce(`+`, lapply(seq_len(nrow(a)), function(i) a0_of_draw(model, i))) / nrow(a)

  means <- summary(model)[["a"]][["means"]]
  expect_equal(unname(means[, ncol(means) - 3:0]), a0_mean)

  selected <- structural_var_k4(varsel = "bvs")
  lambda <- as.matrix(selected[["posterior"]][["a"]][["lambda"]])
  inclusion <- diag(4)
  inclusion[lower.tri(inclusion)] <- colMeans(lambda)[ncol(lambda) - 5:0]
  reported <- summary(selected)[["a"]][["lambda"]]
  expect_equal(unname(reported[, ncol(reported) - 3:0]), inclusion)
})

test_that("a structural VEC model can be summarised, with A0 column by column", {
  set.seed(1)
  data <- us_four_variables()[, c("Dp", "u", "r")]
  # A fourth series of its own, since a transformation of the other three puts
  # collinear columns into the error correction term.
  data <- cbind(data, noise = stats::ts(cumsum(stats::rnorm(nrow(data))),
                                        start = stats::start(data), frequency = 4))
  colnames(data) <- c("Dp", "u", "r", "noise")

  model <- create_bvecmodel(data, p = 2, r = 1, const = "unrestricted", structural = TRUE,
                            error = "gamma", iterations = 40, burnin = 10)
  model <- suppressWarnings(add_priors(model, coef = list(v_i = 1, v_i_det = 0.1),
                                       coint = list(v_i = 0, p_tau_i = 1),
                                       sigma = list(shape = 3, rate = 0.01)))
  model <- add_initial_values(model)
  model <- add_posterior_coefficients(model)

  a <- as.matrix(model[["posterior"]][["a"]][["coeffs"]])
  a0_mean <- Reduce(`+`, lapply(seq_len(nrow(a)), function(i) a0_of_draw(model, i))) / nrow(a)

  s <- summary(model)
  means <- s[["a"]][["means"]]
  expect_identical(colnames(means)[ncol(means) - 3:0], c("Dp", "u", "r", "noise"))
  expect_equal(unname(means[, ncol(means) - 3:0]), a0_mean)
  expect_output(print(s))
})

test_that("a structural model with four variables can be plotted", {
  model <- structural_var_k4()
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(model))
})
