# The exported building blocks of the cointegration sampler of Koop,
# Leon-Gonzalez and Strachan (2010). VecKlgs2010 does not call them -- it runs
# the vendored core -- so these are the only place they are exercised.

coint_inputs <- function(p) {
  data("e6", package = "bvartools", envir = environment())
  temp <- create_bvecmodel(e6, p = p, r = 1)
  y <- t(temp$data$train$y)
  list(y = y,
       w = t(temp$data$train$w),
       x = if (is.null(temp$data$train$x)) NULL else t(temp$data$train$x),
       sigma_i = solve(tcrossprod(y) / ncol(y)))
}


test_that("the second reparameterisation keeps Pi and makes beta semiorthogonal", {

  alpha <- matrix(c(-0.07, 0.17), 2)
  beta <- matrix(c(1, -4), 2)

  out <- coint_kls2010_reparameterise_two(alpha, beta)

  expect_named(out, c("alpha", "beta"))
  expect_equal(out$alpha %*% t(out$beta), alpha %*% t(beta))
  expect_equal(crossprod(out$beta), diag(1, 1))

  # Two cointegration vectors: the product survives, and the columns of beta
  # come back orthonormal.
  alpha2 <- matrix(c(0.1, -0.2, 0.3, 0.05, 0.4, -0.1), 3)
  beta2 <- matrix(c(1, 0, -2, 0, 1, 3), 3)
  out2 <- coint_kls2010_reparameterise_two(alpha2, beta2)
  expect_equal(out2$alpha %*% t(out2$beta), alpha2 %*% t(beta2))
  expect_equal(crossprod(out2$beta), diag(1, 2))
})


test_that("post_coint_kls draws alpha, beta, Pi and Gamma of matching shapes", {

  d <- coint_inputs(p = 2)
  k <- nrow(d$y)
  n_beta <- nrow(d$w)
  n_x <- nrow(d$x)
  n_ag <- k * (1 + n_x)

  draw <- function() {
    post_coint_kls(y = d$y, beta = matrix(c(1, -4), n_beta), w = d$w, x = d$x,
                   sigma_i = d$sigma_i, v_i = 0, p_tau_i = diag(1, n_beta),
                   g_i = d$sigma_i,
                   gamma_mu_prior = matrix(0, n_ag),
                   gamma_v_i_prior = diag(0, n_ag))
  }

  set.seed(11)
  out <- draw()
  expect_true(all(c("alpha", "beta", "Pi", "Gamma") %in% names(out)))
  expect_equal(dim(out$alpha), c(k, 1L))
  expect_equal(dim(out$beta), c(n_beta, 1L))
  expect_equal(dim(out$Pi), c(k, n_beta))
  # vec(Gamma), K N x 1, which is how bvec() and the vignette consume it
  expect_equal(dim(out$Gamma), c(k * n_x, 1L))
  expect_equal(dim(matrix(out$Gamma, k)), c(k, n_x))
  expect_equal(out$Pi, out$alpha %*% t(out$beta))
  expect_true(all(is.finite(unlist(out))))

  # The draw comes from R's generator, so the seed repeats it.
  set.seed(11)
  expect_identical(draw(), out)
})


test_that("post_coint_kls_sur draws alpha, beta and Pi of matching shapes", {

  d <- coint_inputs(p = 1)
  k <- nrow(d$y)
  n_beta <- nrow(d$w)

  set.seed(12)
  out <- post_coint_kls_sur(y = d$y, beta = matrix(c(1, -4), k), w = d$w,
                            sigma_i = d$sigma_i, v_i = 0,
                            p_tau_i = diag(1, n_beta), g_i = d$sigma_i)

  expect_equal(dim(out$alpha), c(k, 1L))
  expect_equal(dim(out$beta), c(n_beta, 1L))
  expect_equal(dim(out$Pi), c(k, n_beta))
  expect_equal(out$Pi, out$alpha %*% t(out$beta))
  expect_true(all(is.finite(out$Pi)))
})
