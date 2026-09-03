# bvar() and bvec() collect the output of a user-written sampler into the
# standardised object the rest of the package works with.

sampler_output <- function() {
  model <- fx_var_fitted()
  k <- model[["model"]][["k"]]
  coefficients <- t(model[["posterior"]][["a"]][["coeffs"]])
  precision <- t(model[["posterior"]][["u_sigma_inv"]][["coeffs"]])

  list(
    y = model[["data"]][["train"]][["y"]],
    x = model[["data"]][["train"]][["x"]],
    A = coefficients[seq_len(k * k), ],
    C = coefficients[k * k + seq_len(k), ],
    Sigma = apply(precision, 2, function(z) c(solve(matrix(z, k))))
  )
}

test_that("bvar assembles a model object from posterior draws", {
  draws <- sampler_output()
  object <- bvar(y = draws[["y"]], x = draws[["x"]], A = draws[["A"]],
                 C = draws[["C"]], Sigma = draws[["Sigma"]])

  expect_s3_class(object, "bvarmodel")
  expect_named(object, c("model", "data", "posterior"))
  expect_identical(object[["model"]][["type"]], "VAR")
  expect_identical(object[["model"]][["k"]], 3L)
  expect_identical(object[["model"]][["p"]], 1L)
  expect_identical(object[["model"]][["n"]], 1L)
  expect_identical(object[["model"]][["endogen"]], c("invest", "income", "cons"))
})

test_that("bvar puts the draws back in the internal layout", {
  draws <- sampler_output()
  object <- bvar(y = draws[["y"]], x = draws[["x"]], A = draws[["A"]],
                 C = draws[["C"]], Sigma = draws[["Sigma"]])
  coefficients <- object[["posterior"]][["a"]][["coeffs"]]

  expect_s3_class(coefficients, "mcmc")
  expect_identical(dim(coefficients), c(fx_iterations, 12L))
  # Lag coefficients come first, deterministic ones after.
  expect_equal(unname(as.matrix(coefficients[, 1:9])),
               unname(t(draws[["A"]])))
  expect_equal(unname(as.matrix(coefficients[, 10:12])),
               unname(t(draws[["C"]])))
  # Sigma is stored as its inverse.
  expect_equal(
    matrix(object[["posterior"]][["u_sigma_inv"]][["coeffs"]][1, ], 3),
    solve(matrix(draws[["Sigma"]][, 1], 3))
  )
})

test_that("an object built by bvar is usable downstream", {
  draws <- sampler_output()
  object <- bvar(y = draws[["y"]], x = draws[["x"]], A = draws[["A"]],
                 C = draws[["C"]], Sigma = draws[["Sigma"]])

  expect_s3_class(irf(object, impulse = "income", response = "cons",
                      n_ahead = 3), "bvarirf")
  expect_s3_class(fevd(object, response = "cons", n_ahead = 3), "bvarfevd")
  expect_no_error(summary(object))
})

test_that("bvar requires time series input", {
  draws <- sampler_output()

  expect_error(bvar(y = unclass(draws[["y"]]), x = draws[["x"]],
                    A = draws[["A"]], Sigma = draws[["Sigma"]]),
               "class time-series")
  expect_error(bvar(y = draws[["y"]], x = unclass(draws[["x"]]),
                    A = draws[["A"]], Sigma = draws[["Sigma"]]),
               "class time-series")
})

test_that("bvec assembles a VEC model object from posterior draws", {
  model <- fx_vec_fitted()
  k <- model[["model"]][["k"]]
  r <- model[["model"]][["rank"]]
  coefficients <- t(model[["posterior"]][["a"]][["coeffs"]])
  precision <- t(model[["posterior"]][["u_sigma_inv"]][["coeffs"]])

  object <- bvec(y = model[["data"]][["train"]][["y"]],
                 w = model[["data"]][["train"]][["w"]],
                 x = model[["data"]][["train"]][["x"]],
                 r = r,
                 alpha = coefficients[seq_len(k * r), ],
                 beta = t(model[["posterior"]][["beta"]][["coeffs"]]),
                 C = coefficients[k * r + seq_len(k), ],
                 Sigma = apply(precision, 2,
                               function(z) c(solve(matrix(z, k)))))

  expect_s3_class(object, "bvecmodel")
  expect_identical(object[["model"]][["type"]], "VEC")
  expect_identical(object[["model"]][["k"]], k)
  expect_equal(object[["model"]][["rank"]], r)
  expect_identical(nrow(object[["posterior"]][["beta"]][["coeffs"]]),
                   fx_iterations)
})

test_that("bvec requires time series input", {
  model <- fx_vec_fitted()

  expect_error(bvec(y = unclass(model[["data"]][["train"]][["y"]])),
               "must be of class 'ts'")
})
