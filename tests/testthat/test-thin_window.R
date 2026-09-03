test_that("thin keeps every nth draw", {
  model <- fx_var_fitted()
  thinned <- thin(model, thin = 3)
  draws <- model[["posterior"]][["a"]][["coeffs"]]

  expect_s3_class(thinned, "bvarmodel")
  expect_identical(nrow(thinned[["posterior"]][["a"]][["coeffs"]]),
                   as.integer(fx_iterations / 3))
  expect_equal(unname(as.matrix(thinned[["posterior"]][["a"]][["coeffs"]])),
               unname(as.matrix(draws[seq(3, fx_iterations, 3), ])))
})

test_that("thin applies to every posterior element", {
  thinned <- thin(fx_var_fitted(), thin = 2)
  expected <- as.integer(fx_iterations / 2)

  expect_identical(nrow(thinned[["posterior"]][["a"]][["coeffs"]]), expected)
  expect_identical(nrow(thinned[["posterior"]][["u_sigma_inv"]][["coeffs"]]),
                   expected)
  expect_identical(nrow(thinned[["posterior"]][["loglik"]]), expected)
})

test_that("thinning by one leaves the draws untouched", {
  expect_equal(thin(fx_var_fitted(), thin = 1)[["posterior"]],
               fx_var_fitted()[["posterior"]])
})

test_that("thin works for VEC models and model lists", {
  vec <- thin(fx_vec_fitted(), thin = 2)
  models <- thin(fx_var_modellist(), thin = 2)

  expect_identical(nrow(vec[["posterior"]][["beta"]][["coeffs"]]),
                   as.integer(fx_iterations / 2))
  expect_s3_class(models, "modellist")
  expect_true(all(vapply(models,
                         function(x) nrow(x[["posterior"]][["a"]][["coeffs"]]),
                         integer(1)) == fx_iterations / 2))
})

test_that("window subsets the estimation sample in time", {
  model <- stats::window(fx_var_model(), start = c(1965, 1), end = c(1970, 4))
  y <- model[["data"]][["train"]][["y"]]

  expect_equal(stats::tsp(y)[1], 1965)
  expect_equal(stats::tsp(y)[2], 1970.75)
  expect_identical(nrow(y), 24L)
  # x is cut alongside y.
  expect_identical(nrow(model[["data"]][["train"]][["x"]]), nrow(y))
})

test_that("window cuts the SUR matrix in matching row blocks", {
  full <- fx_var_model()
  cut <- stats::window(full, start = c(1965, 1), end = c(1970, 4))
  k <- full[["model"]][["k"]]

  expect_identical(nrow(cut[["data"]][["train"]][["z"]]),
                   as.integer(k * nrow(cut[["data"]][["train"]][["y"]])))

  # The retained rows are the ones belonging to the retained periods.
  keep <- which(stats::time(full[["data"]][["train"]][["y"]]) %in%
                  stats::time(cut[["data"]][["train"]][["y"]]))
  rows <- rep((keep - 1) * k, each = k) + seq_len(k)
  expect_equal(cut[["data"]][["train"]][["z"]],
               full[["data"]][["train"]][["z"]][rows, ])
})

test_that("window works for VEC models and model lists", {
  vec <- stats::window(fx_vec_model(), start = c(1980, 1))
  models <- stats::window(fx_var_modellist(), start = c(1965, 1))

  expect_equal(stats::tsp(vec[["data"]][["train"]][["y"]])[1], 1980)
  expect_s3_class(models, "modellist")
  expect_true(all(vapply(
    models,
    function(x) stats::tsp(x[["data"]][["train"]][["y"]])[1],
    numeric(1)
  ) == 1965))
})
