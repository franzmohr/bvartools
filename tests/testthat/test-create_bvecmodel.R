test_that("a bvecmodel records the cointegration specification", {
  model <- fx_vec_model()

  expect_s3_class(model, "bvecmodel")
  expect_named(model, c("model", "data"))
  expect_identical(model[["model"]][["type"]], "VEC")
  expect_identical(model[["model"]][["k"]], 2L)
  expect_identical(model[["model"]][["rank"]], 1L)
  expect_identical(model[["model"]][["endogen"]], c("lr", "Dp"))
})

test_that("differences and lagged levels are built from the data", {
  model <- fx_vec_model()
  train <- model[["data"]][["train"]]
  endogen <- model[["data"]][["original"]][["endogen"]]
  nobs <- nrow(endogen) - model[["model"]][["p"]]

  expect_equal(nrow(train[["y"]]), nobs)
  expect_identical(colnames(train[["y"]]), c("d.lr", "d.Dp"))
  expect_identical(colnames(train[["w"]]), c("l.lr", "l.Dp"))

  strip <- function(x) matrix(as.numeric(x), nrow = NROW(x))
  # y holds the first differences and w the levels lagged once.
  expect_equal(strip(train[["y"]]), strip(diff(endogen)))
  expect_equal(strip(train[["w"]]), strip(endogen[-nrow(endogen), ]))
})

test_that("restricted deterministic terms join the cointegration term", {
  restricted <- create_bvecmodel(vec_data(), p = 2, r = 1, const = "restricted",
                                 iterations = 10, burnin = 5)
  unrestricted <- create_bvecmodel(vec_data(), p = 2, r = 1,
                                   const = "unrestricted",
                                   iterations = 10, burnin = 5)

  # A restricted constant is part of beta, so it widens the cointegration term.
  expect_identical(restricted[["model"]][["n_restricted"]], 1L)
  expect_identical(restricted[["model"]][["k_beta"]],
                   restricted[["model"]][["k"]] + 1L)
  expect_true("const" %in% colnames(restricted[["data"]][["train"]][["w"]]))

  # An unrestricted constant is an ordinary regressor.
  expect_identical(unrestricted[["model"]][["n_restricted"]], 0L)
  expect_identical(unrestricted[["model"]][["k_beta"]],
                   unrestricted[["model"]][["k"]])
  expect_identical(unrestricted[["model"]][["n"]], 1L)
  expect_true("const" %in% colnames(unrestricted[["data"]][["train"]][["x"]]))
})

test_that("lagged differences are added for p greater than one", {
  model <- create_bvecmodel(vec_data(), p = 3, r = 1, const = "unrestricted",
                            iterations = 10, burnin = 5)

  # p is the lag order of the model in levels, so p - 1 lagged differences enter.
  expect_identical(colnames(model[["data"]][["train"]][["x"]]),
                   c("d.lr.l01", "d.Dp.l01", "d.lr.l02", "d.Dp.l02", "const"))
})

test_that("several ranks produce a modellist", {
  models <- create_bvecmodel(vec_data(), p = 1, r = 0:2, const = "unrestricted",
                             iterations = 10, burnin = 5)

  expect_s3_class(models, "modellist")
  expect_length(models, 3)
  expect_identical(vapply(models, function(x) x[["model"]][["rank"]], integer(1)),
                   0:2)
})

test_that("invalid input is rejected", {
  expect_error(create_bvecmodel(unclass(vec_data()), p = 1, r = 1),
               "must be an object of class 'ts'")
  expect_error(create_bvecmodel(vec_data(), p = 1, r = 1, const = "wrong"),
               "not valid")
  expect_error(create_bvecmodel(vec_data(), p = 1, r = 1, trend = "wrong"),
               "not valid")
})

test_that("lag orders of the VAR in levels must be given", {
  expect_error(create_bvecmodel(vec_data(), r = 1), "Argument 'p' must be specified")
  expect_error(create_bvecmodel(vec_data(), p = 0, r = 1), "at least 1")
  expect_error(create_bvecmodel(vec_data(), p = 1.5, r = 1), "integer vector")

  data <- vec_data()
  exo <- stats::window(at_macrodata[["foreign"]][, "r.s", drop = FALSE], end = c(2005, 4))
  expect_error(create_bvecmodel(data, p = 2, exogen = exo, r = 1),
               "Argument 's' must be specified")
  expect_error(create_bvecmodel(data, p = 2, exogen = exo, s = 0, r = 1), "at least 1")

  # Without exogenous variables s plays no role.
  expect_s3_class(create_bvecmodel(vec_data(), p = 1, s = 0, r = 1,
                                   iterations = 10, burnin = 5), "bvecmodel")

  # s = 1 is the contemporaneous difference of the exogenous variables only.
  model <- create_bvecmodel(data, p = 2, exogen = exo, s = 1, r = 1,
                            iterations = 10, burnin = 5)
  expect_identical(model[["model"]][["s"]], 1L)
  expect_identical(colnames(model[["data"]][["train"]][["x"]]),
                   c("d.lr.l01", "d.Dp.l01", "d.r.s.l00"))
})
