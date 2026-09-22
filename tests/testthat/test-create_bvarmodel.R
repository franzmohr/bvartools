test_that("a single lag order produces one bvarmodel", {
  model <- fx_var_model()

  expect_s3_class(model, "bvarmodel")
  expect_named(model, c("model", "data"))
  expect_identical(model[["model"]][["type"]], "VAR")
  expect_identical(model[["model"]][["k"]], 3L)
  expect_identical(model[["model"]][["p"]], 1L)
  expect_identical(model[["model"]][["endogen"]], c("y", "Dp", "r"))
  expect_identical(model[["model"]][["iterations"]], fx_iterations)
  expect_identical(model[["model"]][["burnin"]], fx_burnin)
})

test_that("several lag orders produce a modellist", {
  models <- create_bvarmodel(var_data(), p = 0:2, deterministic = "const",
                             iterations = 10, burnin = 5)

  expect_s3_class(models, "modellist")
  expect_length(models, 3)
  expect_true(all(vapply(models, inherits, logical(1), "bvarmodel")))
  expect_identical(vapply(models, function(x) x[["model"]][["p"]], integer(1)),
                   0:2)
})

test_that("data matrices have the dimensions implied by the specification", {
  model <- fx_var_model()
  spec <- model[["model"]]
  train <- model[["data"]][["train"]]

  nobs <- nrow(model[["data"]][["original"]][["endogen"]]) - spec[["p"]]
  n_regressors <- spec[["k"]] * spec[["p"]] + spec[["n"]]

  expect_identical(dim(train[["y"]]), c(nobs, spec[["k"]]))
  expect_identical(dim(train[["x"]]), c(nobs, n_regressors))
  # z is the SUR form of x: one row block per observation, one column block per
  # equation.
  expect_identical(dim(train[["z"]]),
                   c(spec[["k"]] * nobs, spec[["k"]] * n_regressors))
})

test_that("y and x are aligned in time", {
  model <- fx_var_model()
  y <- model[["data"]][["train"]][["y"]]
  x <- model[["data"]][["train"]][["x"]]
  endogen <- model[["data"]][["original"]][["endogen"]]

  expect_identical(stats::tsp(y), stats::tsp(x))
  # The first lag block of x is the endogenous series shifted by one period.
  strip <- function(x) matrix(as.numeric(x), nrow = NROW(x))
  expect_equal(strip(x[, 1:3]), strip(endogen[-nrow(endogen), ]))
  expect_equal(strip(y), strip(endogen[-1, ]))
})

test_that("regressors are named after variable and lag", {
  model <- fx_var_model()

  expect_identical(colnames(model[["data"]][["train"]][["x"]]),
                   c("y.01", "Dp.01", "r.01", "const"))
  expect_identical(colnames(model[["data"]][["train"]][["y"]]),
                   c("y", "Dp", "r"))
})

test_that("deterministic terms drive the number of deterministic regressors", {
  n_deterministic <- function(...) {
    create_bvarmodel(var_data(), p = 1, iterations = 10, burnin = 5,
                     ...)[["model"]][["n"]]
  }

  expect_identical(n_deterministic(deterministic = "none"), 0L)
  expect_identical(n_deterministic(deterministic = "const"), 1L)
  expect_identical(n_deterministic(deterministic = "trend"), 1L)
  expect_identical(n_deterministic(deterministic = "both"), 2L)
  # Quarterly data, so seasonal dummies add frequency - 1 columns.
  expect_identical(n_deterministic(deterministic = "const", seasonal = TRUE), 4L)
})

test_that("exogenous variables enter with the requested number of lags", {
  endogen <- var_data()
  exogen <- endogen[, 1, drop = FALSE]
  colnames(exogen) <- "ex"
  exogen <- stats::ts(exogen, start = stats::start(endogen),
                      frequency = stats::frequency(endogen))

  model <- create_bvarmodel(endogen[, -1], p = 1, exogen = exogen, s = 1,
                            deterministic = "const",
                            iterations = 10, burnin = 5)

  expect_identical(model[["model"]][["m"]], 1L)
  expect_identical(model[["model"]][["s"]], 1L)
  expect_true(any(grepl("^ex", colnames(model[["data"]][["train"]][["x"]]))))
})

test_that("model options are recorded in the specification", {
  model <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                            structural = TRUE, tvp = TRUE, error = "sv",
                            varsel = "ssvs", iterations = 10, burnin = 5)

  expect_true(model[["model"]][["structural"]])
  expect_true(model[["model"]][["tvp"]])
  expect_identical(model[["model"]][["error"]], "sv")
  expect_identical(model[["model"]][["varsel"]], "ssvs")
})

test_that("invalid input is rejected", {
  expect_error(create_bvarmodel(unclass(var_data()), p = 1),
               "must be an object of class 'ts'")
  expect_error(create_bvarmodel(var_data(), p = 1, exogen = 1:10),
               "must be an object of class 'ts'")
  expect_error(create_bvarmodel(var_data(), p = 1, deterministic = "none",
                                seasonal = TRUE),
               "must be either 'const' or 'both'")
})

test_that("exogenous variables enter with two lags by default", {
  data <- bvartools::at_macrodata
  model <- create_bvarmodel(data = data[["domestic"]], exogen = data[["foreign"]],
                            iterations = 10, burnin = 5)

  expect_s3_class(model, "bvarmodel")
  expect_identical(model[["model"]][["m"]], NCOL(data[["foreign"]]))
  expect_identical(model[["model"]][["s"]], 2L)
  x_names <- colnames(model[["data"]][["train"]][["x"]])
  expect_true(all(paste0("poil.l", c("00", "01", "02")) %in% x_names))
  expect_false("poil.l03" %in% x_names)
})

test_that("an invalid lag order of exogenous variables is rejected", {
  data <- bvartools::at_macrodata
  for (s in list(NULL, -1, 1.5, NA_real_)) {
    expect_error(create_bvarmodel(data = data[["domestic"]], exogen = data[["foreign"]],
                                  s = s, iterations = 10, burnin = 5),
                 "Argument 's' must contain non-negative integers")
  }
})

test_that("the trend does not depend on where the exogenous series starts", {
  # The trend counted rows of the data and 'exogen' together, so an exogenous
  # series reaching further back shifted it.
  endogen <- stats::window(at_macrodata[["domestic"]][, c("y", "Dp")], start = c(1998, 1)) * 100
  foreign <- at_macrodata[["foreign"]][, "Dp.s", drop = FALSE] * 100
  short <- stats::window(foreign, start = c(1998, 1))

  x_trend <- function(exogen) {
    model <- create_bvarmodel(endogen, p = 1, exogen = exogen, s = 1,
                              deterministic = "both", iterations = 10, burnin = 5)
    as.numeric(model[["data"]][["train"]][["x"]][, "trend"])
  }
  expect_identical(x_trend(foreign), x_trend(short))
  expect_identical(x_trend(short)[1], 1)
})

test_that("the lag order is validated", {
  expect_error(create_bvarmodel(var_data(), p = 1.5), "non-negative whole numbers")
  expect_error(create_bvarmodel(var_data(), p = -1), "non-negative whole numbers")
  expect_error(create_bvarmodel(var_data(), p = NA), "non-negative whole numbers")
})
