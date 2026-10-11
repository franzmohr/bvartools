# Student t errors of VarNormalWishart: create_bvarmodel(student_t = TRUE) and
# sigma$student_t in add_priors(). The sampler is BayesTS's and is tested there
# against its conditionals; what is checked here is that the R side reaches it,
# that its draws come back where they are documented, and that a period with
# outlying shocks is the one the mixing scales single out.

student_t_data <- function() {
  data <- var_data()
  # Shocks fifteen times their size in one period, as BayesTS's own test does.
  data[40, ] <- data[40, ] + 15 * apply(data, 2, stats::sd)
  data
}

student_t_model <- function(prior = list(shape = 2, rate = 0.1), data = student_t_data(),
                            iterations = 200, burnin = 100) {
  model <- create_bvarmodel(data, p = 1, deterministic = "const", student_t = TRUE,
                            iterations = iterations, burnin = burnin)
  add_priors(model, coef = list(v_i = 1, v_i_det = 0.1),
             sigma = list(df = "k", scale = 1, student_t = prior))
}

student_t_fit <- function() {
  cached_fixture("student_t_fit", {
    add_posterior_coefficients(add_seed(add_initial_values(student_t_model()), 4))
  })
}

test_that("Student t errors are offered where the sampler draws them and refused elsewhere", {
  expect_true(create_bvarmodel(var_data(), p = 1, student_t = TRUE)[["model"]][["student_t"]])
  expect_null(create_bvarmodel(var_data(), p = 1)[["model"]][["student_t"]])
  expect_error(create_bvarmodel(var_data(), p = 1, error = "gamma", student_t = TRUE),
               "wishart")
  expect_error(create_bvarmodel(var_data(), p = 1, tvp = TRUE, student_t = TRUE), "constant")
  expect_error(create_bvarmodel(var_data(), p = 1, student_t = NA), "TRUE or FALSE")
})

test_that("the prior of the degrees of freedom has to be stated, once", {
  model <- create_bvarmodel(var_data(), p = 1, student_t = TRUE)
  coef <- list(v_i = 1, v_i_det = 0.1)
  expect_error(add_priors(model, coef = coef, sigma = list(df = "k", scale = 1)),
               "sigma\\$student_t")
  expect_error(add_priors(model, coef = coef,
                          sigma = list(df = "k", scale = 1,
                                       student_t = list(df = 5, shape = 2, rate = 0.1))),
               "not both")
  expect_error(add_priors(model, coef = coef,
                          sigma = list(df = "k", scale = 1, student_t = list(shape = 2))),
               "rate")
  expect_error(add_priors(model, coef = coef,
                          sigma = list(df = "k", scale = 1, student_t = list(df = -1))),
               "positive")
  expect_identical(add_priors(model, coef = coef,
                              sigma = list(df = "k", scale = 1, student_t = list(df = 5)))
                   [["priors"]][["student_t"]], list(df = 5))

  # A model with normal errors, VAR or VEC, is refused the prior.
  normal <- create_bvarmodel(var_data(), p = 1)
  expect_error(add_priors(normal, coef = coef,
                          sigma = list(df = "k", scale = 1, student_t = list(df = 5))),
               "student_t = TRUE")
  vec <- create_bvecmodel(vec_data(), p = 2, r = 1)
  expect_error(add_priors(vec, coef = coef, coint = list(v_i = 0, p_tau_i = 1),
                          sigma = list(df = "k", scale = 1, student_t = list(df = 5))),
               "student_t = TRUE")
})

test_that("the mixing scales single out the period whose shocks were blown up", {
  fit <- student_t_fit()
  lambda <- as.matrix(fit[["posterior"]][["student_t"]][["lambda"]])
  tt <- nrow(fit[["data"]][["train"]][["y"]])
  expect_s3_class(fit[["posterior"]][["student_t"]][["lambda"]], "mcmc")
  expect_identical(dim(lambda), c(200L, tt))
  expect_identical(dim(as.matrix(fit[["posterior"]][["student_t"]][["df"]])), c(200L, 1L))

  # Period 40 of the data is period 39 of the sample after one lag.
  medians <- apply(lambda, 2, stats::median)
  expect_identical(unname(which.max(medians)), 39L)
  expect_gt(medians[39], 5 * stats::median(medians))
})

test_that("held degrees of freedom stay where they were put", {
  model <- student_t_model(prior = list(df = 7), iterations = 20, burnin = 10)
  fit <- add_posterior_coefficients(add_seed(add_initial_values(model), 2))
  expect_true(all(as.matrix(fit[["posterior"]][["student_t"]][["df"]]) == 7))
})

test_that("the likelihood, forecasts, window, thinning and model file carry the t errors", {
  fit <- student_t_fit()
  fit <- add_posterior_loglik(fit)
  expect_true(all(is.finite(as.matrix(fit[["posterior"]][["loglik"]]))))
  forecast <- add_posterior_forecasts(add_forecast_input(fit, n_ahead = 3))
  expect_identical(dim(forecast[["posterior"]][["forecast"]][["forecasts"]]), c(200L, 9L))

  y <- fit[["data"]][["train"]][["y"]]
  windowed <- window(fit, start = stats::time(y)[11])
  expect_identical(ncol(windowed[["posterior"]][["student_t"]][["lambda"]]),
                   nrow(windowed[["data"]][["train"]][["y"]]))
  expect_identical(nrow(thin(fit, 2)[["posterior"]][["student_t"]][["df"]]), 100L)

  path <- temp_h5_file()
  write_to_hdf5(fit, filename = path)
  back <- read_model_from_hdf5(path)
  expect_equal(back[["priors"]][["student_t"]][["shape"]], 2)
  expect_equal(unclass(as.matrix(back[["posterior"]][["student_t"]][["lambda"]])),
               unclass(as.matrix(fit[["posterior"]][["student_t"]][["lambda"]])),
               ignore_attr = TRUE)
})

test_that("a scenario is refused under Student t errors", {
  fit <- add_forecast_input(student_t_fit(), n_ahead = 2)
  scenario <- data.frame(period = 1, variable = 1, value = 0)
  expect_error(add_posterior_forecasts(fit, scenario = scenario))
})
