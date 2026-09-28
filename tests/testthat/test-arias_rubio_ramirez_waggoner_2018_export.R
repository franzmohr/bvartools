# The exported identification step, arias_rubio_ramirez_waggoner_2018(). It is
# what add_sign_zero_restrictions() runs, so the checks are that the two agree
# under the same seed, and that what only the exported function allows -- a
# table of signs alone, draws that did not come from a 'bvarmodel' -- works.

ex_model <- function() {
  cached_fixture("szr_model", {
    object <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                               iterations = 60, burnin = 10)
    object <- add_priors(object, coef = list(v_i = 1, v_i_det = 0.1),
                         sigma = list(df = "k", scale = 1))
    object <- add_initial_values(object)
    add_posterior_coefficients(add_seed(object, 216))
  })
}

ex_draws <- function(object) {
  bvartools:::.collect_draws(object, need_A0 = FALSE, need_Sigma = TRUE,
                             all_regressors = TRUE)
}

test_that("it identifies the draws add_sign_zero_restrictions() identifies", {
  object <- ex_model()
  restrictions <- data.frame(impulse = "y", response = c("Dp", "r"), sign = c(0, 1))

  set.seed(99)
  direct <- arias_rubio_ramirez_waggoner_2018(ex_draws(object), restrictions,
                                              object[["model"]][["endogen"]], lags = 1)
  set.seed(99)
  method <- add_sign_zero_restrictions(object, restrictions)
  record <- method[["model"]][["sign_zero_restrictions"]]

  expect_identical(direct[["accepted"]], record[["accepted"]])
  expect_identical(direct[["effective_sample_size"]], record[["effective_sample_size"]])
  expect_identical(direct[["pareto_k"]], record[["pareto_k"]])
  expect_equal(sum(direct[["weights"]]), 1)
  expect_identical(which(direct[["weights"]] > 0), which(!is.na(direct[["q"]][, 1])))
})

test_that("the rotation maps the Choleski factor into the restricted impact", {
  object <- ex_model()
  restrictions <- data.frame(impulse = "y", response = c("Dp", "r"), sign = c(0, 1))
  draws <- ex_draws(object)
  set.seed(3)
  out <- arias_rubio_ramirez_waggoner_2018(draws, restrictions,
                                           object[["model"]][["endogen"]], lags = 1)
  for (i in which(!is.na(out[["q"]][, 1]))[1:5]) {
    impact <- t(chol(draws[[i]][["Sigma"]])) %*% matrix(out[["q"]][i, ], 3)
    expect_equal(impact[2, 1], 0, tolerance = 1e-10)
    expect_gt(impact[3, 1], 0)
  }
})

test_that("a table of signs alone is accepted and weighted equally", {
  object <- ex_model()
  set.seed(4)
  out <- arias_rubio_ramirez_waggoner_2018(
    ex_draws(object), data.frame(impulse = "y", response = "y", sign = 1),
    object[["model"]][["endogen"]], lags = 1)
  kept <- out[["weights"]][out[["weights"]] > 0]
  expect_equal(kept, rep(1 / length(kept), length(kept)))
  # The method keeps its own refusal, since add_sign_restrictions() is cheaper.
  expect_error(add_sign_zero_restrictions(object, data.frame(impulse = "y", response = "y",
                                                             sign = 1)),
               "add_sign_restrictions")
})

test_that("the draws and the arguments are checked", {
  object <- ex_model()
  draws <- ex_draws(object)
  v <- object[["model"]][["endogen"]]
  r <- data.frame(impulse = "y", response = "Dp", sign = 0)
  expect_error(arias_rubio_ramirez_waggoner_2018(list(), r, v, 1), "non-empty list")
  expect_error(arias_rubio_ramirez_waggoner_2018(draws, r, v[1:2], 1), "one row and column")
  expect_error(arias_rubio_ramirez_waggoner_2018(draws, r, v, 5), "lag blocks")
  expect_error(arias_rubio_ramirez_waggoner_2018(draws, r, v, 1, max_tries = 0),
               "single positive integer")
  expect_error(arias_rubio_ramirez_waggoner_2018(draws, data.frame(impulse = "z",
                                                                   response = "y", sign = 1),
                                              v, 1),
               "does not contain")
})
