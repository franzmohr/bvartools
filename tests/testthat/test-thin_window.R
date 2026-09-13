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

test_that("thin keeps every block of draws in step", {
  keep <- seq(3, fx_iterations, 3)

  for (model in list(fx_var_tvp_fitted("gamma"), fx_vec_tvp_fitted("sv"))) {
    before <- draw_blocks(model[["posterior"]])
    after <- draw_blocks(thin(model, thin = 3)[["posterior"]])

    expect_named(after, names(before))
    for (name in names(before)) {
      expect_equal(unname(as.matrix(after[[name]])),
                   unname(as.matrix(before[[name]])[keep, , drop = FALSE]), label = name)
      expect_equal(attr(after[[name]], "mcpar"), c(3, max(keep), 3), label = name)
    }
  }

  # The blocks a list of names used to miss are among the ones checked.
  expect_true(all(c("a$sigma", "beta$rho") %in%
                    names(draw_blocks(fx_vec_tvp_fitted("sv")[["posterior"]]))))
})

test_that("window cuts the posterior paths to the periods it keeps", {
  specs <- list(list(model = fx_var_tvp_fitted("sv"), start = c(1965, 1), end = c(1975, 4)),
                list(model = fx_vec_tvp_fitted("sv"), start = c(1980, 1), end = c(1990, 4)))

  for (spec in specs) {
    full <- spec[["model"]]
    cut <- stats::window(full, start = spec[["start"]], end = spec[["end"]])
    k <- full[["model"]][["k"]]
    keep <- which(stats::time(full[["data"]][["train"]][["y"]]) %in%
                    stats::time(cut[["data"]][["train"]][["y"]]))
    columns <- function(width) rep((keep - 1) * width, each = width) + seq_len(width)

    widths <- c("a$coeffs" = ncol(full[["data"]][["train"]][["z"]]),
                "u_sigma_inv$coeffs" = k^2, "u_omega_inv$coeffs" = k, "loglik" = 1)
    if (!is.null(full[["posterior"]][["beta"]])) {
      widths[["beta$coeffs"]] <- ncol(full[["data"]][["train"]][["w"]]) * full[["model"]][["rank"]]
    }

    before <- draw_blocks(full[["posterior"]])
    after <- draw_blocks(cut[["posterior"]])
    for (name in names(widths)) {
      expect_equal(unname(as.matrix(after[[name]])),
                   unname(as.matrix(before[[name]])[, columns(widths[[name]]), drop = FALSE]),
                   label = name)
      expect_equal(attr(after[[name]], "mcpar"), attr(before[[name]], "mcpar"), label = name)
    }

    # Draws that do not vary by period are left as they are.
    for (name in setdiff(names(before), names(widths))) {
      expect_identical(after[[name]], before[[name]], label = name)
    }

    # A period of the window is the period of the original sample it came from.
    expect_equal(summary(cut)[["a"]][["means"]],
                 summary(full, period = max(keep))[["a"]][["means"]])
    expect_equal(summary(cut, period = 1)[["sigma"]][["means"]],
                 summary(full, period = min(keep))[["sigma"]][["means"]])
  }
})

test_that("the thinning interval is a positive integer within the draws", {
  model <- fx_at_vec_tvp()

  expect_error(thin(model, thin = fx_iterations + 1),
               paste0("only ", fx_iterations, " draws"))
  expect_error(thin(model, thin = 0), "single positive integer")
  expect_error(thin(model, thin = 2.5), "single positive integer")
  expect_error(thin(model, thin = c(2, 3)), "single positive integer")
  expect_identical(nrow(thin(model, thin = fx_iterations)[["posterior"]][["a"]][["coeffs"]]), 1L)
})
