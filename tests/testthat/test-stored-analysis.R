# Impulse responses and variance decompositions of a model read in pieces.
#
# What a stored model gives has to be what the same call on the model in memory
# gives: the pieces are a way of reading the chain, not a different estimator.

stored_pair <- function() {
  model <- fx_var_model()
  set.seed(9001)
  model <- add_posterior_coefficients(add_initial_values(
    add_priors(model, coef = list(v_i = 1), sigma = list(df = 3, scale = 1e-8))))

  file <- file.path(tempdir(), "bvartools-stored-analysis.h5")
  unlink(file)
  write_to_hdf5(model, filename = file)

  list(memory = model, file = open_model(file))
}

test_that("a stored model says what it holds", {
  pair <- stored_pair()
  stored <- pair[["file"]]

  expect_s3_class(stored, "bvarfile")
  expect_identical(stored[["draws"]],
                   nrow(pair[["memory"]][["posterior"]][["u_sigma_inv"]][["coeffs"]]))
  expect_identical(stored[["model"]][["k"]], pair[["memory"]][["model"]][["k"]])
  expect_equal(stored[["data"]], pair[["memory"]][["data"]])
  expect_output(print(stored), "draws of a")
  expect_error(open_model(file.path(tempdir(), "no-such-model.h5")), "does not exist")
})

test_that("map_draws walks the chain in pieces", {
  pair <- stored_pair()
  stored <- pair[["file"]]

  sums <- map_draws(stored, function(model) {
    colSums(unclass(model[["posterior"]][["a"]][["coeffs"]]))
  }, chunk = 7)
  expect_gt(length(sums), 1)
  expect_equal(Reduce(`+`, sums),
               colSums(unclass(pair[["memory"]][["posterior"]][["a"]][["coeffs"]])))
  expect_error(map_draws(pair[["memory"]], function(model) model), "'bvarfile'")
})

test_that("the impulse responses of a stored model are those of the model", {
  pair <- stored_pair()
  variables <- pair[["memory"]][["model"]][["endogen"]]

  for (type in c("feir", "oir", "gir")) {
    for (chunk in c(6, 1000)) {
      from_memory <- irf(pair[["memory"]], impulse = variables[1],
                         response = variables[2], n_ahead = 4, type = type)
      from_file <- irf(pair[["file"]], impulse = variables[1],
                       response = variables[2], n_ahead = 4, type = type,
                       chunk = chunk)
      expect_equal(from_file, from_memory, info = paste(type, chunk))
    }
  }

  # And the draws themselves, where they are asked for.
  draws_memory <- irf(pair[["memory"]], impulse = variables[1], response = variables[2],
                      n_ahead = 3, keep_draws = TRUE)
  draws_file <- irf(pair[["file"]], impulse = variables[1], response = variables[2],
                    n_ahead = 3, keep_draws = TRUE, chunk = 5)
  expect_equal(unclass(draws_file), unclass(draws_memory), ignore_attr = TRUE)
})

test_that("the variance decomposition of a stored model is that of the model", {
  pair <- stored_pair()
  variables <- pair[["memory"]][["model"]][["endogen"]]

  for (type in c("oir", "gir")) {
    for (chunk in c(6, 1000)) {
      from_memory <- fevd(pair[["memory"]], response = variables[2], n_ahead = 4,
                          type = type)
      from_file <- fevd(pair[["file"]], response = variables[2], n_ahead = 4,
                        type = type, chunk = chunk)
      expect_equal(from_file, from_memory, info = paste(type, chunk))
    }
  }

  # Normalising and collapsing groups describe the whole chain, not a piece.
  expect_equal(fevd(pair[["file"]], response = variables[2], n_ahead = 4, type = "gir",
                    normalise_gir = TRUE, chunk = 6),
               fevd(pair[["memory"]], response = variables[2], n_ahead = 4, type = "gir",
                    normalise_gir = TRUE))
  expect_equal(fevd(pair[["file"]], response = variables[2], n_ahead = 4, type = "oir",
                    max_groups = 2, chunk = 6),
               fevd(pair[["memory"]], response = variables[2], n_ahead = 4, type = "oir",
                    max_groups = 2))
})
