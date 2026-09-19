# Working on a folder of models one at a time.
#
# What the methods for a 'bvarfolder' have to deliver is what the same steps on
# a list of models in the session deliver: the same priors, the same starting
# values, the same seeds and so the same draws, with nothing but one model per
# worker in memory while they run.

folder_models <- function(type = "var") {
  data("e1", envir = environment())
  train <- diff(log(e1)) * 100
  if (type == "var") {
    create_bvarmodel(data = train, p = 1:2, deterministic = "const",
                     iterations = 20, burnin = 5)
  } else {
    create_bvecmodel(data = train, p = 1, r = 0:1, const = "unrestricted",
                     iterations = 20, burnin = 5)
  }
}

folder_priors <- function(object, type) {
  if (type == "var") {
    add_priors(object, coef = list(v_i = 1), sigma = list(df = 3, scale = 1))
  } else {
    add_priors(object, coef = list(v_i = 0), coint = list(v_i = 0, p_tau_i = 1),
               sigma = list(df = 3, scale = 1))
  }
}

written_folder <- function(object, name) {
  folder <- file.path(tempdir(), paste0("bvartools-folder-", name))
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  write_to_hdf5(object, folder = folder)
  open_models(folder)
}

test_that("a folder is opened with a row per model and no draws", {

  object <- folder_models("var")
  stored <- written_folder(object, "open")

  expect_s3_class(stored, "bvarfolder")
  manifest <- stored[["manifest"]]
  expect_equal(nrow(manifest), length(object))
  expect_setequal(manifest[["p"]], vapply(object, function(m) m[["model"]][["p"]], numeric(1)))
  expect_true(all(manifest[["iterations"]] == 20))
  expect_true(all(manifest[["group"]] == ""))
  expect_output(print(stored), "models in")

  # The files, which is what an external sampler is pointed at.
  files <- model_files(stored)
  expect_equal(length(files), nrow(manifest))
  expect_true(all(file.exists(files)))
  expect_equal(length(model_files(stored, models = manifest[["model"]][1])), 1)
  expect_error(model_files(stored, models = "nothing"), "does not hold the model")

  expect_error(open_models(file.path(tempdir(), "no-such-folder")), "does not exist")
})

test_that("a folder is estimated to the posterior of the models in the session", {
  for (type in c("var", "vec")) {

    object <- folder_models(type)
    stored <- written_folder(object, paste0("estimate-", type))

    in_session <- folder_priors(object, type)
    in_session <- add_seed(add_initial_values(in_session), 1234)
    in_session <- add_posterior_loglik(add_posterior_coefficients(in_session))

    on_disk <- folder_priors(stored, type)
    on_disk <- add_initial_values(on_disk)
    on_disk <- add_seed(on_disk, 1234)
    on_disk <- add_posterior_coefficients(on_disk)
    on_disk <- add_posterior_loglik(on_disk)
    expect_s3_class(on_disk, "bvarfolder")

    back <- read_models_from_folder(on_disk[["folder"]])
    expect_equal(length(back), length(in_session))
    for (i in seq_along(in_session)) {
      expect_equal(back[[i]][["model"]][["seed"]],
                   in_session[[i]][["model"]][["seed"]], info = i)
      expect_equal(unclass(back[[i]][["posterior"]][["a"]][["coeffs"]]),
                   unclass(in_session[[i]][["posterior"]][["a"]][["coeffs"]]),
                   info = i)
      expect_equal(unclass(back[[i]][["posterior"]][["loglik"]]),
                   unclass(in_session[[i]][["posterior"]][["loglik"]]), info = i)
    }
  }
})

test_that("the draws do not depend on the number of workers", {
  skip_on_cran()

  object <- folder_models("var")
  one <- written_folder(object, "cores-one")
  two <- written_folder(object, "cores-two")

  estimate <- function(stored, cores) {
    stored <- add_priors(stored, coef = list(v_i = 1),
                         sigma = list(df = 3, scale = 1), cores = cores)
    stored <- add_initial_values(stored, cores = cores)
    stored <- add_seed(stored, 99, cores = cores)
    add_posterior_coefficients(stored, cores = cores)
  }
  one <- estimate(one, 1)
  two <- estimate(two, 2)

  first <- read_models_from_folder(one[["folder"]])
  second <- read_models_from_folder(two[["folder"]])
  for (i in seq_along(first)) {
    expect_equal(unclass(second[[i]][["posterior"]][["a"]][["coeffs"]]),
                 unclass(first[[i]][["posterior"]][["a"]][["coeffs"]]), info = i)
  }
})

test_that("a step over one model at a time is the step over all of them", {
  skip_on_cran()

  object <- folder_models("var")
  whole <- written_folder(object, "subset-whole")
  piecewise <- written_folder(object, "subset-piecewise")

  estimate <- function(stored, models) {
    stored <- add_priors(stored, coef = list(v_i = 1),
                         sigma = list(df = 3, scale = 1), models = models)
    stored <- add_initial_values(stored, models = models)
    stored <- add_seed(stored, 99, models = models)
    add_posterior_coefficients(stored, models = models)
  }

  whole <- estimate(whole, NULL)
  for (model in piecewise[["manifest"]][["model"]]) {
    piecewise <- estimate(piecewise, model)
  }

  first <- read_models_from_folder(whole[["folder"]])
  second <- read_models_from_folder(piecewise[["folder"]])
  for (i in seq_along(first)) {
    expect_equal(unclass(second[[i]][["posterior"]][["a"]][["coeffs"]]),
                 unclass(first[[i]][["posterior"]][["a"]][["coeffs"]]), info = i)
  }
  expect_equal(open_models(piecewise[["folder"]])[["manifest"]],
               open_models(whole[["folder"]])[["manifest"]])
})

test_that("a step leaves the models it was not asked for alone", {

  object <- folder_models("var")
  stored <- written_folder(object, "untouched")
  names_of <- stored[["manifest"]][["model"]]

  before <- file.info(model_files(stored))[["mtime"]]
  Sys.sleep(1)
  stored <- add_priors(stored, coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 1), models = names_of[1])
  after <- file.info(model_files(stored))[["mtime"]]

  expect_true(after[1] > before[1])
  expect_equal(after[-1], before[-1])

  back <- read_models_from_folder(stored[["folder"]])
  expect_false(is.null(back[[1]][["priors"]]))
  expect_null(back[[2]][["priors"]])

  expect_error(map_models(stored, add_initial_values, models = "nothing"),
               "does not hold the model")
  expect_error(map_models(unclass(stored), add_initial_values), "must be of class")
  expect_error(map_models(stored, function(model) "not a model"),
               "rather than a model")
})

test_that("a folder is read without its draws and a piece at a time", {

  object <- folder_models("var")
  stored <- written_folder(object, "partial")
  stored <- add_priors(stored, coef = list(v_i = 1), sigma = list(df = 3, scale = 1))
  stored <- add_posterior_coefficients(add_initial_values(stored))

  none <- read_models_from_folder(stored[["folder"]], draws = integer(0))
  expect_equal(length(none), nrow(stored[["manifest"]]))
  expect_equal(nrow(none[[1]][["posterior"]][["a"]][["coeffs"]]), 0)
  expect_false(is.null(none[[1]][["priors"]]))

  piece <- read_models_from_folder(stored[["folder"]], draws = 1:5)
  expect_equal(nrow(piece[[1]][["posterior"]][["a"]][["coeffs"]]), 5)

  whole <- read_models_from_folder(stored[["folder"]])
  expect_equal(unclass(piece[[1]][["posterior"]][["a"]][["coeffs"]]),
               unclass(whole[[1]][["posterior"]][["a"]][["coeffs"]])[1:5, , drop = FALSE],
               ignore_attr = TRUE)
})

test_that("a folder is read without being written to", {

  object <- folder_models("var")
  stored <- written_folder(object, "read-only")
  stored <- add_priors(stored, coef = list(v_i = 1), sigma = list(df = 3, scale = 1))
  stored <- add_posterior_coefficients(add_initial_values(stored))
  stored <- add_posterior_loglik(stored)

  before <- file.info(model_files(stored))[["mtime"]]
  Sys.sleep(1)

  sizes <- map_models(stored, function(model) nrow(model[["data"]][["train"]][["y"]]),
                      write = FALSE)
  expect_identical(names(sizes), stored[["manifest"]][["model"]])
  expect_true(all(unlist(sizes) > 0))
  expect_equal(file.info(model_files(stored))[["mtime"]], before)

  # What a list of models in the session gives, from the files instead.
  in_session <- add_priors(object, coef = list(v_i = 1),
                           sigma = list(df = 3, scale = 1))
  in_session <- add_posterior_coefficients(add_initial_values(in_session))
  in_session <- add_posterior_loglik(in_session)
  criteria <- selection_criteria(stored)
  expect_s3_class(criteria, "selcritlist")
  expect_equal(length(criteria), length(in_session))
  best <- choose_best_model(criteria)
  expect_true(length(best) > 0)
})
