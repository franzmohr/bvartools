test_that("add_initial_values() sets a seed that set.seed() reproduces", {
  set.seed(1)
  first <- add_initial_values(fx_var_priors())
  set.seed(1)
  again <- add_initial_values(fx_var_priors())

  expect_type(first[["model"]][["seed"]], "integer")
  expect_gte(first[["model"]][["seed"]], 0L)
  expect_identical(first[["model"]][["seed"]], again[["model"]][["seed"]])
})

test_that("add_initial_values() keeps a seed the model already has", {
  model <- add_seed(fx_var_priors(), 42)
  expect_identical(add_initial_values(model)[["model"]][["seed"]], 42L)

  vec <- add_seed(fx_vec_priors(), 43)
  expect_identical(add_initial_values(vec)[["model"]][["seed"]], 43L)
})

test_that("add_seed() stores a seed as an integer and refuses anything else", {
  expect_identical(add_seed(fx_var_initial(), 20260916)[["model"]][["seed"]], 20260916L)
  expect_identical(add_seed(fx_vec_initial(), 0)[["model"]][["seed"]], 0L)

  for (bad in list(-1, 1.5, NA, "1", c(1, 2), Inf, .Machine$integer.max + 1)) {
    expect_error(add_seed(fx_var_initial(), bad), "whole number")
  }
})

test_that("a list gets one seed per model, counted through nested lists", {
  models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                             iterations = 10, burnin = 5)
  models <- add_priors(models, coef = list(v_i = 0, v_i_det = 0),
                       sigma = list(df = 1, scale = 0.0001))

  seeded <- add_seed(models, 100)
  expect_identical(vapply(seeded, function(m) m[["model"]][["seed"]], integer(1)),
                   c(100L, 101L))

  # use_expanding_window() on a 'modellist' returns a 'modellist' of
  # 'expandingwindow' lists.
  nested <- add_seed(use_expanding_window(models, start = c(1997, 2)), 1)
  seeds <- unlist(lapply(nested, function(windows) {
    vapply(windows, function(m) m[["model"]][["seed"]], integer(1))
  }))
  expect_identical(seeds, seq_along(seeds))
})

test_that("models that are not estimated are skipped and not counted", {
  forecast <- structure(list(), class = "externalforecast")
  models <- structure(list(fx_var_initial(), forecast, fx_var_initial()),
                      class = c("modellist", "list"))

  seeded <- add_seed(models, 5)
  expect_identical(seeded[[1]][["model"]][["seed"]], 5L)
  expect_identical(seeded[[2]], forecast)
  expect_identical(seeded[[3]][["model"]][["seed"]], 6L)
  expect_identical(add_seed(forecast, 1), forecast)
})

test_that("a list seeds the models of other packages that have an add_seed() method", {
  # As dfmtools registers add_seed() methods for its dynamic factor models.
  registerS3method("add_seed", "bvartools_seed_test_model",
                   function(object, seed, ...) {
                     object[["model"]][["seed"]] <- as.integer(seed)
                     object
                   },
                   envir = asNamespace("bvartools"))
  other <- structure(list(model = list()),
                     class = c("bvartools_seed_test_model", "list"))
  unknown <- structure(list(model = list()),
                       class = c("bvartools_seed_test_unknown", "list"))

  models <- structure(list(fx_var_initial(), other, unknown, other),
                      class = c("modellist", "list"))
  seeded <- add_seed(models, 5)
  expect_identical(seeded[[1]][["model"]][["seed"]], 5L)
  expect_identical(seeded[[2]][["model"]][["seed"]], 6L)
  expect_identical(seeded[[3]], unknown)
  expect_identical(seeded[[4]][["model"]][["seed"]], 7L)
})

test_that("an external forecast in an expanding window list is not walked into", {
  forecast <- structure(list(fx_var_initial()),
                        class = c("externalforecast", "expandingwindow", "list"))
  models <- structure(list(forecast, fx_var_initial()),
                      class = c("modellist", "list"))

  seeded <- add_seed(models, 3)
  expect_identical(seeded[[1]], forecast)
  expect_identical(seeded[[2]][["model"]][["seed"]], 3L)
})

test_that("the seed of a model decides its draws and leaves R's generator as it was", {
  model <- add_seed(fx_var_initial(), 7)

  set.seed(1)
  state <- get(".Random.seed", envir = globalenv())
  first <- add_posterior_coefficients(model)
  expect_identical(get(".Random.seed", envir = globalenv()), state)

  set.seed(2)
  second <- add_posterior_coefficients(model)
  expect_equal(first[["posterior"]][["a"]][["coeffs"]],
               second[["posterior"]][["a"]][["coeffs"]])
})

test_that("a seeded model draws the same under another kind of generator", {
  model <- add_seed(fx_var_initial(), 7)
  reference <- add_posterior_coefficients(model)[["posterior"]][["a"]][["coeffs"]]

  # As on a cluster after parallel::clusterSetRNGStream().
  old <- RNGkind("L'Ecuyer-CMRG")
  set.seed(3)
  other <- add_posterior_coefficients(model)[["posterior"]][["a"]][["coeffs"]]
  kind_after <- RNGkind()[1]
  RNGkind(old[1], old[2], old[3])

  expect_equal(other, reference)
  expect_identical(kind_after, "L'Ecuyer-CMRG")
})

test_that("a model without a seed draws from R's generator as it stands", {
  model <- fx_var_initial()
  model[["model"]][["seed"]] <- NULL
  run <- function(state) {
    set.seed(state)
    add_posterior_coefficients(model)[["posterior"]][["a"]][["coeffs"]]
  }

  expect_equal(run(5), run(5))
  expect_false(isTRUE(all.equal(run(5), run(6))))
})

test_that("the seed survives a round trip through an HDF5 file", {
  var_file <- temp_h5_file()
  write_to_hdf5(add_seed(fx_var_initial(), 123), filename = var_file)
  expect_identical(read_model_from_hdf5(filename = var_file)[["model"]][["seed"]], 123L)

  vec_file <- temp_h5_file()
  write_to_hdf5(add_seed(fx_vec_initial(), 456), filename = vec_file)
  expect_identical(read_model_from_hdf5(filename = vec_file)[["model"]][["seed"]], 456L)
})

test_that("expanding windows count a model's seed up from window to window", {
  windows <- use_expanding_window(add_seed(fx_var_priors(), 10), start = c(1997, 2))
  seeds <- vapply(windows, function(m) m[["model"]][["seed"]], integer(1))
  expect_identical(seeds, 10L + seq_along(seeds) - 1L)

  vec_windows <- use_expanding_window(add_seed(fx_vec_priors(), 20), start = c(2004, 1))
  vec_seeds <- vapply(vec_windows, function(m) m[["model"]][["seed"]], integer(1))
  expect_identical(vec_seeds, 20L + seq_along(vec_seeds) - 1L)

  unseeded <- use_expanding_window(fx_var_priors(), start = c(1997, 2))
  expect_true(all(vapply(unseeded, function(m) is.null(m[["model"]][["seed"]]), logical(1))))
})

test_that("bayests_posterior() needs an executable that exists", {
  old_option <- options(bvartools.bayests_executable = NULL)
  old_env <- Sys.getenv("BAYESTS_EXECUTABLE", unset = NA)
  Sys.unsetenv("BAYESTS_EXECUTABLE")

  expect_error(bayests_posterior(), "No BayesTS executable")
  expect_error(bayests_posterior(executable = tempfile()), "not found")

  options(old_option)
  if (!is.na(old_env)) {
    Sys.setenv(BAYESTS_EXECUTABLE = old_env)
  }
})

# Runs only where BAYESTS_EXECUTABLE names a BayesTS build, and
# BAYESTS_LIBRARY_PATH, if needed, the directories of its runtime libraries.
bayests_for_tests <- function() {
  executable <- Sys.getenv("BAYESTS_EXECUTABLE")
  skip_if(!nzchar(executable) || !file.exists(executable),
          "BAYESTS_EXECUTABLE does not name a BayesTS executable")
  library_path <- Sys.getenv("BAYESTS_LIBRARY_PATH")
  library_path <- if (nzchar(library_path)) {
    strsplit(library_path, .Platform$path.sep, fixed = TRUE)[[1]]
  }
  bayests_posterior(executable = executable, library_path = library_path)
}

test_that("BayesTS draws have the structure of the internal ones", {
  simulate <- bayests_for_tests()

  for (model in list(add_seed(fx_var_initial(), 11), add_seed(fx_vec_initial(), 12))) {
    internal <- draw_blocks(add_posterior_coefficients(model)[["posterior"]])
    result <- add_posterior_coefficients(model, posterior_function = simulate)
    external <- draw_blocks(result[["posterior"]])

    expect_setequal(names(external), names(internal))
    for (block in names(internal)) {
      expect_identical(dim(external[[block]]), dim(internal[[block]]))
      expect_equal(attr(external[[block]], "mcpar"), attr(internal[[block]], "mcpar"))
    }
    # Only the posterior is added.
    expect_identical(result[["model"]], model[["model"]])
    expect_identical(result[["data"]], model[["data"]])
    expect_identical(result[["priors"]], model[["priors"]])
  }
})

test_that("BayesTS draws with the model's seed and gives an unseeded model one", {
  simulate <- bayests_for_tests()

  model <- add_seed(fx_var_initial(), 11)
  draws <- function() {
    add_posterior_coefficients(model, posterior_function = simulate)[["posterior"]][["a"]][["coeffs"]]
  }
  expect_equal(draws(), draws())

  unseeded <- fx_var_initial()
  unseeded[["model"]][["seed"]] <- NULL
  result <- add_posterior_coefficients(unseeded, posterior_function = simulate)
  expect_type(result[["model"]][["seed"]], "integer")
})

test_that("bayests_files() needs an executable and a path that exist", {
  old_option <- options(bvartools.bayests_executable = NULL)
  old_env <- Sys.getenv("BAYESTS_EXECUTABLE", unset = NA)
  Sys.unsetenv("BAYESTS_EXECUTABLE")

  expect_error(bayests_files(), "No BayesTS executable")
  expect_error(bayests_files(executable = tempfile()), "not found")

  options(old_option)
  if (!is.na(old_env)) {
    Sys.setenv(BAYESTS_EXECUTABLE = old_env)
  }

  executable <- Sys.getenv("BAYESTS_EXECUTABLE")
  skip_if(!nzchar(executable) || !file.exists(executable),
          "BAYESTS_EXECUTABLE does not name a BayesTS executable")
  run <- bayests_files(executable = executable)
  expect_error(run(file.path(tempdir(), "no-such-directory")), "No such file")
  expect_invisible(run(character(0)))
})

test_that("BayesTS run on the files draws what it draws on a model", {
  simulate <- bayests_for_tests()
  executable <- Sys.getenv("BAYESTS_EXECUTABLE")
  library_path <- Sys.getenv("BAYESTS_LIBRARY_PATH")
  library_path <- if (nzchar(library_path)) {
    strsplit(library_path, .Platform$path.sep, fixed = TRUE)[[1]]
  }
  run <- bayests_files(executable = executable, library_path = library_path)

  for (model in list(add_seed(fx_var_initial(), 11), add_seed(fx_vec_initial(), 12))) {

    # The same model, once drawn through the session and once in its file.
    in_session <- add_posterior_coefficients(model, posterior_function = simulate)

    folder <- file.path(tempdir(), "bvartools-bayests-files")
    unlink(folder, recursive = TRUE)
    dir.create(folder, recursive = TRUE)
    write_to_hdf5(model, filename = file.path(folder, "model.h5"))

    expect_invisible(run(folder, command = "check"))
    run(folder, args = c("--no-loglik", "--no-forecasts"))
    in_file <- read_model_from_hdf5(file.path(folder, "model.h5"))

    expect_false(is.null(in_file[["posterior"]]))
    for (block in names(in_session[["posterior"]])) {
      first <- in_session[["posterior"]][[block]]
      second <- in_file[["posterior"]][[block]]
      if (is.list(first)) {
        first <- first[["coeffs"]]
        second <- second[["coeffs"]]
      }
      expect_equal(unclass(second), unclass(first), info = block)
    }
    unlink(folder, recursive = TRUE)
  }
})

test_that("BayesTS runs several paths at once and reports the ones that fail", {
  bayests_for_tests()
  executable <- Sys.getenv("BAYESTS_EXECUTABLE")
  library_path <- Sys.getenv("BAYESTS_LIBRARY_PATH")
  library_path <- if (nzchar(library_path)) {
    strsplit(library_path, .Platform$path.sep, fixed = TRUE)[[1]]
  }
  run <- bayests_files(executable = executable, library_path = library_path)

  # Three directories of one model each, checked two at a time.
  root <- file.path(tempdir(), "bvartools-bayests-jobs")
  unlink(root, recursive = TRUE)
  folders <- file.path(root, c("one", "two", "three"))
  for (i in seq_along(folders)) {
    dir.create(folders[i], recursive = TRUE)
    write_to_hdf5(add_seed(fx_var_initial(), 20 + i),
                  filename = file.path(folders[i], "model.h5"))
  }

  expect_invisible(run(folders, command = "check", jobs = 2))
  run(folders, args = c("--no-loglik", "--no-forecasts"), jobs = 2, poll = 0.2)
  for (folder in folders) {
    drawn <- read_model_from_hdf5(file.path(folder, "model.h5"))
    expect_false(is.null(drawn[["posterior"]]))
  }

  # A path that holds nothing BayesTS can read fails, and says which.
  empty <- file.path(root, "empty")
  dir.create(empty)
  writeLines("not a model", file.path(empty, "model.h5"))
  expect_error(run(c(folders[1], empty), command = "check", jobs = 2),
               "BayesTS failed on 1 path")

  unlink(root, recursive = TRUE)
})
