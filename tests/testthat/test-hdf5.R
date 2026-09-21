skip_if_not_installed("hdf5r")

test_that("a VAR model survives a write and read round trip", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path)

  expect_true(file.exists(path))
  restored <- read_model_from_hdf5(path)

  expect_s3_class(restored, "bvarmodel")
  expect_equal(unclass(restored[["posterior"]][["a"]][["coeffs"]]),
               unclass(fx_var_fitted()[["posterior"]][["a"]][["coeffs"]]),
               ignore_attr = TRUE)
  expect_equal(unclass(restored[["posterior"]][["u_sigma_inv"]][["coeffs"]]),
               unclass(fx_var_fitted()[["posterior"]][["u_sigma_inv"]][["coeffs"]]),
               ignore_attr = TRUE)
})

test_that("the model specification survives the round trip", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path)
  restored <- read_model_from_hdf5(path)
  original <- fx_var_fitted()[["model"]]

  for (field in c("type", "k", "p", "n", "endogen", "error", "varsel")) {
    expect_equal(restored[["model"]][[field]], original[[field]],
                 info = field)
  }
})

test_that("the estimation data survives the round trip", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path)
  restored <- read_model_from_hdf5(path)
  train <- fx_var_fitted()[["data"]][["train"]]

  expect_equal(unclass(restored[["data"]][["train"]][["y"]]),
               unclass(train[["y"]]), ignore_attr = TRUE)
  expect_equal(unclass(restored[["data"]][["train"]][["x"]]),
               unclass(train[["x"]]), ignore_attr = TRUE)
  expect_equal(stats::tsp(restored[["data"]][["train"]][["y"]]),
               stats::tsp(train[["y"]]))
})

test_that("a restored model can be used downstream", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path)
  restored <- read_model_from_hdf5(path)

  expect_equal(irf(restored, impulse = "Dp", response = "r",
                   n_ahead = 3),
               irf(fx_var_fitted(), impulse = "Dp", response = "r",
                   n_ahead = 3))
})

test_that("a VEC model survives a write and read round trip", {
  path <- temp_h5_file()
  write_to_hdf5(fx_vec_fitted(), filename = path)
  restored <- read_model_from_hdf5(path)

  expect_s3_class(restored, "bvecmodel")
  expect_equal(unclass(restored[["posterior"]][["beta"]][["coeffs"]]),
               unclass(fx_vec_fitted()[["posterior"]][["beta"]][["coeffs"]]),
               ignore_attr = TRUE)
  expect_equal(restored[["model"]][["rank"]],
               fx_vec_fitted()[["model"]][["rank"]])
})

test_that("a modellist is written to a folder and read back", {
  folder <- temp_model_dir()
  write_to_hdf5(fx_var_modellist(), folder = folder)

  expect_length(list.files(folder, pattern = "[.]h5$"),
                length(fx_var_modellist()))

  restored <- read_models_from_folder(folder)
  expect_s3_class(restored, "modellist")
  expect_length(restored, length(fx_var_modellist()))
  expect_equal(
    unname(sort(vapply(restored, function(x) x[["model"]][["p"]], integer(1)))),
    unname(sort(vapply(fx_var_modellist(), function(x) x[["model"]][["p"]],
                       integer(1))))
  )
  # The models come back named after the files they were read from, which the
  # list written here has no counterpart for.
  expect_false(is.null(names(restored)))
})

test_that("an expanding window is written into its own subfolder", {
  folder <- temp_model_dir()
  write_to_hdf5(fx_expanding_window(), folder = folder)

  # The writer names a subfolder after the model specification and puts one
  # file per window into it.
  subfolders <- list.dirs(folder, recursive = FALSE)
  expect_length(subfolders, 1)
  expect_match(basename(subfolders), "ExpWind")
  expect_length(list.files(subfolders), length(fx_expanding_window()))

  # The individual files are ordinary model files.
  restored <- read_model_from_hdf5(list.files(subfolders, full.names = TRUE)[1])
  expect_s3_class(restored, "bvarmodel")
})

test_that("an expanding window can be read back from its folder", {
  folder <- temp_model_dir()
  write_to_hdf5(fx_expanding_window(), folder = folder)
  restored <- read_expanding_window_model_from_folder(
    list.dirs(folder, recursive = FALSE)
  )

  expect_s3_class(restored, "expandingwindow")
  expect_length(restored, length(fx_expanding_window()))
  expect_true(all(vapply(restored, inherits, logical(1), "bvarmodel")))
})

# A write that cannot finish used to be wrapped in try(): the error was
# discarded, the function returned as though it had worked, and a half-written
# file was left where a complete one was expected. The three tests below are
# about the failure path, which the round trips above never reach.

test_that("a write that cannot finish raises the error", {
  path <- temp_h5_file()

  # An element hdf5r has no way to store. Everything before it writes, and then
  # the write fails part way through.
  broken <- fx_var_fitted()
  broken[["data"]][["train"]][["y"]] <- function() NULL

  expect_error(write_to_hdf5(broken, filename = path))
})

test_that("a write that cannot finish leaves no file behind", {
  path <- temp_h5_file()

  broken <- fx_var_fitted()
  broken[["data"]][["train"]][["y"]] <- function() NULL
  try(write_to_hdf5(broken, filename = path), silent = TRUE)

  # Nothing on disk, so the "already exists" guard does not turn the next
  # attempt into a second, misleading error.
  expect_false(file.exists(path))
  expect_no_error(write_to_hdf5(fx_var_fitted(), filename = path))
  expect_true(file.exists(path))
})

test_that("a failed write does not leave the file open", {
  path <- temp_h5_file()

  broken <- fx_var_fitted()
  broken[["data"]][["train"]][["y"]] <- function() NULL
  try(write_to_hdf5(broken, filename = path), silent = TRUE)

  # An HDF5 handle left open would keep the file locked for the rest of the
  # session, so a fresh write to the same path has to work.
  expect_no_error(write_to_hdf5(fx_var_fitted(), filename = path))
  expect_s3_class(read_model_from_hdf5(path), "bvarmodel")
})

# What is still open in the file at `path`, counted through a second handle on
# it, which is not counted itself. HDF5 keeps a file open for as long as
# anything in it is, even once the handle of the file is closed, and on Windows
# the file is then locked against every other process -- BayesTS among them.
# Counting straight after the write is what makes this independent of when the
# garbage collector runs: it would release what was left open, but later.
open_in_h5_file <- function(path) {
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close(), add = TRUE)
  h5$get_obj_count() - 1
}

test_that("a VAR write leaves nothing open in its file", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path)

  expect_equal(open_in_h5_file(path), 0)
})

test_that("a VEC write leaves nothing open in its file", {
  path <- temp_h5_file()
  write_to_hdf5(fx_vec_fitted(), filename = path)

  expect_equal(open_in_h5_file(path), 0)
})

test_that("a write of sign restrictions leaves nothing open in its file", {
  # The settings of the restrictions are attributes of a group of their own.
  path <- temp_h5_file()
  write_to_hdf5(fx_var_sign(), filename = path)

  expect_equal(open_in_h5_file(path), 0)
})

test_that("a write into a group leaves nothing open in its file", {
  path <- temp_h5_file()
  write_to_hdf5(fx_vec_fitted(), filename = path, group = "/models/vec")

  expect_equal(open_in_h5_file(path), 0)
})

test_that("a successful write returns the path invisibly", {
  path <- temp_h5_file()

  expect_invisible(result <- write_to_hdf5(fx_var_fitted(), filename = path))
  expect_equal(result, path)
})

test_that("a VEC write that cannot finish behaves the same way", {
  path <- temp_h5_file()

  broken <- fx_vec_fitted()
  broken[["data"]][["train"]][["y"]] <- function() NULL

  expect_error(write_to_hdf5(broken, filename = path))
  expect_false(file.exists(path))
  expect_no_error(write_to_hdf5(fx_vec_fitted(), filename = path))
})

test_that("an existing file is refused and left untouched", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path)
  before <- file.info(path)[["size"]]

  expect_error(write_to_hdf5(fx_var_fitted(), filename = path), "already exists")
  # Refused before the handle was opened, so the first file is still whole.
  expect_equal(file.info(path)[["size"]], before)
  expect_s3_class(read_model_from_hdf5(path), "bvarmodel")
})

# Groups: one file holding several models, addressed the way the BayesTS
# command line addresses them.

test_that("group names are normalised the way BayesTS spells them", {
  expect_equal(bvartools:::.normalize_hdf5_group(""), "")
  expect_equal(bvartools:::.normalize_hdf5_group("/"), "")
  expect_equal(bvartools:::.normalize_hdf5_group(NULL), "")
  expect_equal(bvartools:::.normalize_hdf5_group("/models/3"), "/models/3")
  expect_equal(bvartools:::.normalize_hdf5_group("models/3"), "/models/3")
  expect_equal(bvartools:::.normalize_hdf5_group("/models/3/"), "/models/3")
  # Idempotent, so a name that has been through it can go through again.
  expect_equal(bvartools:::.normalize_hdf5_group(
    bvartools:::.normalize_hdf5_group("models/3/")), "/models/3")

  for (bad in c("/models//3", "/models/./3", "/models/../3", "//")) {
    expect_error(bvartools:::.normalize_hdf5_group(bad), info = bad)
  }
  expect_error(bvartools:::.normalize_hdf5_group(c("a", "b")))
})

test_that("a model survives a round trip through a group", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/models/3")

  restored <- read_model_from_hdf5(path, group = "/models/3")

  expect_s3_class(restored, "bvarmodel")
  expect_equal(unclass(restored[["posterior"]][["a"]][["coeffs"]]),
               unclass(fx_var_fitted()[["posterior"]][["a"]][["coeffs"]]),
               ignore_attr = TRUE)
  expect_equal(restored[["model"]][["p"]], fx_var_fitted()[["model"]][["p"]])
})

test_that("a model written under a group is at that group and nowhere else", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/models/3")

  h5 <- hdf5r::h5file(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  # The raw paths, not only what comes back out: a prefix that was silently
  # dropped would write to the root and read back the same numbers.
  expect_true(h5$exists("/models/3/model"))
  expect_true(h5$exists("/models/3/data/train/y"))
  expect_false(h5$exists("/model"))
  expect_false(h5$exists("/data"))
  # The intermediate group is created on the way.
  expect_true(h5$exists("/models"))
})

test_that("several models live side by side in one file", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/submodels/US")
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/submodels/JP")
  write_to_hdf5(fx_vec_fitted(), filename = path, group = "/submodels/CA")

  expect_equal(list_models_in_hdf5(path),
               c("/submodels/CA", "/submodels/JP", "/submodels/US"))

  # Each is read back as the model it is, VEC included.
  expect_s3_class(read_model_from_hdf5(path, "/submodels/US"), "bvarmodel")
  expect_s3_class(read_model_from_hdf5(path, "/submodels/CA"), "bvecmodel")
})

test_that("writing into a group that is taken is refused", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/models/3")

  expect_error(write_to_hdf5(fx_var_fitted(), filename = path, group = "/models/3"),
               "already exists")
  # An existing file is still refused when no group is given: without one the
  # model is the whole file.
  expect_error(write_to_hdf5(fx_var_fitted(), filename = path), "already exists")
})

test_that("a failed write into a group leaves the models beside it alone", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/submodels/US")

  broken <- fx_var_fitted()
  broken[["data"]][["train"]][["y"]] <- function() NULL
  expect_error(write_to_hdf5(broken, filename = path, group = "/submodels/JP"))

  # The file is still there, the model that was already in it is intact, and
  # the half-written group is gone so the name can be used again.
  expect_true(file.exists(path))
  expect_equal(list_models_in_hdf5(path), "/submodels/US")
  expect_s3_class(read_model_from_hdf5(path, "/submodels/US"), "bvarmodel")
  expect_no_error(write_to_hdf5(fx_var_fitted(), filename = path,
                                group = "/submodels/JP"))
})

test_that("reading a group that is not there is an error", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/models/3")

  expect_error(read_model_from_hdf5(path, group = "/models/9"),
               "does not contain group")
  expect_error(list_models_in_hdf5(path, group = "/models/9"),
               "does not contain group")
})

test_that("list_models_in_hdf5 finds models and nothing else", {
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/submodels/US")

  # A model at the root of its own file is reported as the root.
  root_path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = root_path)
  expect_equal(list_models_in_hdf5(root_path), "")

  # A group that is not a model is not reported, and the walk does not descend
  # into a model's own subtree.
  h5 <- hdf5r::h5file(path, mode = "a")
  h5$create_group("global")
  h5$close_all()

  found <- list_models_in_hdf5(path)
  expect_equal(found, "/submodels/US")
  expect_false(any(grepl("/data|/priors|/posterior", found)))

  # A root that restricts the walk.
  expect_equal(list_models_in_hdf5(path, group = "/submodels"), "/submodels/US")
  expect_equal(list_models_in_hdf5(path, group = "/global"), character(0))
})

test_that("read_models_from_folder names every model it returns", {
  folder <- temp_model_dir()
  dir.create(file.path(folder, "US"))
  dir.create(file.path(folder, "JP"))
  write_to_hdf5(fx_var_fitted(), filename = file.path(folder, "US", "001.h5"))
  write_to_hdf5(fx_var_fitted(), filename = file.path(folder, "JP", "001.h5"))

  restored <- read_models_from_folder(folder)

  # Flat and named, so a caller can tell which sub-model is which. The nesting
  # this used to return dropped the names entirely.
  expect_s3_class(restored, "modellist")
  expect_length(restored, 2)
  expect_equal(names(restored), c("JP/001", "US/001"))
  expect_true(all(vapply(restored, inherits, logical(1), "bvarmodel")))
})

test_that("read_models_from_folder reads every model of a multi-model file", {
  folder <- temp_model_dir()
  path <- file.path(folder, "models.h5")
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/submodels/US")
  write_to_hdf5(fx_var_fitted(), filename = path, group = "/submodels/JP")

  restored <- read_models_from_folder(folder)

  expect_length(restored, 2)
  expect_equal(names(restored),
               c("models:/submodels/JP", "models:/submodels/US"))
})

test_that("the shape of the result does not depend on how deep the caller points", {
  folder <- temp_model_dir()
  dir.create(file.path(folder, "US"))
  write_to_hdf5(fx_var_fitted(), filename = file.path(folder, "US", "001.h5"))

  outer <- read_models_from_folder(folder)
  inner <- read_models_from_folder(file.path(folder, "US"))

  # One flat list either way; only the names differ, because the models sit at
  # different depths below what was asked for.
  expect_length(outer, 1)
  expect_length(inner, 1)
  expect_equal(names(outer), "US/001")
  expect_equal(names(inner), "001")
  expect_equal(outer[[1]][["model"]], inner[[1]][["model"]])
})

test_that("an expanding window is recognised from the files, not the path", {
  folder <- temp_model_dir()
  write_to_hdf5(fx_expanding_window(), folder = folder)

  restored <- read_models_from_folder(folder)
  expect_s3_class(restored, "expandingwindow")

  # The models say so themselves, so a directory renamed away from "ExpWind"
  # is still read as one.
  written <- list.dirs(folder, recursive = FALSE)
  renamed <- file.path(folder, "plain-name")
  file.rename(written, renamed)

  expect_s3_class(read_models_from_folder(folder), "expandingwindow")
})

test_that("a file without the class attribute is classed by its algorithm", {
  # The writer stores the class of the object, so the algorithm behind this
  # fallback is only consulted for a file written before it did. Every
  # algorithm has to be named there: one that is not comes back as a bare list
  # that no method of the package applies to.
  path <- temp_h5_file()
  write_to_hdf5(fx_var_fitted(), filename = path)

  handle <- hdf5r::H5File$new(path, mode = "r+")
  handle[["model"]]$attr_delete("rclass")
  handle$close_all()

  restored <- read_model_from_hdf5(path)
  expect_null(restored[["model"]][["rclass"]])
  expect_s3_class(restored, "bvarmodel")
})

# A quantile VAR. Local to this file: the fixtures of test-quantile_var.R are
# not shared between test files.
ald_fitted_h5 <- function(tvp = FALSE) {
  cached_fixture(paste0("ald_fitted_hdf5_", tvp), {
    object <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                               error = "ald", quantile = 0.25, tvp = tvp,
                               iterations = fx_iterations, burnin = fx_burnin)
    object <- add_priors(object,
                         coef = if (tvp) list(v_i = 1, shape = 3, rate = 1e-8)
                                else list(v_i = 1),
                         sigma = list(shape = 3, rate = 0.01))
    set.seed(987654)
    add_posterior_coefficients(add_initial_values(object))
  })
}

test_that("a quantile VAR survives a write and read round trip", {
  object <- ald_fitted_h5()
  expect_identical(object[["model"]][["algorithm"]], "VarNormalAld")

  path <- temp_h5_file()
  write_to_hdf5(object, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_s3_class(restored, "bvarmodel")
  expect_identical(restored[["model"]][["algorithm"]], "VarNormalAld")
  # The estimand itself. A model that loses its quantile estimates the median
  # instead and looks entirely healthy, which is why it is checked by value.
  expect_equal(restored[["model"]][["quantile"]], 0.25)
  expect_equal(unclass(restored[["posterior"]][["a"]][["coeffs"]]),
               unclass(object[["posterior"]][["a"]][["coeffs"]]),
               ignore_attr = TRUE)
})

test_that("the scale of a quantile VAR survives the round trip", {
  object <- ald_fitted_h5()

  path <- temp_h5_file()
  write_to_hdf5(object, filename = path)
  restored <- read_model_from_hdf5(path)

  # An ald error keeps its shape and rate under u_scale rather than u_sigma,
  # because the scale it describes is not a prior on Sigma. Both the prior and
  # the draws of the scale used to be dropped: the writer stopped on the error
  # specification before reaching them.
  expect_equal(as.numeric(restored[["priors"]][["u_scale"]][["shape"]]),
               as.numeric(object[["priors"]][["u_scale"]][["shape"]]))
  expect_equal(as.numeric(restored[["priors"]][["u_scale"]][["rate"]]),
               as.numeric(object[["priors"]][["u_scale"]][["rate"]]))
  expect_equal(unclass(restored[["posterior"]][["u_scale"]][["coeffs"]]),
               unclass(object[["posterior"]][["u_scale"]][["coeffs"]]),
               ignore_attr = TRUE)
})

test_that("a time varying quantile VAR is written as well", {
  object <- ald_fitted_h5(tvp = TRUE)
  expect_identical(object[["model"]][["algorithm"]], "VarTvpAld")

  path <- temp_h5_file()
  write_to_hdf5(object, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_s3_class(restored, "bvarmodel")
  expect_identical(restored[["model"]][["algorithm"]], "VarTvpAld")
  expect_equal(restored[["model"]][["quantile"]], 0.25)
  expect_equal(unclass(restored[["posterior"]][["u_scale"]][["coeffs"]]),
               unclass(object[["posterior"]][["u_scale"]][["coeffs"]]),
               ignore_attr = TRUE)
})

test_that("a quantile VAR without the class attribute is classed by algorithm", {
  # VarNormalAld and VarTvpAld were missing from the list the fallback reads,
  # so a quantile VAR written before the class attribute existed lost its class
  # and came back as a bare list that no method applies to.
  for (tvp in c(FALSE, TRUE)) {
    path <- temp_h5_file()
    write_to_hdf5(ald_fitted_h5(tvp = tvp), filename = path)

    handle <- hdf5r::H5File$new(path, mode = "r+")
    handle[["model"]]$attr_delete("rclass")
    handle$close_all()

    expect_s3_class(read_model_from_hdf5(path), "bvarmodel")
  }
})

# A time varying model, whose coefficient draws come with the variance of the
# state equation beside them. Local to this file.
tvp_fitted_h5 <- function() {
  cached_fixture("tvp_fitted_hdf5", {
    object <- create_bvarmodel(var_data(), p = 1, deterministic = "const",
                               tvp = TRUE, error = "wishart",
                               iterations = fx_iterations, burnin = fx_burnin)
    object <- add_priors(object, coef = list(v_i = 1, shape = 3, rate = 1e-8),
                         sigma = list(df = 3, scale = 1))
    set.seed(987654)
    add_posterior_coefficients(add_initial_values(object))
  })
}

test_that("the state variance of a time varying model survives the round trip", {
  object <- tvp_fitted_h5()
  # The thing that is easy to lose: a second set of draws in the same group as
  # the coefficients.
  expect_true("sigma" %in% names(object[["posterior"]][["a"]]))

  path <- temp_h5_file()
  write_to_hdf5(object, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_identical(names(restored[["posterior"]][["a"]]),
                   names(object[["posterior"]][["a"]]))
  expect_equal(unclass(restored[["posterior"]][["a"]][["sigma"]]),
               unclass(object[["posterior"]][["a"]][["sigma"]]),
               ignore_attr = TRUE)
  # Read with the start, end and thinning interval it was written with, so that
  # window() and thin() still mean the same thing on the restored draws.
  expect_identical(coda::mcpar(restored[["posterior"]][["a"]][["sigma"]]),
                   coda::mcpar(object[["posterior"]][["a"]][["sigma"]]))
})

test_that("draws kept outside a group are still read", {
  object <- add_posterior_loglik(tvp_fitted_h5())

  path <- temp_h5_file()
  write_to_hdf5(object, filename = path)
  restored <- read_model_from_hdf5(path)

  # loglik is a dataset of its own rather than a group of draws, and is told
  # apart by that rather than by its name.
  expect_equal(unclass(restored[["posterior"]][["loglik"]]),
               unclass(object[["posterior"]][["loglik"]]),
               ignore_attr = TRUE)
})

# A fitted VAR carrying forecasts and the errors they made against the periods
# held back from the estimation sample. Local to this file.
forecast_errors_fitted_h5 <- function() {
  cached_fixture("forecast_errors_hdf5", {
    data <- var_data()
    train <- stats::window(data, end = c(1997, 1))
    test <- stats::window(data, start = c(1997, 2))

    object <- create_bvarmodel(train, p = 1, deterministic = "const",
                               iterations = fx_iterations, burnin = fx_burnin)
    object <- add_priors(object, coef = list(v_i = 1),
                         sigma = list(df = 3, scale = 1))
    set.seed(987654)
    object <- add_posterior_coefficients(add_initial_values(object))
    object <- add_forecast_input(object, n_ahead = 4)
    object <- add_posterior_forecasts(object)
    add_forecast_errors(object, test_sample = test)
  })
}

test_that("forecast errors survive the round trip", {
  object <- forecast_errors_fitted_h5()
  # The name the rest of the package uses. The writers looked for the singular
  # of it, which nothing produces, so the errors were dropped on the way out
  # and the loss was visible only on reading the file back. They are a member of
  # the forecast group now, which is the name the file has to carry.
  expect_true("errors" %in% names(object[["posterior"]][["forecast"]]))

  path <- temp_h5_file()
  write_to_hdf5(object, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_true("errors" %in% names(restored[["posterior"]][["forecast"]]))
  expect_equal(unclass(restored[["posterior"]][["forecast"]][["errors"]]),
               unclass(object[["posterior"]][["forecast"]][["errors"]]),
               ignore_attr = TRUE)
  # The forecasts they were computed from come back as well.
  expect_equal(unclass(restored[["posterior"]][["forecast"]][["forecasts"]]),
               unclass(object[["posterior"]][["forecast"]][["forecasts"]]),
               ignore_attr = TRUE)
})

test_that("the values a model was scored against survive the round trip", {
  object <- forecast_errors_fitted_h5()
  expect_false(is.null(object[["data"]][["test"]][["y"]]))

  path <- temp_h5_file()
  write_to_hdf5(object, filename = path)

  # /data/test/y, beside /data/train/y: what a forecast is scored against is
  # input, like everything else under /data, and the errors taken against it are
  # draws and live in the posterior.
  h5 <- hdf5r::H5File$new(path, mode = "r")
  expect_true("test" %in% names(h5[["data"]]))
  expect_identical(h5[["data/test/y"]]$dims,
                   c(4L, as.integer(object[["model"]][["k"]])))
  h5$close_all()

  restored <- read_model_from_hdf5(path)
  expect_equal(restored[["data"]][["test"]][["y"]], object[["data"]][["test"]][["y"]],
               ignore_attr = TRUE)
})

test_that("a restored model scores itself without a test sample", {
  path <- temp_h5_file()
  write_to_hdf5(forecast_errors_fitted_h5(), filename = path)
  restored <- read_model_from_hdf5(path)

  # Which is what the values are in the file for: a window of an expanding
  # window exercise carries what it is to be judged by.
  again <- add_forecast_errors(restored)
  expect_equal(unclass(get_forecast_errors(again)),
               unclass(get_forecast_errors(forecast_errors_fitted_h5())),
               ignore_attr = TRUE)
})

test_that("a restored model still reports its forecast errors", {
  path <- temp_h5_file()
  write_to_hdf5(forecast_errors_fitted_h5(), filename = path)
  restored <- read_model_from_hdf5(path)

  # Which is what they are stored for.
  expect_equal(get_forecast_errors(restored),
               get_forecast_errors(forecast_errors_fitted_h5()))
})

test_that("a round trip returns the elements that were written, as they were", {
  # The members of an HDF5 group are listed in the order of their names, not of
  # their creation, so a list is compared element by element under its names.
  by_name <- function(x) {
    if (is.list(x) && !inherits(x, c("mcmc", "ts", "data.frame")) && !is.null(names(x))) {
      cls <- class(x)
      x <- lapply(x[order(names(x))], by_name)
      class(x) <- cls
    }
    x
  }

  # Scalar and character priors, the classes of the series, the class of the
  # model and elements holding NULL all used to come back different.
  for (model in list(fx_at_vec_tvp(), fx_at_var())) {
    path <- temp_h5_file()
    write_to_hdf5(model, filename = path)
    expect_equal(by_name(read_model_from_hdf5(path)), by_name(model))
  }
})

test_that("a file without shape marks is read as it always was", {
  path <- temp_h5_file()
  write_to_hdf5(fx_at_vec_tvp(), filename = path)

  # A file written by BayesTS or from Python carries no marks, and its values
  # come back as matrices.
  h5 <- hdf5r::h5file(path, mode = "r+")
  dataset <- h5[["priors/beta/rho"]]
  dataset$attr_delete("rshape")
  dataset$close()
  h5$close_all()

  expect_true(is.matrix(read_model_from_hdf5(path)[["priors"]][["beta"]][["rho"]]))
})

test_that("the starting precision of a constant gamma model is where BayesTS reads it", {
  # Regression test. BayesTS reads the starting error precision of
  # VarNormalGamma and VecNormalGamma from /initial/u_sigma_inv, while R keeps
  # it as u_omega_inv. The file used to carry the R name, so BayesTS refused a
  # structural VAR with "initial error precision must be 3x3, got 0x0".
  data("e1")
  e1 <- diff(log(e1)) * 100
  model <- create_bvarmodel(e1, p = 1, deterministic = "const",
                            structural = TRUE, error = "gamma",
                            iterations = 10, burnin = 5)
  model <- add_priors(model, coef = list(v_i = 0.1, v_i_det = 0.1),
                      sigma = list(shape = 3, rate = 1e-4))
  model <- add_initial_values(model)

  path <- temp_h5_file()
  write_to_hdf5(model, filename = path)

  h5 <- hdf5r::H5File$new(path, mode = "r")
  in_file <- names(h5[["initial"]])
  h5$close_all()
  expect_true("u_sigma_inv" %in% in_file)
  expect_false("u_omega_inv" %in% in_file)

  restored <- read_model_from_hdf5(path)
  expect_equal(restored[["initial"]][["u_omega_inv"]],
               model[["initial"]][["u_omega_inv"]], ignore_attr = TRUE)
  expect_null(restored[["initial"]][["u_sigma_inv"]])
})

# --- Reading part of a chain -----------------------------------------------------
#
# A caller that works through a long chain in pieces reads the draws it is
# working on rather than the whole posterior.

test_that("read_model_from_hdf5 reads the draws it is asked for", {
  model <- fx_var_model()
  set.seed(4711)
  model <- add_posterior_loglik(add_posterior_coefficients(add_initial_values(
    add_priors(model, coef = list(v_i = 1), sigma = list(df = 3, scale = 1e-8)))))

  file <- file.path(tempdir(), "bvartools-partial-read.h5")
  unlink(file)
  write_to_hdf5(model, filename = file)

  whole <- read_model_from_hdf5(filename = file)
  part <- read_model_from_hdf5(filename = file, draws = c(2L, 5L))

  a <- unclass(whole[["posterior"]][["a"]][["coeffs"]])
  expect_equal(unclass(part[["posterior"]][["a"]][["coeffs"]]),
               a[c(2, 5), , drop = FALSE], ignore_attr = TRUE)
  expect_equal(nrow(part[["posterior"]][["u_sigma_inv"]][["coeffs"]]), 2L)
  expect_equal(unclass(part[["posterior"]][["loglik"]]),
               unclass(whole[["posterior"]][["loglik"]])[c(2, 5), , drop = FALSE],
               ignore_attr = TRUE)

  # Everything that is not a draw comes back as it does from a full read.
  expect_identical(part[["model"]], whole[["model"]])
  expect_equal(part[["data"]], whole[["data"]])
  expect_equal(part[["priors"]], whole[["priors"]])
  expect_s3_class(part, class(whole)[1])
})

test_that("no draws gives the model without its chain", {
  model <- fx_var_model()
  set.seed(4712)
  model <- add_posterior_coefficients(add_initial_values(
    add_priors(model, coef = list(v_i = 1), sigma = list(df = 3, scale = 1e-8))))

  file <- file.path(tempdir(), "bvartools-no-draws.h5")
  unlink(file)
  write_to_hdf5(model, filename = file)

  none <- read_model_from_hdf5(filename = file, draws = integer(0))
  expect_equal(nrow(none[["posterior"]][["a"]][["coeffs"]]), 0L)
  expect_equal(ncol(none[["posterior"]][["a"]][["coeffs"]]),
               ncol(model[["posterior"]][["a"]][["coeffs"]]))
  expect_equal(none[["data"]], read_model_from_hdf5(filename = file)[["data"]])
})

test_that("draws that are not in the chain are refused", {
  model <- fx_var_model()
  set.seed(4713)
  model <- add_posterior_coefficients(add_initial_values(
    add_priors(model, coef = list(v_i = 1), sigma = list(df = 3, scale = 1e-8))))

  file <- file.path(tempdir(), "bvartools-bad-draws.h5")
  unlink(file)
  write_to_hdf5(model, filename = file)

  draws <- nrow(model[["posterior"]][["a"]][["coeffs"]])
  expect_error(read_model_from_hdf5(filename = file, draws = draws + 1L),
               "the chain holds")
  expect_error(read_model_from_hdf5(filename = file, draws = 0L), "at least one")
  expect_error(read_model_from_hdf5(filename = file, draws = 1.5), "whole numbers")
})


test_that("the writer's explicit types are the ones hdf5r would guess", {

  # write_to_hdf5() hands hdf5r a cached type and dataspace instead of letting it
  # guess them, which is most of what a write costs. The file must not change
  # for it, so each cached type is checked against the guess hdf5r makes with
  # the string length it passes, and each dataspace against guess_space().
  values <- list(character = c("VecNormalWishart", "a"), scalar_character = "none",
                 logical = TRUE, integer = 3L, integers = 1:4, double = 0.98,
                 doubles = c(1990, 2023.75, 4), matrix = matrix(1:6 * 1.0, 2))
  for (name in names(values)) {
    value <- values[[name]]
    guessed <- hdf5r::guess_dtype(value, scalar = FALSE, string_len = Inf)
    expect_equal(.hdf5_dtype(value)$to_text(), guessed$to_text(), info = name)
    expected <- hdf5r::guess_space(value, dtype = guessed, chunked = FALSE)
    expect_equal(.hdf5_attr_space(value)$get_simple_extent_dims(),
                 expected$get_simple_extent_dims(), info = name)
  }

  # What it has no type for is left to hdf5r.
  expect_null(.hdf5_dtype(factor("a")))
  expect_null(.hdf5_dtype(list(1)))
  expect_null(.hdf5_attr_space(list(1)))
})
