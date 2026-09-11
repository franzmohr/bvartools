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

  expect_equal(irf(restored, impulse = "income", response = "cons",
                   n_ahead = 3),
               irf(fx_var_fitted(), impulse = "income", response = "cons",
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
    train <- stats::window(data, end = c(1977, 4))
    test <- stats::window(data, start = c(1978, 1))

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
  # and the loss was visible only on reading the file back.
  expect_true("forecast_errors" %in% names(object[["posterior"]]))

  path <- temp_h5_file()
  write_to_hdf5(object, filename = path)
  restored <- read_model_from_hdf5(path)

  expect_true("forecast_errors" %in% names(restored[["posterior"]]))
  expect_equal(unclass(restored[["posterior"]][["forecast_errors"]]),
               unclass(object[["posterior"]][["forecast_errors"]]),
               ignore_attr = TRUE)
  # The forecasts they were computed from come back as well.
  expect_equal(unclass(restored[["posterior"]][["forecast"]]),
               unclass(object[["posterior"]][["forecast"]]),
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
