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
    sort(vapply(restored, function(x) x[["model"]][["p"]], integer(1))),
    sort(vapply(fx_var_modellist(), function(x) x[["model"]][["p"]],
                integer(1)))
  )
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
