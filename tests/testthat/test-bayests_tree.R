# write_bayests_tree(), read_bayests_tree() and from_bayests_tree(): the layer
# beneath write_to_hdf5() that another package's models are written and read
# through, and a 'modellist' holding such a model beside a VAR. The VAR and VEC
# writers themselves are test-hdf5.R's.

skip_if_not_installed("hdf5r")

# A model class of another package, as dfmtools' factor models are: written as
# a tree in BayesTS's names, and read back by undoing the renaming. Registered
# rather than defined here, because the generics dispatch from the namespace.
local({
  registerS3method("write_to_hdf5", "testmodel", function(object, filename, group = "", ...) {
    write_bayests_tree(list(
      "model" = list(".attributes" = c(object[["model"]],
                                       list("rclass" = class(object)))),
      "data" = list("train" = list("y" = object[["data"]][["x"]]))
    ), filename = filename, group = group)
  }, envir = asNamespace("bvartools"))
  registerS3method("from_bayests_tree", "testmodel", function(tree, ...) {
    model <- tree[["model"]][[".attributes"]]
    model[["rclass"]] <- NULL
    structure(list("data" = list("x" = tree[["data"]][["train"]][["y"]]),
                   "model" = model),
              class = class(tree))
  }, envir = asNamespace("bvartools"))
})

test_model <- function(n = 2L) {
  x <- stats::ts(matrix(seq_len(20) / 10, 10, 2), start = c(2000, 1), frequency = 4)
  colnames(x) <- c("a", "b")
  structure(list("data" = list("x" = x),
                 "model" = list("algorithm" = "TestFactor", "k" = 2L,
                                "n" = n, "p" = 1L)),
            class = c("testmodel", "list"))
}

test_that("a tree survives a write and read round trip as what it was", {
  path <- temp_h5_file()
  y <- stats::ts(matrix(rnorm(20), 10, 2), start = c(1990, 2), frequency = 4)
  colnames(y) <- c("gdp", "cpi")
  draws <- coda::mcmc(matrix(rnorm(12), 4, 3), start = 11, thin = 2)
  tree <- list(
    "model" = list(".attributes" = list("algorithm" = "Example", "k" = 2L)),
    "data" = list("train" = list("y" = y)),
    "priors" = list("lambda" = list("mu" = c(0, 0, 0), "v_inv" = diag(3)),
                    "unused" = NULL),
    "posterior" = list("lambda" = draws)
  )

  expect_identical(write_bayests_tree(tree, filename = path), path)
  back <- read_bayests_tree(path)

  expect_identical(back[["model"]][[".attributes"]][["algorithm"]], "Example")
  expect_identical(back[["model"]][[".attributes"]][["k"]], 2L)
  expect_equal(back[["data"]][["train"]][["y"]], y)
  expect_identical(back[["priors"]][["lambda"]][["mu"]], c(0, 0, 0))
  expect_identical(back[["priors"]][["lambda"]][["v_inv"]], diag(3))
  expect_false("unused" %in% names(back[["priors"]]))
  expect_equal(unclass(back[["posterior"]][["lambda"]]), unclass(draws),
               ignore_attr = TRUE)
  expect_identical(coda::mcpar(back[["posterior"]][["lambda"]]), coda::mcpar(draws))
})

test_that("a partial read reads only the draws asked for, and only of the posterior", {
  path <- temp_h5_file()
  draws <- coda::mcmc(matrix(as.numeric(1:12), 4, 3))
  write_bayests_tree(list("model" = list(".attributes" = list("algorithm" = "Example")),
                          "initial" = list("a" = matrix(1, 4, 3)),
                          "posterior" = list("a" = list("coeffs" = draws))),
                     filename = path)

  back <- read_bayests_tree(path, draws = c(2, 4))

  expect_equal(unclass(back[["posterior"]][["a"]][["coeffs"]]),
               unclass(draws)[c(2, 4), ], ignore_attr = TRUE)
  expect_identical(dim(back[["initial"]][["a"]]), c(4L, 3L))
})

test_that("the writer refuses an existing model and keeps the models beside a new one", {
  path <- temp_h5_file()
  tree <- list("model" = list(".attributes" = list("algorithm" = "Example")))
  write_bayests_tree(tree, filename = path, group = "/one")

  expect_error(write_bayests_tree(tree, filename = path, group = "/one"),
               "already exists")
  write_bayests_tree(tree, filename = path, group = "two")
  expect_identical(list_models_in_hdf5(path), c("/one", "/two"))

  expect_error(write_bayests_tree(tree, filename = path), "already exists")
})

test_that("a write that fails leaves nothing behind", {
  path <- temp_h5_file()
  bad <- list("model" = list(".attributes" = list("algorithm" = "Example")),
              "data" = list("y" = function() NULL))
  expect_error(write_bayests_tree(bad, filename = path))
  expect_false(file.exists(path))

  good <- list("model" = list(".attributes" = list("algorithm" = "Example")))
  write_bayests_tree(good, filename = path, group = "/kept")
  expect_error(write_bayests_tree(bad, filename = path, group = "/failed"))
  expect_identical(list_models_in_hdf5(path), "/kept")
})

test_that("a tree must be named all the way down", {
  path <- temp_h5_file()
  expect_error(write_bayests_tree(list(1, 2), filename = path), "named list")
  expect_error(write_bayests_tree(list("data" = list(1)), filename = path),
               "no name")
  expect_error(write_bayests_tree(list("model" = list(".attributes" = list(1))),
                                  filename = path),
               "'.attributes' of '/model'")
})

test_that("read_model_from_hdf5() hands a model of another package to its method", {
  path <- temp_h5_file()
  model <- test_model()
  write_to_hdf5(model, filename = path)

  back <- read_model_from_hdf5(path)

  expect_s3_class(back, "testmodel")
  expect_equal(back[["data"]][["x"]], model[["data"]][["x"]])
  expect_identical(back[["model"]][["algorithm"]], "TestFactor")
  expect_identical(back[["model"]][["n"]], 2L)
})

test_that("a class without a from_bayests_tree() method is refused by name", {
  path <- temp_h5_file()
  write_bayests_tree(list("model" = list(".attributes" = list(
    "algorithm" = "Example", "rclass" = c("orphanmodel", "list")))),
    filename = path)

  expect_error(read_model_from_hdf5(path), "class 'orphanmodel'")
  expect_type(read_bayests_tree(path), "list")
})

test_that("a modellist holding another package's model beside a VAR goes through a folder", {
  folder <- temp_model_dir()
  models <- structure(list(fx_var_fitted(), test_model(n = 1L), test_model(n = 3L)),
                      class = c("modellist", "list"))

  paths <- write_to_hdf5(models, folder = folder)
  expect_length(paths, 3)
  # The factor-like model has no 's' and its 'm' is not an exogenous count,
  # so its file name says only what it has.
  expect_false(any(grepl("-s=-", basename(paths), fixed = TRUE)))

  back <- read_models_from_folder(folder)

  expect_s3_class(back, "modellist")
  expect_length(back, 3)
  expect_s3_class(back[[1]], "bvarmodel")
  expect_s3_class(back[[2]], "testmodel")
  expect_s3_class(back[[3]], "testmodel")
  expect_identical(back[[2]][["model"]][["n"]], 1L)
  expect_identical(back[[3]][["model"]][["n"]], 3L)
  expect_null(back[[2]][["model"]][["rindex_collection"]])
})

test_that("the package a file names is loaded, and a missing one is named", {
  path <- temp_h5_file()
  write_bayests_tree(list("model" = list(".attributes" = list(
    "algorithm" = "Example", "rclass" = c("orphanmodel", "list"),
    "rpackage" = "notapackageanywhere"))),
    filename = path)

  expect_error(read_model_from_hdf5(path), "Install notapackageanywhere")
})
