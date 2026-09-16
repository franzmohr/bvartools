# Simulation of the models of a list on several worker processes.
#
# The workers load the installed bvartools, so these tests only say something
# about the code under test when that is what is installed -- under R CMD check,
# not under devtools::load_all().
skip_if_workers_load_other_version <- function() {
  skip_on_cran()
  skip_if(isNamespaceLoaded("pkgload") && pkgload::is_dev_package("bvartools"),
          "workers would load the installed bvartools, not the one under test")
}

# Two lag orders over three expanding windows: a 'modellist' of
# 'expandingwindow' lists with six models.
parallel_fixture <- function() {
  cached_fixture("parallel_models", {
    models <- create_bvarmodel(var_data(), p = 1:2, deterministic = "const",
                               iterations = 20, burnin = 5)
    models <- add_priors(models, coef = list(v_i = 0, v_i_det = 0),
                         sigma = list(df = 1, scale = 0.0001))
    models <- use_expanding_window(models, start = c(1997, 4))
    set.seed(31)
    add_initial_values(models)
  })
}

draws_of <- function(object, element) {
  lapply(bvartools:::.model_paths(object),
         function(path) object[[path]][["posterior"]][[element]])
}

test_that("'cores' must be a whole number of at least one", {
  models <- parallel_fixture()
  for (cores in list(0, 1.5, NA_real_, "2", c(1, 2))) {
    expect_error(add_posterior_coefficients(models, cores = cores), "'cores'")
  }
})

test_that("models are found through nested lists and external forecasts are kept whole", {
  models <- parallel_fixture()
  expect_s3_class(models, "modellist")
  expect_s3_class(models[[1]], "expandingwindow")
  expect_identical(bvartools:::.model_paths(models),
                   list(c(1L, 1L), c(1L, 2L), c(1L, 3L), c(2L, 1L), c(2L, 2L), c(2L, 3L)))

  external <- structure(list(list(), list()),
                        class = c("externalforecast", "expandingwindow", "list"))
  nested <- structure(list(models[[1]], external), class = c("modellist", "list"))
  expect_identical(bvartools:::.model_paths(nested),
                   list(c(1L, 1L), c(1L, 2L), c(1L, 3L), 2L))
})

# Evaluates 'fun(object, ...)' one model after the other, in a process that runs
# its BLAS on one thread like the workers do. An optimised BLAS rounds some
# results differently on one thread than on several, and a chain carries that
# forward, so this session, which may run on several, is no reference.
on_one_thread <- function(fun, object, ...) {
  cl <- bvartools:::.start_model_cluster(1)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, fun, object, ...)[[1]]
}

test_that("coefficient draws on two workers equal those of one process", {
  skip_if_workers_load_other_version()
  models <- parallel_fixture()

  sequential <- on_one_thread(add_posterior_coefficients, models)
  parallel <- add_posterior_coefficients(models, cores = 2)

  expect_s3_class(parallel, "modellist")
  expect_s3_class(parallel[[2]], "expandingwindow")
  expect_identical(draws_of(parallel, "a"), draws_of(sequential, "a"))
  expect_identical(draws_of(parallel, "u_sigma_inv"), draws_of(sequential, "u_sigma_inv"))

  # They do not depend on the number of workers either.
  expect_identical(draws_of(add_posterior_coefficients(models, cores = 3), "a"),
                   draws_of(parallel, "a"))

  # An expanding window on its own is shared out as well.
  window <- add_posterior_coefficients(models[[1]], cores = 2)
  expect_identical(draws_of(window, "a"), draws_of(sequential[[1]], "a"))
})

test_that("log-likelihoods on two workers equal those of one process", {
  skip_if_workers_load_other_version()
  models <- add_posterior_coefficients(parallel_fixture())

  expect_identical(draws_of(add_posterior_loglik(models, cores = 2), "loglik"),
                   draws_of(on_one_thread(add_posterior_loglik, models), "loglik"))
})

test_that("forecasts on two workers are reproduced by set.seed()", {
  skip_if_workers_load_other_version()
  models <- add_forecast_input(add_posterior_coefficients(parallel_fixture()),
                               n_ahead = 2)

  set.seed(5)
  first <- add_posterior_forecasts(models, cores = 2)
  set.seed(5)
  second <- add_posterior_forecasts(models, cores = 2)

  expect_identical(draws_of(first, "forecast"), draws_of(second, "forecast"))
  expect_identical(lapply(draws_of(first, "forecast"), dim),
                   lapply(draws_of(add_posterior_forecasts(models), "forecast"), dim))
})

test_that("workers run with one BLAS thread and the session keeps its own", {
  skip_if_workers_load_other_version()
  report_threads <- function(object) {
    object[["threads"]] <- Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                                        "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))
    object
  }
  environment(report_threads) <- globalenv()

  old <- Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"), unset = NA)
  on.exit({
    for (v in names(old)) {
      if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(old[v]))
    }
  }, add = TRUE)
  Sys.setenv(OPENBLAS_NUM_THREADS = "4")
  Sys.unsetenv("OMP_NUM_THREADS")

  result <- add_posterior_coefficients(parallel_fixture(),
                                       posterior_function = report_threads, cores = 2)

  for (path in bvartools:::.model_paths(result)) {
    expect_identical(unname(result[[path]][["threads"]]), rep("1", 4))
  }
  expect_identical(Sys.getenv("OPENBLAS_NUM_THREADS"), "4")
  expect_identical(Sys.getenv("OMP_NUM_THREADS", unset = NA), NA_character_)
})

test_that("workers load bvartools from the library of this session", {
  skip_if_workers_load_other_version()
  report_library <- function(object) {
    object[["library"]] <- getNamespaceInfo("bvartools", "path")
    object
  }
  environment(report_library) <- globalenv()

  result <- add_posterior_coefficients(parallel_fixture(),
                                       posterior_function = report_library, cores = 2)

  for (path in bvartools:::.model_paths(result)) {
    expect_identical(normalizePath(result[[path]][["library"]]),
                     normalizePath(getNamespaceInfo("bvartools", "path")))
  }
})

test_that("a model without a seed is given one before it is sent off", {
  skip_if_workers_load_other_version()
  models <- parallel_fixture()
  models[[1]][[2]][["model"]][["seed"]] <- NULL

  set.seed(8)
  first <- add_posterior_coefficients(models, cores = 2)
  set.seed(8)
  second <- add_posterior_coefficients(models, cores = 2)

  expect_true(is.integer(first[[1]][[2]][["model"]][["seed"]]))
  expect_identical(draws_of(first, "a"), draws_of(second, "a"))
})

test_that("an error on a worker is raised as it would be on one core", {
  skip_if_workers_load_other_version()
  models <- parallel_fixture()
  models[[2]][[3]][["model"]][["algorithm"]] <- "NoSuchSampler"

  message_of <- function(expr) tryCatch(expr, error = conditionMessage)
  message <- message_of(add_posterior_coefficients(models, cores = 2))
  expect_type(message, "character")
  expect_identical(message, message_of(add_posterior_coefficients(models)))
})
