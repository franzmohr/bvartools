# Simulation of the models of a list on several worker processes.
#
# The methods of add_posterior_coefficients(), add_posterior_forecasts() and
# add_posterior_loglik() for 'modellist' and 'expandingwindow' objects take an
# argument 'cores'. With one core they keep their plain lapply over the list.
# With more, the list is flattened into its models -- a 'modellist' built by
# use_expanding_window() holds 'expandingwindow' lists, and handing out whole
# inner lists would leave all but one worker idle when there is only one -- the
# models are simulated on a PSOCK cluster and put back where they came from.

# The BLAS thread counts a worker is started with. An optimised BLAS starts one
# thread per core in every process, so a cluster of n workers would run n times
# as many threads as there are cores, and the samplers, whose sweeps are made
# of small factorisations, gain nothing from threads in the first place. The
# variables are read once, when a process loads its BLAS, which is why they
# are set around the start of the workers and cannot be changed afterwards.
.worker_thread_variables <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                              "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")

.check_cores <- function(cores) {
  if (length(cores) != 1 || !is.numeric(cores) || is.na(cores) ||
      cores < 1 || cores != floor(cores) || cores > .Machine$integer.max) {
    stop("Argument 'cores' must be a single whole number of at least 1.", call. = FALSE)
  }
  as.integer(cores)
}

# Lists whose elements are simulated one by one. An 'externalforecast' is an
# 'expandingwindow' as well, but it has methods of its own that leave it as it
# is, so it is handed on whole like a model.
.is_model_container <- function(x) {
  class(x)[1] %in% c("modellist", "expandingwindow")
}

# The index path of every model in 'x', in the order of its elements, counting
# through nested lists. A path is suitable for x[[path]].
.model_paths <- function(x, path = integer()) {
  if (!.is_model_container(x)) {
    return(list(path))
  }
  unlist(lapply(seq_along(x), function(i) .model_paths(x[[i]], c(path, i))),
         recursive = FALSE)
}

# TRUE if 'object' is to be simulated on a cluster: more than one core asked
# for and more than one model to share out.
.use_cluster <- function(object, cores) {
  .check_cores(cores) > 1 && length(.model_paths(object)) > 1
}

# Runs on a worker. It lives in the namespace so that a cluster call serialises
# it by reference rather than together with the list it was called for. An
# error is returned rather than raised, so that the caller can raise it as the
# sequential methods would.
.simulate_on_worker <- function(model, .fun, ...) {
  tryCatch(.fun(model, ...), error = function(e) e)
}

# Starts 'cores' workers with one BLAS thread each, the library paths of this
# session, and independent random number streams seeded from R's generator.
.start_model_cluster <- function(cores) {
  old <- Sys.getenv(.worker_thread_variables, unset = NA, names = TRUE)
  on.exit({
    was_set <- !is.na(old)
    if (any(was_set)) {
      do.call(Sys.setenv, as.list(old[was_set]))
    }
    if (any(!was_set)) {
      Sys.unsetenv(names(old)[!was_set])
    }
  }, add = TRUE)
  do.call(Sys.setenv, stats::setNames(as.list(rep("1", length(old))), names(old)))

  cl <- parallel::makeCluster(cores)
  ok <- FALSE
  on.exit(if (!ok) parallel::stopCluster(cl), add = TRUE)

  # A worker loads bvartools when it receives the first function of it, from
  # its default library paths unless it is told this session's. .libPaths
  # itself cannot be sent: it keeps the paths in an environment of its own,
  # which would travel as a copy and be set there. A function of the base
  # environment is sent by reference and calls the worker's .libPaths.
  set_library_paths <- function(paths) invisible(.libPaths(paths))
  environment(set_library_paths) <- baseenv()
  parallel::clusterCall(cl, set_library_paths, .libPaths())
  parallel::clusterSetRNGStream(cl, sample.int(.Machine$integer.max, 1L))

  ok <- TRUE
  cl
}

# Applies '.fun' to every model of 'object' on 'cores' workers and returns
# 'object' with each model replaced by the result. With 'seed = TRUE' a model
# without a seed is given one from R's generator first, as its draws would
# otherwise depend on the stream of the worker it happens to land on.
.simulate_models_in_parallel <- function(object, .fun, cores, ..., seed = FALSE) {
  paths <- .model_paths(object)
  models <- lapply(paths, function(path) object[[path]])

  if (seed) {
    for (i in seq_along(models)) {
      if (inherits(models[[i]], c("bvarmodel", "bvecmodel")) &&
          is.null(models[[i]][["model"]][["seed"]])) {
        models[[i]][["model"]][["seed"]] <- .draw_model_seed()
      }
    }
  }

  cl <- .start_model_cluster(min(.check_cores(cores), length(models)))
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # parLapply rather than its load balancing variant: it hands the models out
  # in fixed chunks, so which stream a model draws from depends only on the
  # number of workers, and set.seed() reproduces draws that are not seeded per
  # model, such as forecasts.
  results <- parallel::parLapply(cl, models, .simulate_on_worker, .fun = .fun, ...)
  models <- NULL

  for (i in seq_along(results)) {
    if (inherits(results[[i]], "error")) {
      stop(results[[i]])
    }
  }
  for (i in seq_along(paths)) {
    object[[paths[[i]]]] <- results[[i]]
  }

  object
}
