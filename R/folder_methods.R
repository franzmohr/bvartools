#' Estimation Steps on a Folder of Models
#'
#' The steps of an estimation, applied to the models of a folder one at a time
#' rather than to a list of models in the session. Each reads a model, changes
#' it, writes it back and drops it, so a step costs one model per worker.
#'
#' @param object an object of class 'bvarfolder', from
#' \code{\link{open_models}}.
#' @param x an object of class 'bvarfolder'.
#' @param seed an integer, the seed of the first model. The models are numbered
#' through in the order of the manifest, as the models of a list are.
#' @param cores the number of worker processes. Defaults to one.
#' @param export character vector of names of objects the workers need, for
#' instance the sampler that \code{posterior_function} refers to. See
#' \code{\link{map_models}}.
#' @param ... further arguments passed to the method for a single model.
#'
#' @details
#' The methods do what their counterparts for a 'bvarmodel' or a 'bvecmodel'
#' do, model by model, and they take the same arguments.
#'
#' Each model is drawn with a seed of its own, so the draws do not depend on how
#' many workers a step runs on, nor on whether the models were estimated in one
#' run or in several. The seeds are the ones a list of models gets from
#' \code{\link{add_seed}}: \code{seed}, \code{seed + 1}, and so on in the order
#' of the manifest.
#'
#' A step takes \code{models}, which \code{\link{map_models}} understands: a
#' character vector of the models it is applied to, the rest of the folder being
#' left as it is. An estimation of many expensive models is run in parts that
#' way, and one that was interrupted is picked up where it stopped.
#'
#' \code{add_posterior_coefficients} is the step to think twice about. It writes
#' the largest thing a model has, and a cluster of R workers carries every draw
#' into the session and out again to do it. Where the sampler is the BayesTS
#' executable, \code{\link{bayests_files}} runs it on the files instead:
#' \preformatted{run <- bayests_files(executable = "/opt/bayests/bin/bayests")
#' run(model_files(stored), jobs = 6)}
#' which gives the same draws, since a model is drawn with the seed in its file.
#'
#' @return The object in \code{object}, with its manifest brought up to date.
#'
#' @examples
#'
#' data("e1")
#' train <- diff(log(e1)) * 100
#'
#' models <- create_bvarmodel(data = train, p = 1:2, deterministic = "const",
#'                            iterations = 20, burnin = 10)
#'
#' folder <- file.path(tempdir(), "bvartools-example-steps")
#' unlink(folder, recursive = TRUE)
#' dir.create(folder, recursive = TRUE)
#' write_to_hdf5(models, folder = folder)
#'
#' stored <- open_models(folder)
#' stored <- add_priors(stored, coef = list(v_i = 1),
#'                      sigma = list(df = 3, scale = 1))
#' stored <- add_initial_values(stored)
#' stored <- add_seed(stored, 20260919)
#' stored <- add_posterior_coefficients(stored)
#' stored <- add_posterior_loglik(stored)
#'
#' @name folder_steps
#' @family model comparison
NULL

#' @rdname folder_steps
#' @export
#' @method add_priors bvarfolder
add_priors.bvarfolder <- function(object, ..., cores = 1) {
  map_models(object, function(model, ...) add_priors(model, ...), ...,
             cores = cores)
}

#' @rdname folder_steps
#' @export
#' @method add_initial_values bvarfolder
add_initial_values.bvarfolder <- function(object, ..., cores = 1) {
  map_models(object, function(model, ...) add_initial_values(model, ...), ...,
             cores = cores)
}

#' @rdname folder_steps
#' @export
#' @method add_seed bvarfolder
add_seed.bvarfolder <- function(object, seed, ..., cores = 1) {
  # The seeds are worked out here rather than on the workers: what a worker
  # runs is sent to it as it stands, so it must not call anything the package
  # installed there may not have.
  seeds <- .offset_seed(.check_seed(seed),
                        seq_len(nrow(object[["manifest"]])) - 1L)
  map_models(object, function(model, index, seeds, ...) {
    bvartools::add_seed(model, seeds[[index]], ...)
  }, seeds = seeds, ..., cores = cores)
}

#' @rdname folder_steps
#' @export
#' @method add_posterior_coefficients bvarfolder
add_posterior_coefficients.bvarfolder <- function(object, ..., cores = 1,
                                                  export = NULL) {
  map_models(object, function(model, ...) add_posterior_coefficients(model, ...),
             ..., cores = cores, export = export)
}

#' @rdname folder_steps
#' @export
#' @method add_posterior_loglik bvarfolder
add_posterior_loglik.bvarfolder <- function(object, ..., cores = 1) {
  map_models(object, function(model, ...) add_posterior_loglik(model, ...), ...,
             cores = cores)
}

#' @rdname folder_steps
#' @export
#' @method thin bvarfolder
thin.bvarfolder <- function(x, ..., cores = 1) {
  map_models(x, function(model, ...) thin(model, ...), ..., cores = cores)
}

#' @rdname folder_steps
#' @details
#' \code{selection_criteria} is the one method that reads rather than writes:
#' it hands back a criterion per model and leaves the folder as it is. It is
#' what \code{\link{choose_best_model}} compares, so a grid of specifications
#' too large to hold is still chosen from.
#'
#' @export
#' @method selection_criteria bvarfolder
selection_criteria.bvarfolder <- function(object, ..., cores = 1) {

  result <- map_models(object, function(model, ...) selection_criteria(model, ...),
                       ..., write = FALSE, cores = cores)
  class(result) <- append("selcritlist", class(result))
  result
}
