# Registry of the transition messages that have already been shown in this
# session. The environment lives in the package namespace, so the record is
# reset when the session ends and every message is shown once again in the next
# one. The leading dot keeps it out of the exports, which 'NAMESPACE' otherwise
# hands out to everything whose name begins with a letter.
.bvartools_transition <- new.env(parent = emptyenv())

#' Announce a Function That Will Be Removed in bvartools 1.0.0
#'
#' Emits a message the first time a function is used that does not exist in
#' bvartools 1.0.0 any more. Subsequent calls in the same session are silent.
#'
#' @param fn a character specifying the name of the function that is used.
#' @param successor a character specifying the name of the function that takes
#' its place in bvartools 1.0.0. If \code{NULL} (default), the function is
#' reported as removed without a replacement.
#' @param note a character with an additional sentence that is appended to the
#' message. If \code{NULL} (default), no sentence is appended.
#' @param also a character vector of further function names that are recorded as
#' announced without a message of their own. Used where one function calls
#' another that points at the same successor, so that a single call does not
#' produce a series of messages. If \code{NULL} (default), nothing else is
#' recorded.
#'
#' @details Messages are used instead of warnings so that the announcement can
#' be silenced with \code{\link[base]{suppressMessages}} and does not become an
#' error under \code{options(warn = 2)}. Set
#' \code{options(bvartools.transition.messages = FALSE)} to switch them off.
#'
#' @return \code{NULL} invisibly, called for its side effect.
#'
#' @keywords internal
#' @noRd
.transition_message <- function(fn, successor = NULL, note = NULL, also = NULL) {

  if (isFALSE(getOption("bvartools.transition.messages", TRUE))) {
    return(invisible(NULL))
  }

  # Only announce a function once per session
  if (!is.null(.bvartools_transition[[fn]])) {
    return(invisible(NULL))
  }
  .bvartools_transition[[fn]] <- TRUE

  # Functions that this one calls on its way and that point at the same
  # successor are recorded as announced, so that one call gives one message
  for (i in also) {
    .bvartools_transition[[i]] <- TRUE
  }

  msg <- paste0("bvartools: '", fn, "()' will be removed in bvartools 1.0.0.")
  if (is.null(successor)) {
    msg <- paste0(msg, " It has no replacement.")
  } else {
    msg <- paste0(msg, " Use '", successor, "()' instead.")
  }
  if (!is.null(note)) {
    msg <- paste0(msg, " ", note)
  }

  message(msg,
          "\n  See vignette(\"transition\", package = \"bvartools\") for the full list of changes.",
          "\n  This message is shown once per session. Use",
          " options(bvartools.transition.messages = FALSE) to silence it.")

  invisible(NULL)
}
