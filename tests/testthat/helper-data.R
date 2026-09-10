# A small, fast data set shared by the tests. The West German investment,
# income and consumption series of Luetkepohl (2006) in first log differences,
# which is what the examples and vignettes of this package use throughout.
bvartools_test_data <- function(n = 40) {
  data("e1", package = "bvartools", envir = environment())
  y <- diff(log(get("e1", envir = environment()))) * 100
  stats::ts(y[seq_len(min(n, nrow(y))), , drop = FALSE],
            start = stats::start(y), frequency = stats::frequency(y))
}

# The announcements are noise in tests that are not about them. The tests in
# test-transition-messages.R switch them back on where they are the subject.
options(bvartools.transition.messages = FALSE)
