# The documentation in inst/agents/ is what a coding agent copies into someone
# else's analysis, so an example there that no longer matches the package fails
# in their session rather than here -- unless it is run here. Every ```r block of
# every file under the skill is evaluated, top to bottom in one environment per
# file, and the stopifnot() calls in the examples assert the shapes the text
# states. Blocks in any other language (```text for shell commands) are not run.
#
# Skipped on CRAN: the examples run real samplers, if briefly.

skip_on_cran()

agent_doc_blocks <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  starts <- which(lines == "```r")
  closes <- which(lines == "```")
  lapply(starts, function(start) {
    end <- closes[closes > start][1]
    lines[seq.int(start + 1, end - 1)]
  })
}

skill_dir <- system.file("agents", "skills", "bvartools", package = "bvartools")
skill_files <- list.files(skill_dir, pattern = "[.]md$", recursive = TRUE, full.names = TRUE)

test_that("the agent documentation is installed with the package", {
  expect_true(nzchar(skill_dir))
  expect_true(file.exists(file.path(skill_dir, "SKILL.md")))
})

for (path in skill_files) {
  blocks <- agent_doc_blocks(path)
  if (length(blocks) == 0) {
    next
  }

  test_that(paste("the R examples in", basename(path), "run"), {
    env <- new.env(parent = globalenv())
    for (i in seq_along(blocks)) {
      failure <- tryCatch({
        suppressWarnings(eval(parse(text = blocks[[i]]), envir = env))
        NULL
      }, error = function(e) conditionMessage(e))
      expect(is.null(failure),
             sprintf("%s: R example %d of %d failed: %s",
                     basename(path), i, length(blocks), failure))
    }
  })
}
