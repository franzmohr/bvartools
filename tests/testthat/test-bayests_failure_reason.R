test_that("the failing file is reported, not the tail of the log", {
  # What a directory walk looks like when the third of five files fails and the
  # rest are already drawn: the real reason is nowhere near the end.
  log <- c(
    "Processing: models/AT_r1/2021.00.h5",
    "Posterior data already exists in file. Skipping simulation.",
    "Processing: models/AT_r1/2021.25.h5",
    "Posterior data already exists in file. Skipping simulation.",
    "Processing: models/AT_r1/2021.50.h5",
    "Error processing models/AT_r1/2021.50.h5: Failed to open HDF5 file: File has been truncated",
    "Processing: models/AT_r1/2021.75.h5",
    "Posterior data already exists in file. Skipping simulation.",
    "Processing: models/AT_r1/2022.00.h5",
    "Posterior data already exists in file. Skipping simulation.")

  reason <- bvartools:::.bayests_failure_reason(log)

  expect_length(reason, 1L)
  expect_match(reason, "2021.50.h5")
  expect_match(reason, "File has been truncated")
  expect_false(any(grepl("already exists", reason)))
})

test_that("every failing file is named, up to the limit", {
  log <- c("Processing: a.h5", "Error processing a.h5: one",
           "Processing: b.h5", "Error processing b.h5: two",
           "Processing: c.h5", "Posterior data already exists in file. Skipping simulation.")

  reason <- bvartools:::.bayests_failure_reason(log)
  expect_length(reason, 2L)
  expect_match(reason[1], "a.h5")
  expect_match(reason[2], "b.h5")

  many <- paste0("Error processing ", 1:9, ".h5: no")
  expect_length(bvartools:::.bayests_failure_reason(many), 5L)
  expect_length(bvartools:::.bayests_failure_reason(many, n = 2), 2L)
})

test_that("without a marked error the tail is still reported", {
  # A process killed from outside writes no "Error processing" line, and the
  # last thing it managed to say is the best guess there is.
  log <- c("Processing: a.h5", "Progress: 50%", "Progress: 60%")
  expect_identical(bvartools:::.bayests_failure_reason(log), log)
  expect_identical(bvartools:::.bayests_failure_reason(character()), character())
  expect_identical(bvartools:::.bayests_failure_reason(as.character(1:8)),
                   as.character(4:8))
})
