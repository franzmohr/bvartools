# The seven FRED-MD/FRED-QD transformation codes. Each is checked against the
# arithmetic it names, and each keeps the time index of its input: the periods a
# transformation reaches back over come back as NA rather than being dropped.

test_that("every code computes what it names and keeps the index", {

  x <- stats::ts(cbind(a = c(1, 2, 4, 7, 11, 16), b = c(2, 3, 5, 8, 13, 21)),
                 start = c(2000, 1), frequency = 4)
  a <- as.numeric(x[, "a"])

  expected <- list(
    `1` = a,
    `2` = c(NA, diff(a)),
    `3` = c(NA, NA, diff(a, differences = 2)),
    `4` = log(a),
    `5` = c(NA, diff(log(a))),
    `6` = c(NA, NA, diff(log(a), differences = 2)),
    # The difference of a growth rate reaches back two periods, not one.
    `7` = c(NA, NA, diff(a[-1] / a[-length(a)] - 1)))

  for (code in names(expected)) {
    out <- transform_variables(x, c(a = as.integer(code)))
    expect_equal(stats::tsp(out), stats::tsp(x))
    expect_equal(as.numeric(out[, "a"]), expected[[code]])
    # A column without a code is left alone.
    expect_equal(as.numeric(out[, "b"]), as.numeric(x[, "b"]))
  }

  expect_equal(sum(is.na(transform_variables(x, c(a = 7L))[, "a"])), 2L)
})


test_that("a single series keeps its shape and takes an unnamed code", {

  x <- stats::ts(c(1, 2, 4, 7, 11), start = 2000)

  out <- transform_variables(x, 2)
  expect_null(dim(out))
  expect_equal(as.numeric(out), c(NA, 1, 2, 3, 4))
  expect_equal(stats::tsp(out), stats::tsp(x))

  # A named code on a series without a name used to reach code[[NA]].
  expect_equal(as.numeric(transform_variables(x, c(y = 2))), c(NA, 1, 2, 3, 4))

  named <- stats::ts(matrix(c(1, 2, 4, 7, 11), dimnames = list(NULL, "y")), start = 2000)
  expect_equal(as.numeric(transform_variables(named, c(y = 2))), c(NA, 1, 2, 3, 4))
  expect_error(transform_variables(named, c(z = 2)), "not found")
})


test_that("invalid input is refused with a reason", {

  x <- stats::ts(cbind(a = 1:4, b = 2:5))
  expect_error(transform_variables(1:4, 1), "class 'ts'")
  expect_error(transform_variables(x, c(a = 8)), "1 to 7")
  expect_error(transform_variables(x, c(a = 1.5)), "1 to 7")
  expect_error(transform_variables(x, 2), "named after the columns")
  expect_error(transform_variables(x, c(a = 2, a = 3)), "duplicate")
  expect_error(transform_variables(stats::ts(1:4), c(2, 3)), "single code")
})
