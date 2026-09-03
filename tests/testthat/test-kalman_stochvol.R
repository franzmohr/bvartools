test_that("the simulation smoother returns one state per period plus one", {
  set.seed(1)
  k <- 1
  tt <- 20
  m <- 2

  z <- matrix(stats::rnorm(k * tt * m), k * tt, m)
  y <- matrix(stats::rnorm(k * tt), k, tt)

  states <- kalman_durbin_koopman_2002(y, z, diag(1, k), diag(0.01, m),
                                       diag(1, m), matrix(0, m), diag(1, m))

  expect_equal(dim(states), c(m, tt + 1))
  expect_true(all(is.finite(states)))
})

test_that("a vanishing state variance keeps the state path flat", {
  set.seed(2)
  k <- 1
  tt <- 30
  m <- 1

  z <- matrix(1, k * tt, m)
  y <- matrix(stats::rnorm(k * tt), k, tt)

  states <- kalman_durbin_koopman_2002(y, z, diag(1, k), diag(1e-12, m),
                                       diag(1, m), matrix(0, m), diag(1e-12, m))

  # Without state innovations the random walk cannot move.
  expect_lt(max(abs(diff(as.numeric(states)))), 1e-4)
})

test_that("the smoother is reproducible from the seed", {
  k <- 1
  tt <- 15
  m <- 2
  z <- matrix(stats::rnorm(k * tt * m), k * tt, m)
  y <- matrix(stats::rnorm(k * tt), k, tt)
  draw <- function() {
    kalman_durbin_koopman_2002(y, z, diag(1, k), diag(0.01, m), diag(1, m),
                               matrix(0, m), diag(1, m))
  }

  set.seed(3)
  first <- draw()
  set.seed(3)
  second <- draw()
  set.seed(4)
  third <- draw()

  expect_identical(first, second)
  expect_false(isTRUE(all.equal(first, third)))
})

test_that("the smoother tracks a state that moves", {
  set.seed(5)
  tt <- 300
  m <- 1

  # A slow random walk observed with little noise.
  true_state <- cumsum(stats::rnorm(tt, 0, 0.1))
  z <- matrix(1, tt, m)
  y <- matrix(true_state + stats::rnorm(tt, 0, 0.05), 1, tt)

  states <- kalman_durbin_koopman_2002(y, z, diag(0.05^2, 1), diag(0.1^2, m),
                                       diag(1, m), matrix(0, m), diag(1, m))

  expect_gt(stats::cor(as.numeric(states)[seq_len(tt)], true_state), 0.9)
})

for (sampler in c("stochvol_ksc_1998", "stochvol_ocsn_2007")) {

  test_that(paste(sampler, "returns a volatility path of the right shape"), {
    set.seed(6)
    tt <- 40
    y <- matrix(stats::rnorm(tt))
    h <- matrix(0, tt, 1)

    draw <- do.call(sampler,
                    list(y, h, matrix(0.05), matrix(0), matrix(0.0001)))

    expect_equal(dim(draw), c(tt, 1))
    expect_true(all(is.finite(draw)))
  })

  test_that(paste(sampler, "keeps the path near its start without innovations"), {
    set.seed(7)
    tt <- 40
    y <- matrix(stats::rnorm(tt) * exp(0.5 * 1.5))
    h <- matrix(1.5, tt, 1)

    # A state variance this small leaves the log volatility no room to move.
    draw <- do.call(sampler,
                    list(y, h, matrix(1e-10), matrix(1.5), matrix(0.0001)))

    expect_equal(as.numeric(draw), rep(1.5, tt), tolerance = 1e-3)
  })

  test_that(paste(sampler, "picks up a difference in scale"), {
    set.seed(8)
    tt <- 400
    # The second half of the sample is far more volatile than the first.
    y <- matrix(c(stats::rnorm(tt / 2, 0, 0.1),
                  stats::rnorm(tt / 2, 0, 10)))
    h <- matrix(0, tt, 1)

    draw <- do.call(sampler,
                    list(y, h, matrix(0.5), matrix(0), matrix(0.0001)))

    expect_gt(mean(draw[(tt / 2 + 1):tt]), mean(draw[seq_len(tt / 2)]))
  })
}
