# Dynamic multipliers: the response of the endogenous variables to a change in
# a weakly exogenous one.

varx_model <- function(iterations = 30) {
  data("e1", envir = environment())
  e1 <- diff(log(e1)) * 100
  model <- create_bvarmodel(data = e1[, c("invest", "income")],
                            exogen = e1[, "cons", drop = FALSE],
                            p = 2, s = 1, deterministic = "const",
                            iterations = iterations, burnin = 10)
  model <- add_priors(model, coef = list(v_i = 1, v_i_det = 1 / 10),
                      sigma = list(df = "k", scale = 1))
  add_posterior_coefficients(add_initial_values(model))
}

# The recursion of the documentation, written out for one draw.
by_hand <- function(model, draw, impulse, response, n_ahead, permanent = TRUE) {
  k <- model[["model"]][["k"]]
  p <- model[["model"]][["p"]]
  m <- model[["model"]][["m"]]
  s <- model[["model"]][["s"]]
  a <- model[["posterior"]][["a"]][["coeffs"]][draw, ]
  A <- matrix(a[seq_len(k * k * p)], k)
  B <- matrix(a[k * k * p + seq_len(k * m * (s + 1))], k)
  j_imp <- which(model[["model"]][["exogen"]] == impulse)
  path <- matrix(0, n_ahead + 1, k)
  for (h in 0:n_ahead) {
    v <- rep(0, k)
    for (l in seq_len(min(h, p))) {
      v <- v + A[, (l - 1) * k + 1:k] %*% path[h - l + 1, ]
    }
    if (permanent) {
      for (j in 0:min(h, s)) v <- v + B[, j * m + j_imp]
    } else {
      if (h <= s) v <- v + B[, h * m + j_imp]
    }
    path[h + 1, ] <- v
  }
  path[, which(model[["model"]][["endogen"]] == response)]
}

test_that("the multipliers of a VARX model follow the documented recursion", {

  model <- varx_model()
  drawn <- multipliers(model, impulse = "cons", response = "income",
                       n_ahead = 6, keep_draws = TRUE)

  expect_equal(nrow(drawn), nrow(model[["posterior"]][["a"]][["coeffs"]]))
  expect_equal(ncol(drawn), 7)
  for (i in c(1, 5, nrow(drawn))) {
    expect_equal(as.numeric(drawn[i, ]),
                 as.numeric(by_hand(model, i, "cons", "income", 6)), info = i)
  }

  # The impact multiplier is B_0 alone, whatever the lag structure does later.
  k <- model[["model"]][["k"]]
  p <- model[["model"]][["p"]]
  b0 <- model[["posterior"]][["a"]][["coeffs"]][, k * k * p + 1:k]
  expect_equal(as.numeric(drawn[, 1]),
               as.numeric(b0[, which(model[["model"]][["endogen"]] == "income")]))
})

test_that("a transitory change is the permanent one without the earlier lags", {

  model <- varx_model()
  perm <- multipliers(model, impulse = "cons", response = "invest", n_ahead = 5,
                      keep_draws = TRUE)
  tran <- multipliers(model, impulse = "cons", response = "invest", n_ahead = 5,
                      type = "transitory", keep_draws = TRUE)

  # They agree at impact, where only B_0 has come into range.
  expect_equal(as.numeric(perm[, 1]), as.numeric(tran[, 1]))
  expect_false(isTRUE(all.equal(as.numeric(perm[, 3]), as.numeric(tran[, 3]))))
  for (i in c(1, nrow(tran))) {
    expect_equal(as.numeric(tran[i, ]),
                 as.numeric(by_hand(model, i, "cons", "invest", 5, permanent = FALSE)),
                 info = i)
  }
})

test_that("the summary, the shock and the cumulation behave", {

  model <- varx_model()

  summarised <- multipliers(model, impulse = "cons", response = "income",
                            n_ahead = 4, ci = 0.9)
  drawn <- multipliers(model, impulse = "cons", response = "income",
                       n_ahead = 4, keep_draws = TRUE)
  expect_s3_class(summarised, "bvarirf")
  expect_equal(dim(unclass(summarised)), c(5L, 3L))
  expect_equal(as.numeric(unclass(summarised)[, 2]),
               as.numeric(apply(drawn, 2, stats::median)))
  expect_equal(as.numeric(unclass(summarised)[, 1]),
               as.numeric(apply(drawn, 2, stats::quantile, probs = 0.05)))

  # The multipliers are linear in the size of the change.
  doubled <- multipliers(model, impulse = "cons", response = "income",
                         n_ahead = 4, shock = 2, keep_draws = TRUE)
  expect_equal(unclass(doubled), unclass(drawn) * 2)

  cumulated <- multipliers(model, impulse = "cons", response = "income",
                           n_ahead = 4, cumulative = TRUE, keep_draws = TRUE)
  expect_equal(as.numeric(cumulated[1, ]), cumsum(as.numeric(drawn[1, ])))

  expect_equal(ncol(multipliers(model, impulse = "cons", response = "income",
                                n_ahead = 0, keep_draws = TRUE)), 1)
})

test_that("multipliers refuse what they are not defined for", {

  model <- varx_model()
  expect_error(multipliers(model, impulse = "nothing", response = "income"),
               "Impulse variable not available")
  expect_error(multipliers(model, impulse = "cons", response = "nothing"),
               "Response variable not available")
  expect_error(multipliers(model, impulse = "cons", response = "income", n_ahead = -1),
               "at least 0")
  expect_error(multipliers(model, impulse = "cons", response = "income", type = "other"),
               "permanent")

  data("e1", envir = environment())
  e1 <- diff(log(e1)) * 100
  no_exogen <- create_bvarmodel(data = e1, p = 1, deterministic = "const",
                                iterations = 20, burnin = 10)
  expect_error(multipliers(no_exogen, impulse = "cons", response = "income"),
               "weakly exogenous")
})

test_that("a VEC model is multiplied through its levels form", {

  data("e6", envir = environment())
  set.seed(1)
  exogen <- stats::ts(cbind(g = as.numeric(e6[, "R"]) * 0.3 +
                              stats::rnorm(nrow(e6), sd = 0.1)),
                      start = stats::start(e6), frequency = stats::frequency(e6))

  model <- create_bvecmodel(data = e6, exogen = exogen, p = 2, s = 1, r = 1,
                            const = "unrestricted", iterations = 30, burnin = 10)
  model <- add_priors(model, coef = list(v_i = 0),
                      coint = list(v_i = 0, p_tau_i = 1),
                      sigma = list(df = 3, scale = 0.0001))
  model <- add_posterior_coefficients(add_initial_values(model))

  from_vec <- multipliers(model, impulse = "g", response = "R", n_ahead = 6,
                          keep_draws = TRUE)
  from_var <- multipliers(vec_to_var(model), impulse = "g", response = "R",
                          n_ahead = 6, keep_draws = TRUE)
  expect_equal(unclass(from_vec), unclass(from_var))

  # The levels form the multipliers are taken from is a VAR with the exogenous
  # variable the error correction model had.
  as_var <- vec_to_var(model)
  expect_s3_class(as_var, "bvarmodel")
  expect_identical(as_var[["model"]][["exogen"]], "g")

  # And the recursion is the same one the VAR method documents.
  for (i in c(1, nrow(from_vec))) {
    expect_equal(as.numeric(from_vec[i, ]),
                 as.numeric(by_hand(as_var, i, "g", "R", 6)), info = i)
  }
})
