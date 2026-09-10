# The announcements that make this a transition release: every function that
# bvartools 1.0.0 does not have any more says so once per session, names its
# successor where there is one, and returns exactly what it always returned.

# Forget which functions have already been announced, so that each expectation
# starts from the state of a fresh session.
reset_transition <- function() {
  registry <- bvartools:::.bvartools_transition
  rm(list = ls(envir = registry, all.names = TRUE), envir = registry)
  # The helper switches the announcements off for the rest of the suite; here
  # they are the subject, so they are switched back on.
  options(bvartools.transition.messages = TRUE)
}

# The announcement of every function that is removed or renamed in 1.0.0,
# together with the successor the message has to name. NA where there is none.
transition_functions <- list(
  gen_var = "create_bvarmodel",
  gen_vec = "create_bvecmodel",
  bvec_to_bvar = "vec_to_var",
  draw_posterior = "add_posterior_coefficients",
  bvarpost = "add_posterior_coefficients",
  bvecpost = "add_posterior_coefficients",
  kalman_dk = "kalman_durbin_koopman_2002",
  stochvol_ksc1998 = "stochvol_ksc_1998",
  stochvol_ocsn2007 = "stochvol_ocsn_2007",
  stoch_vol = "stochvol_ksc_1998",
  bvs = "post_bvs",
  post_normal_covar_const = NA_character_,
  post_normal_covar_tvp = NA_character_,
  dfm = NA_character_,
  dfmpost = NA_character_,
  gen_dfm = NA_character_
)

test_that("every removed function announces itself and names its successor", {
  for (fn in names(transition_functions)) {
    reset_transition()
    msg <- capture_messages(bvartools:::.transition_message(
      fn, successor = if (is.na(transition_functions[[fn]])) NULL else transition_functions[[fn]]))

    expect_length(msg, 1)
    expect_match(msg, paste0("'", fn, "()' will be removed in bvartools 1.0.0"), fixed = TRUE)

    if (is.na(transition_functions[[fn]])) {
      expect_match(msg, "It has no replacement.", fixed = TRUE)
    } else {
      expect_match(msg, paste0("Use '", transition_functions[[fn]], "()' instead."), fixed = TRUE)
    }
  }
})

test_that("an announcement is made once per session", {
  reset_transition()

  expect_message(gen_var(bvartools_test_data(), p = 1, deterministic = "const",
                         iterations = 10, burnin = 5),
                 "gen_var")
  expect_silent(suppressWarnings(gen_var(bvartools_test_data(), p = 1, deterministic = "const",
                                         iterations = 10, burnin = 5)))
})

test_that("announcements can be switched off", {
  reset_transition()
  options(bvartools.transition.messages = FALSE)
  expect_silent(gen_var(bvartools_test_data(), p = 1, deterministic = "const",
                        iterations = 10, burnin = 5))
})

test_that("one call to draw_posterior gives one message", {
  reset_transition()

  object <- suppressMessages(gen_var(bvartools_test_data(), p = 1, deterministic = "const",
                                     iterations = 20, burnin = 10))
  object <- add_priors(object)

  msg <- capture_messages(draw_posterior(object))
  # 'bvarpost' is an internal step of 'draw_posterior' and points at the same
  # successor, so it must not produce a second announcement.
  expect_false(any(grepl("bvarpost", msg, fixed = TRUE)))
  expect_true(any(grepl("draw_posterior", msg, fixed = TRUE)))
})

test_that("the wrappers of the compiled functions return what they always returned", {
  reset_transition()
  options(bvartools.transition.messages = FALSE)

  y <- bvartools_test_data()
  h_init <- log(diag(stats::var(y)))
  h <- t(matrix(h_init, ncol(y), nrow(y)))
  sigma <- rep(0.05, ncol(y))
  constant <- rep(0.0001, ncol(y))

  # The wrapper must not consume any random numbers of its own: from one seed
  # the draw has to be the draw of the compiled routine.
  set.seed(99)
  through_wrapper <- stochvol_ksc1998(y, h, sigma, h_init, constant)
  set.seed(99)
  direct <- bvartools:::.stochvol_ksc1998_cpp(y, h, sigma, h_init, constant)

  expect_identical(through_wrapper, direct)
  expect_equal(dim(through_wrapper), dim(h))
})

test_that("the compiled samplers do not announce anything of their own", {
  reset_transition()

  # 'bvaralg' reaches the stochastic volatility draw through the registered C++
  # symbol rather than through the R binding, so estimating a model must not
  # produce the announcement of 'stochvol_ksc1998'.
  object <- suppressMessages(gen_var(bvartools_test_data(), p = 1, deterministic = "const",
                                     sv = TRUE, iterations = 20, burnin = 10))
  object <- add_priors(object, sigma = list(mu = 0, v_i = 0.01, shape = 3, rate = 0.0001,
                                            sigma_h = 0.05, constant = 0.0001))

  msg <- capture_messages(draw_posterior(object))
  expect_false(any(grepl("stochvol", msg, fixed = TRUE)))
  expect_false(any(grepl("stoch_vol", msg, fixed = TRUE)))
})

# Leave the suite as the helper set it up, so that the order of the files cannot
# change what the other tests see.
options(bvartools.transition.messages = FALSE)
