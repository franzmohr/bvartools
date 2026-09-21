## Ad hoc check of vec_to_var.bvecmodel(): the VAR representation has to
## reproduce the VEC model's fit, draw for draw, and its data matrices have to
## be the level data the transformed coefficients expect.

library("bvartools", lib.loc = Sys.getenv("BVARTOOLS_LIB", .libPaths()[1]))

ok <- function(what, value) {
  cat(sprintf("  %-52s %s\n", what, if (isTRUE(value)) "ok" else "FAILED"))
  if (!isTRUE(value)) stop("check failed: ", what, call. = FALSE)
}

## Positionally, not by time index: the data matrices are time-series objects
## and '-' would otherwise align them on their tsp attributes, which is exactly
## the thing under test.
close_enough <- function(a, b, tol = 1e-10) {
  a <- matrix(as.numeric(a), NROW(a))
  b <- matrix(as.numeric(b), NROW(b))
  identical(dim(a), dim(b)) && max(abs(a - b)) < tol
}

## Fitted values of a model in SUR form, draw by draw: k x tt per draw.
fitted_draws <- function(object) {
  k <- object[["model"]][["k"]]
  z <- object[["data"]][["train"]][["z"]]
  a <- t(object[["posterior"]][["a"]][["coeffs"]])
  matrix(z %*% a, k)
}

## The same for a VEC model, whose SUR form leaves the columns of the loadings
## to be filled in from the error correction term and the current draw of beta.
fitted_draws_vec <- function(object) {
  k <- object[["model"]][["k"]]
  rank <- object[["model"]][["rank"]]
  k_beta <- object[["model"]][["k_beta"]]
  w <- as.matrix(object[["data"]][["train"]][["w"]])
  z <- object[["data"]][["train"]][["z"]]
  a <- t(object[["posterior"]][["a"]][["coeffs"]])
  beta <- if (rank > 0) t(object[["posterior"]][["beta"]][["coeffs"]]) else NULL

  out <- NULL
  for (draw in 1:ncol(a)) {
    z_draw <- z
    if (rank > 0) {
      z_draw[, 1:(k * rank)] <- kronecker(w %*% matrix(beta[, draw], k_beta), diag(1, k))
    }
    out <- cbind(out, matrix(z_draw %*% a[, draw], k))
  }
  out
}

check_model <- function(label, model) {

  cat("\n", label, "\n", sep = "")

  model <- add_initial_values(model)
  set.seed(123)
  model <- add_posterior_coefficients(model)

  object <- vec_to_var(model)

  ok("class", identical(class(object), c("bvarmodel", "list")))

  k <- model[["model"]][["k"]]
  m <- model[["model"]][["m"]]
  p <- model[["model"]][["p"]]
  s <- model[["model"]][["s"]]
  rank <- model[["model"]][["rank"]]
  n_det <- model[["model"]][["n"]] + ifelse(rank > 0, model[["model"]][["n_restricted"]], 0)

  ok("lag order carried over", object[["model"]][["p"]] == p)
  ok("no cointegration rank left", is.null(object[["model"]][["rank"]]))
  ok("deterministic terms folded together", object[["model"]][["n"]] == n_det)
  ok("same number of observations",
     nrow(object[["data"]][["train"]][["y"]]) == nrow(model[["data"]][["train"]][["y"]]))
  ok("same sample period",
     identical(stats::tsp(object[["data"]][["train"]][["y"]]),
               stats::tsp(model[["data"]][["train"]][["y"]])))

  ## The level data has to be the level data: y_t reconstructed from the
  ## differences must equal the original series over the estimation sample.
  y_original <- stats::window(model[["data"]][["original"]][["endogen"]],
                              start = stats::tsp(object[["data"]][["train"]][["y"]])[1])
  ok("endogenous levels recovered",
     close_enough(object[["data"]][["train"]][["y"]], y_original))

  ## ... and the first lag block must be that series shifted by one period.
  ok("first lag block is the lagged level",
     close_enough(object[["data"]][["train"]][["x"]][, 1:k],
                  stats::window(model[["data"]][["original"]][["endogen"]],
                                start = stats::tsp(object[["data"]][["train"]][["y"]])[1] -
                                  1 / stats::frequency(y_original),
                                end = stats::tsp(object[["data"]][["train"]][["y"]])[2] -
                                  1 / stats::frequency(y_original))))

  if (m > 0) {
    x_original <- stats::window(model[["data"]][["original"]][["exogen"]],
                                start = stats::tsp(object[["data"]][["train"]][["y"]])[1])
    ok("unmodelled variables in levels recovered",
       close_enough(object[["data"]][["train"]][["x"]][, k * p + 1:m], x_original))
  }

  ok("number of coefficients",
     ncol(object[["posterior"]][["a"]][["coeffs"]]) ==
       k * (k * p + m * (s + 1) + n_det))
  ok("chain length preserved",
     nrow(object[["posterior"]][["a"]][["coeffs"]]) ==
       nrow(model[["posterior"]][["u_sigma_inv"]][["coeffs"]]))
  ok("mcpar preserved",
     identical(attr(object[["posterior"]][["a"]][["coeffs"]], "mcpar"),
               attr(model[["posterior"]][["a"]][["coeffs"]], "mcpar")))
  ok("error precision untouched",
     identical(object[["posterior"]][["u_sigma_inv"]], model[["posterior"]][["u_sigma_inv"]]))
  ok("cointegration space dropped", is.null(object[["posterior"]][["beta"]]))

  ## The point of the whole exercise: the two parameterisations describe the
  ## same model, so the fit of the VAR representation is the fit of the VEC
  ## model plus the level it was differenced from.
  fit_vec <- fitted_draws_vec(model) +
    matrix(rep(t(as.matrix(model[["data"]][["train"]][["w"]])[, 1:k]),
               nrow(model[["posterior"]][["a"]][["coeffs"]])), k)
  ok("fitted values agree with the VEC model",
     close_enough(fitted_draws(object), fit_vec, tol = 1e-8))
}

data("e6")
e6 <- e6 * 100

check_model("p = 2, r = 1, restricted constant",
            create_bvecmodel(e6, p = 2, r = 1, const = "restricted",
                             iterations = 20, burnin = 10) |>
              add_priors(coef = list(v_i = 1, v_i_det = 1 / 10),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = "k", scale = 1)))

check_model("p = 1, r = 1, unrestricted constant",
            create_bvecmodel(e6, p = 1, r = 1, const = "unrestricted",
                             iterations = 20, burnin = 10) |>
              add_priors(coef = list(v_i = 1, v_i_det = 1 / 10),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = "k", scale = 1)))

check_model("p = 3, r = 0, unrestricted constant",
            create_bvecmodel(e6, p = 3, r = 0, const = "unrestricted",
                             iterations = 20, burnin = 10) |>
              add_priors(coef = list(v_i = 1, v_i_det = 1 / 10),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = "k", scale = 1)))

check_model("p = 2, r = 1, restricted trend, unrestricted constant",
            create_bvecmodel(e6, p = 2, r = 1, const = "unrestricted", trend = "restricted",
                             iterations = 20, burnin = 10) |>
              add_priors(coef = list(v_i = 1, v_i_det = 1 / 10),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = "k", scale = 1)))

## e6 holds two series only, so the unmodelled, non-deterministic ones come from
## e1: investment and income are modelled, consumption is not.
data("e1")
e1 <- log(e1) * 100
endogen <- e1[, c("invest", "income")]
exogen <- e1[, "cons", drop = FALSE]

check_model("p = 2, s = 2, m = 1, r = 1, restricted constant",
            create_bvecmodel(endogen, p = 2, r = 1, exogen = exogen, s = 2,
                             const = "restricted",
                             iterations = 20, burnin = 10) |>
              add_priors(coef = list(v_i = 1, v_i_det = 1 / 10),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = "k", scale = 1)))

check_model("p = 2, s = 1, m = 1, r = 2, unrestricted constant",
            create_bvecmodel(endogen, p = 2, r = 2, exogen = exogen, s = 1,
                             const = "unrestricted",
                             iterations = 20, burnin = 10) |>
              add_priors(coef = list(v_i = 1, v_i_det = 1 / 10),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = "k", scale = 1)))

check_model("p = 1, s = 1, m = 1, r = 1, restricted trend",
            create_bvecmodel(endogen, p = 1, r = 1, exogen = exogen, s = 1,
                             const = "unrestricted", trend = "restricted",
                             iterations = 20, burnin = 10) |>
              add_priors(coef = list(v_i = 1, v_i_det = 1 / 10),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = "k", scale = 1)))


## The structural case, which no VEC sampler in this package can produce -- a
## structural model needs 'error = "gamma"' or '"sv"' and only VecNormalWishart
## has a VEC implementation. The coefficient transformation is reached directly
## instead, with draws written by hand, so that the one thing that distinguishes
## it is still covered: A_0 stands in the first lag of the VAR representation,
## where a non-structural model has the identity.
cat("\nstructural VEC, coefficient transformation only\n")

k <- 3
rank <- 1
n <- 1
n_structural <- k * (k - 1) / 2
draws <- 5

set.seed(42)
alpha <- matrix(rnorm(k * rank), k)
beta <- matrix(rnorm(k * rank), k)
gamma_1 <- matrix(rnorm(k * k), k)
cc <- matrix(rnorm(k * n), k)
a0_packed <- rnorm(n_structural)

## Unit diagonal, free elements filling the strict lower triangle by column --
## the order create_bvecmodel() leaves the columns of kronecker(-y, I) in.
a_0 <- diag(1, k)
a_0[lower.tri(a_0)] <- a0_packed

object <- list(
  "model" = list("k" = k, "p" = 2, "m" = 0, "s" = 0, "n" = n,
                 "n_restricted" = 0, "rank" = rank, "k_beta" = k,
                 "varsel" = "none", "structural" = TRUE,
                 "iterations" = draws, "burnin" = 0),
  "posterior" = list(
    "a" = list("coeffs" = matrix(rep(c(alpha, cbind(gamma_1, cc), a0_packed), each = draws), draws)),
    "beta" = list("coeffs" = matrix(rep(beta, each = draws), draws)),
    "u_sigma_inv" = list("coeffs" = matrix(rep(diag(1, k), each = draws), draws))))

out <- bvartools:::.VecToVarCoefficients(object)[["a"]][["coeffs"]]

ok("chain length preserved", nrow(out) == draws)
ok("number of coefficients", ncol(out) == k * (k * 2 + n) + n_structural)

level <- matrix(out[1, 1:(k * (k * 2 + n))], k)
ok("A_1 is A_0 + Pi + Gamma_1",
   close_enough(level[, 1:k], a_0 + alpha %*% t(beta) + gamma_1))
ok("A_1 is not I + Pi + Gamma_1",
   !close_enough(level[, 1:k], diag(1, k) + alpha %*% t(beta) + gamma_1))
ok("A_2 is -Gamma_1", close_enough(level[, k + 1:k], -gamma_1))
ok("C carries over", close_enough(level[, 2 * k + 1:n, drop = FALSE], cc))
ok("contemporaneous coefficients carry over",
   close_enough(out[1, ncol(out) - n_structural + 1:n_structural], a0_packed))

cat("\nall checks passed\n")
