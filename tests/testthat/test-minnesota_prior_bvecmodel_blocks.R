# The Minnesota prior of a VEC model with exogenous variables and several
# deterministic terms, recalculated block by block from the formulas in
# ?minnesota_prior.bvecmodel.
#
# The prior of the deterministic terms overwrote the last block of the
# exogenous variables, and of the deterministic terms only the first entered the
# regressions that give the residual standard deviations.

test_that("every block of a VEC model's Minnesota prior follows its formula", {
  data("us_macrodata", package = "bvartools")
  kappa <- list(kappa1 = 2, kappa2 = 0.5, kappa3 = 0.5, kappa4 = 5)

  for (s in 1:2) {
    model <- create_bvecmodel(data = us_macrodata[, 1:2],
                              exogen = us_macrodata[, 3, drop = FALSE],
                              p = 2, s = s, r = 1,
                              const = "unrestricted", trend = "restricted",
                              seasonal = "unrestricted",
                              iterations = 10, burnin = 5)
    prior <- minnesota_prior(model, kappa1 = kappa$kappa1, kappa2 = kappa$kappa2,
                             kappa3 = kappa$kappa3, kappa4 = kappa$kappa4)

    y <- unclass(model[["data"]][["train"]][["y"]])
    w <- unclass(model[["data"]][["train"]][["w"]])
    x <- unclass(model[["data"]][["train"]][["x"]])
    k <- 2
    r <- 1
    vars <- c("Dp", "u")
    det_x <- c("const", "season.1", "season.2", "season.3")

    # Residual standard deviations from regressions of each variable on its own
    # lagged level, its own lagged difference and every deterministic term.
    s_endo <- vapply(1:k, function(i) {
      X <- cbind(w[, paste0("l.", vars[i])], x[, paste0("d.", vars[i], ".l01")],
                 w[, "trend"], x[, det_x])
      e <- qr.resid(qr(X), y[, i])
      sqrt(sum(e^2) / (nrow(X) - ncol(X)))
    }, numeric(1))

    v <- 1 / diag(prior[["v_i"]])
    expect_equal(v[1:(k * r)], rep(kappa$kappa1, k * r), info = s)
    V <- matrix(v[-(1:(k * r))], k, dimnames = list(NULL, colnames(x)))

    s_x <- stats::sd(x[, "d.r.l00"])
    for (l in 1:k) {
      for (lag in 0:(s - 1)) {
        expect_equal(V[l, sprintf("d.r.l%02d", lag)],
                     kappa$kappa1 * kappa$kappa3 / (lag + 1)^2 * s_endo[l]^2 / s_x^2,
                     ignore_attr = TRUE, info = paste("s =", s, "lag", lag))
      }
      for (term in det_x) {
        expect_equal(V[l, term], kappa$kappa1 * kappa$kappa4 * s_endo[l]^2,
                     ignore_attr = TRUE, info = paste("s =", s, term))
      }
    }
  }
})
