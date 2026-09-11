# A caller supplied impact matrix, the hook that identification schemes the
# package does not derive itself reach the recursion through.
#
# The recursion is the part that must not have changed, so most of what follows
# feeds `type = "custom"` the matrix the named types build internally and asks
# for the same answer back.

test_that("a custom impact matrix reproduces the forecast error response", {
  model <- fx_var_fitted()
  k <- model[["model"]][["k"]]

  # A single matrix identifies every draw the same way. The identity is what
  # "feir" uses, so the two have to agree.
  expect_equal(
    irf(model, impulse = "income", response = "cons", n_ahead = 4,
        type = "custom", impact = diag(1, k)),
    irf(model, impulse = "income", response = "cons", n_ahead = 4,
        type = "feir")
  )
})

test_that("a custom impact matrix reproduces the orthogonalised response", {
  model <- fx_var_fitted()

  # Under "oir" the recursion normalises the Choleski factor to a unit shock.
  # The custom path normalises nothing, so the same responses need it done
  # here -- which is the asymmetry worth pinning down.
  impact <- lapply(bvartools:::.collect_draws(model), function(draw) {
    p <- t(chol(draw[["Sigma"]]))
    p %*% diag(1 / diag(p))
  })

  expect_equal(
    irf(model, impulse = "income", response = "cons", n_ahead = 5,
        type = "custom", impact = impact),
    irf(model, impulse = "income", response = "cons", n_ahead = 5,
        type = "oir")
  )
})

test_that("a custom impact matrix reproduces the orthogonalised decomposition", {
  model <- fx_var_fitted()

  # .vardecomp does not normalise the Choleski factor, so here the unmodified
  # one is what reproduces "oir".
  impact <- lapply(bvartools:::.collect_draws(model),
                   function(draw) t(chol(draw[["Sigma"]])))

  expect_equal(
    fevd(model, response = "cons", n_ahead = 5, type = "custom",
         impact = impact),
    fevd(model, response = "cons", n_ahead = 5, type = "oir")
  )
})

test_that("a custom impact matrix reproduces the orthogonalised spillovers", {
  model <- fx_var_fitted()
  impact <- lapply(bvartools:::.collect_draws(model),
                   function(draw) t(chol(draw[["Sigma"]])))

  custom <- spillover(model, n_ahead = 5, type = "custom", impact = impact)
  oir <- spillover(model, n_ahead = 5, type = "oir")

  # Everything but the record of which type was asked for.
  expect_equal(custom[["total"]], oir[["total"]])
  expect_equal(custom[["table"]], oir[["table"]])
  expect_equal(custom[["net"]], oir[["net"]])
})

test_that("a rotation of the Choleski factor still decomposes the variance", {
  model <- fx_var_fitted()
  k <- model[["model"]][["k"]]

  set.seed(4711)
  rotation <- qr.Q(qr(matrix(stats::rnorm(k * k), k)))
  impact <- lapply(bvartools:::.collect_draws(model),
                   function(draw) t(chol(draw[["Sigma"]])) %*% rotation)

  decomp <- fevd(model, response = "cons", n_ahead = 5, type = "custom",
                 impact = impact)

  # P P' = Sigma survives an orthogonal rotation, which is what makes the
  # shares add up. This is the property a sign restricted identification will
  # rely on.
  expect_equal(as.numeric(rowSums(decomp)), rep(1, nrow(decomp)))

  # It moves weight between the shocks, so it is not the orthogonalised
  # decomposition under another name.
  oir <- fevd(model, response = "cons", n_ahead = 5, type = "oir")
  expect_false(isTRUE(all.equal(unclass(decomp), unclass(oir))))
})

test_that("a custom identification is checked before it reaches the recursion", {
  model <- fx_var_fitted()
  k <- model[["model"]][["k"]]

  expect_error(irf(model, impulse = "income", response = "cons",
                   type = "custom"),
               "need an impact matrix")
  expect_error(fevd(model, response = "cons", type = "custom"),
               "needs an impact matrix")
  expect_error(spillover(model, type = "custom"),
               "need an impact matrix")

  expect_error(irf(model, impulse = "income", response = "cons",
                   type = "custom", impact = diag(1, k + 1)),
               "matrices")
  expect_error(irf(model, impulse = "income", response = "cons",
                   type = "custom", impact = list(diag(1, k))),
               "one per posterior draw")
  expect_error(irf(model, impulse = "income", response = "cons",
                   type = "custom", impact = "chol"),
               "matrix, a list of matrices or a function")

  # The size of a custom shock is carried by the impact matrix, so there is no
  # standard deviation left for the recursion to read off Sigma.
  expect_error(irf(model, impulse = "income", response = "cons",
                   type = "custom", impact = diag(1, k), shock = "sd"),
               "must be numeric")
})

test_that("a structural model rejects a custom identification", {
  expect_error(irf(fx_svar_fitted(), impulse = "income", response = "cons",
                   type = "custom", impact = diag(1, 3)),
               "not defined for a structural model")
})

test_that("the workers reject a draw they cannot identify", {
  draw <- bvartools:::.collect_draws(fx_var_fitted())[[1]]
  draw[["shock"]] <- 1

  expect_error(bvartools:::.ir(draw, h = 2, type = "nonsense", impulse = 1,
                               response = 1),
               "unknown type")
  expect_error(bvartools:::.vardecomp(draw, h = 2, type = "nonsense",
                                      response = 1),
               "unknown type")
  expect_error(bvartools:::.spillover_table(draw, h = 2, type = "nonsense"),
               "unknown type")

  # Reachable only by calling a worker directly -- the exported functions check
  # for the matrix first -- but the workers are the last line and say so.
  expect_error(bvartools:::.vardecomp(draw, h = 2, type = "custom",
                                      response = 1),
               "element 'P'")
})
