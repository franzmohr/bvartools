#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "bayests_r_io.h"
#include "bayests/vec_to_var.h"

// The VEC-to-VAR transformation of the vendored BayesTS core, made reachable
// from R. A VEC and its level VAR are the same model in two parameterisations,
// so this is a change of basis rather than an estimation step: no sampler runs
// here, and there is nothing to seed or to report progress on. See
// inst/include/bayests/vec_to_var.h for what the two entry points do and
// src/core/VENDORED.md for where they come from.
//
// bvecmodel_to_bvarmodel() is the only caller. It assembles the level data
// matrices itself -- the core transforms coefficients, not observations -- and
// uses these two to carry the specification and the posterior across.

namespace {

using namespace bayests_r;

/// The VEC posterior as the core keeps it. Only the parts the transformation
/// reads are taken: the loadings and non-cointegration coefficients in `a`, the
/// cointegration space in `beta`, and the error precision, whose column count is
/// what tells the core how long the chain is. Inclusion indicators are left
/// behind on purpose -- they do not survive the transformation.
bayests::VecNormalWishartDraws read_draws(const Rcpp::List &object) {

  bayests::VecNormalWishartDraws draws;

  if (!has(object, "posterior")) {
    return draws;
  }

  const Rcpp::List posterior = object["posterior"];
  if (has(posterior, "a")) {
    const Rcpp::List posterior_a = posterior["a"];
    read_draws_if_present(posterior_a, "coeffs", draws.a);
  }
  if (has(posterior, "beta")) {
    const Rcpp::List posterior_beta = posterior["beta"];
    read_draws_if_present(posterior_beta, "coeffs", draws.beta);
  }
  if (has(posterior, "u_sigma_inv")) {
    const Rcpp::List posterior_u_sigma_inv = posterior["u_sigma_inv"];
    read_draws_if_present(posterior_u_sigma_inv, "coeffs", draws.u_sigma_inv);
  }

  return draws;
}

} // namespace

// The dimensions of the level VAR implied by a VEC. Takes and returns the
// elements of a model object's `model` list rather than the object itself, so
// that a specification can be converted before any data or posterior exists.
//
// [[Rcpp::export(.VecToVarSpecification)]]
Rcpp::List VecToVarSpecification(Rcpp::List model) {

  const bayests::VarSpec spec = bayests_r::read_spec(model);
  const bayests::VarSpec var_spec = bayests::vec_to_var_spec(spec);

  // `rank`, `k_beta` and `n_restricted` are deliberately not returned: a VAR
  // has no cointegration relation, and leaving them in a 'bvarmodel' would make
  // every reader of that object count k * rank coefficients that are not there.
  return Rcpp::List::create(Rcpp::Named("k") = var_spec.k,
                            Rcpp::Named("p") = var_spec.p,
                            Rcpp::Named("m") = var_spec.m,
                            Rcpp::Named("s") = var_spec.s,
                            Rcpp::Named("n") = var_spec.n,
                            Rcpp::Named("varsel") = std::string(bayests::to_string(var_spec.varsel)),
                            Rcpp::Named("structural") = var_spec.structural);
}

// Posterior draws of a VEC, rewritten in the level VAR parameterisation. The
// whole 'bvecmodel' goes in, because the transformation needs its `model` list
// to know how to cut the draws up.
//
// [[Rcpp::export(.VecToVarCoefficients)]]
Rcpp::List VecToVarCoefficients(Rcpp::List object) {

  const Rcpp::List model = object["model"];
  const bayests::VarSpec spec = bayests_r::read_spec(model);
  const bayests::VecNormalWishartDraws draws = read_draws(object);

  // The core validates the draws against the specification and throws
  // std::invalid_argument naming the first inconsistency; Rcpp turns that into
  // an R error.
  const bayests::VarNormalWishartDraws out =
    bayests::vec_to_var_coefficients(spec, draws);

  // `a_lambda` is empty by construction here, so unlike a sampler's output this
  // never carries a `lambda` element.
  return Rcpp::List::create(
    Rcpp::Named("a") = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(out.a)),
    Rcpp::Named("u_sigma_inv") = Rcpp::List::create(Rcpp::Named("coeffs") = draws_to_r(out.u_sigma_inv)));
}

/*** R

data("e6")

object <- create_bvecmodel(data = e6 * 100,
                           p = 2, r = 1, const = "restricted",
                           error = "wishart",
                           iterations = 20, burnin = 10)

object <- add_priors(object,
                     coef = list(v_i = 1),
                     coint = list(v_i = 0, p_tau_i = 1),
                     sigma = list(df = 3, scale = 1))

object <- add_initial_values(object)
object <- .VecNormalWishartCoefficients(object)

.VecToVarSpecification(object[["model"]])
str(.VecToVarCoefficients(object))

*/
