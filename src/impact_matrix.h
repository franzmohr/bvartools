#ifndef BVARTOOLS_IMPACT_MATRIX_H
#define BVARTOOLS_IMPACT_MATRIX_H

#include <RcppArmadillo.h>

#include <string>

// The impact matrix of one posterior draw, supplied by the caller.
//
// The recursion behind an impulse response, a variance decomposition and a
// spillover table is the same in all three: the forecast error responses
// Phi_i, post-multiplied by a k x k matrix P that says what a shock is. Only
// the choice of P separates the types from one another.
//
// For the types this package derives itself the three workers each build their
// own P, and they deliberately disagree: .ir normalises the Cholesky factor to
// a unit shock and .vardecomp does not, .ir scales a generalised response by
// the variance of the impulse variable while .vardecomp scales by the standard
// deviation of the response, and the structural types carry a factor of Sigma
// in one and not in the other. Those blocks are not duplicates of each other
// and stay where they are.
//
// What they did share was the assumption that P can only come from a name, so
// that every further identification scheme meant another branch in three
// files. Reading P from the draw instead opens the recursion to any scheme
// that can produce an impact matrix -- a sign restricted rotation being the
// case this was added for -- without touching the recursion again.
//
// Only the shape is checked here. Whether P factorises Sigma is the caller's
// business: the impulse responses do not care, and the two decompositions say
// in their own documentation what they assume of it.
inline arma::mat bvartools_impact_matrix(const Rcpp::List &A, int k,
                                         const std::string &caller) {

  if (!A.containsElementNamed("P") || Rf_isNull(A["P"])) {
    Rcpp::stop("%s: type \"custom\" needs the impact matrix of the draw in element 'P'.",
               caller);
  }

  arma::mat P = Rcpp::as<arma::mat>(A["P"]);

  if (static_cast<int>(P.n_rows) != k || static_cast<int>(P.n_cols) != k) {
    Rcpp::stop("%s: the impact matrix in element 'P' is %dx%d, but the model has %d endogenous variables.",
               caller, static_cast<int>(P.n_rows), static_cast<int>(P.n_cols), k);
  }

  if (!P.is_finite()) {
    Rcpp::stop("%s: the impact matrix in element 'P' contains non-finite values.", caller);
  }

  return P;
}

#endif
