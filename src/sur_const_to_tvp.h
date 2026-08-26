#ifndef BVARTOOLS_SUR_CONST_TO_TVP_H
#define BVARTOOLS_SUR_CONST_TO_TVP_H

#include <RcppArmadillo.h>

// Declared for the callers inside this package, which have to call it as a
// plain C++ function.
//
// The generated bvartools::sur_const_to_tvp() is for other packages. It reaches
// the same code through R_GetCCallable inside an Rcpp::RNGScope, and entering
// that scope re-reads .Random.seed -- which rewinds the random number stream of
// whatever sampler is running and hands the same draws out twice. That the
// callee itself draws nothing is no protection; the rewind is in the wrapper.
// See src/core/VENDORED.md.
arma::sp_mat sur_const_to_tvp(const arma::mat& z, const arma::uword& k, const int& tt);

#endif // BVARTOOLS_SUR_CONST_TO_TVP_H
