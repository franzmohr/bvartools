// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

//' Connectedness table of one posterior draw
//'
//' The normalised forecast error variance decomposition table behind the
//' spillover index of Diebold and Yilmaz (2012): row \eqn{j} column \eqn{k} is
//' the share of the \eqn{h} step forecast error variance of variable \eqn{j}
//' that is attributed to a shock to variable \eqn{k}, and every row sums to one.
//'
//' Separate from \code{.vardecomp} rather than a loop over it, for two reasons.
//' The forecast error impulse responses do not depend on the response variable,
//' so filling the whole table in one pass costs what \code{.vardecomp} costs for
//' a single row of it. And the generalised decomposition here divides by the
//' variance of the \emph{impulse} variable, as in Pesaran and Shin (1998), which
//' is the quantity the spillover index is defined on; \code{.vardecomp} divides
//' by the standard deviation of the response instead, a factor that cancels
//' when a row is normalised and so cannot serve here.
//'
//' @param A a list with elements \code{A}, the k x kp coefficients of one draw,
//'   and \code{Sigma}, its k x k error covariance.
//' @param h an integer of the forecast horizon, at least one.
//' @param type either \code{"gir"} or \code{"oir"}.
//'
//' @noRd
// [[Rcpp::export(.spillover_table)]]
arma::mat spillover_table(Rcpp::List A, int h, std::string type) {

  arma::mat a = Rcpp::as<arma::mat>(A["A"]);
  arma::mat sigma = Rcpp::as<arma::mat>(A["Sigma"]);

  const int k = a.n_rows;

  // Regressor columns the recursion can reach. Lags beyond the horizon never
  // enter it, which is the economy .vardecomp makes as well.
  int n_use = a.n_cols;
  const int lags = a.n_cols / k;
  if (h < lags) {
    n_use = h * k;
  }

  arma::mat P, sigma_mse;
  if (type == "oir") {
    P = arma::trans(arma::chol(sigma));
  } else {
    P = sigma;
  }
  sigma_mse = sigma;

  arma::mat A_temp = arma::zeros<arma::mat>(k, h * k);
  A_temp.cols(0, n_use - 1) = a.cols(0, n_use - 1);

  // The forecast error impulse responses, kept because the recursion needs
  // every earlier one.
  arma::mat phi_store = arma::zeros<arma::mat>((h + 1) * k, k);
  arma::mat phi = arma::eye<arma::mat>(k, k);
  phi_store.rows(0, k - 1) = phi;

  // Numerator of the whole table at once: element (j, i) accumulates
  // (e_j' Phi P e_i)^2. The denominator is the mean squared error of each
  // response, which is the diagonal of Phi Sigma Phi'.
  arma::mat numerator = arma::square(phi * P);
  arma::mat mse_step = phi * sigma_mse * arma::trans(phi);
  arma::vec mse = mse_step.diag();

  for (int i = 1; i <= h; i++) {
    phi.zeros();
    for (int j = 1; j <= i; j++) {
      phi += phi_store.rows((i - j) * k, (i - j + 1) * k - 1) *
             A_temp.cols((j - 1) * k, j * k - 1);
    }
    phi_store.rows(i * k, (i + 1) * k - 1) = phi;

    numerator += arma::square(phi * P);
    mse_step = phi * sigma_mse * arma::trans(phi);
    mse += mse_step.diag();
  }

  arma::mat theta = numerator;

  // Row j by its own forecast error variance.
  theta.each_col() /= mse;

  // Column i by the variance of the shock it carries. This is the Pesaran and
  // Shin scaling, and the one place this differs from .vardecomp.
  if (type != "oir") {
    arma::rowvec shock_var = arma::trans(sigma.diag());
    theta.each_row() /= shock_var;
  }

  // Normalise each row to one. Under "oir" the shares already add up, so this
  // divides by one and is kept only to leave a single code path.
  arma::vec row_total = arma::sum(theta, 1);
  theta.each_col() /= row_total;

  return theta;
}
