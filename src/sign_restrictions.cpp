// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <vector>

//' One rotation of a draw that satisfies a set of sign restrictions
//'
//' Draws orthogonal matrices \eqn{Q} until the impulse responses of the
//' rotated model \eqn{\Phi_i P Q}, with \eqn{P} the Choleski factor of the
//' error covariance, carry the signs the caller asked for. Every such rotation
//' leaves \eqn{P Q (P Q)^{\prime} = \Sigma}, so the rotated model describes the
//' data exactly as well as the draw it came from: the restrictions choose among
//' models the likelihood cannot tell apart, which is what makes the
//' identification set valued rather than a point.
//'
//' The forecast error responses do not involve \eqn{Q}, so they are built once
//' per draw and every attempt reads them. That is what makes a high rejection
//' rate affordable.
//'
//' @param A a list with elements \code{A}, the k x kp coefficients of one draw,
//'   and \code{Sigma}, its k x k error covariance.
//' @param restrictions a matrix with one row per restriction and four columns:
//'   the shock, the response variable, the sign, and the horizon. The first two
//'   are counted from one, the horizon from zero.
//' @param max_tries the number of rotations to try before giving up on the
//'   draw.
//'
//' @return The accepted k x k rotation, or a 0 x 0 matrix if none was found.
//'
//' @noRd
// [[Rcpp::export(.draw_sign_restricted_q)]]
arma::mat draw_sign_restricted_q(Rcpp::List A, const arma::mat &restrictions,
                                 int max_tries) {

  arma::mat coef = Rcpp::as<arma::mat>(A["A"]);
  arma::mat sigma = Rcpp::as<arma::mat>(A["Sigma"]);

  const int k = coef.n_rows;
  const int n_restrictions = restrictions.n_rows;

  arma::mat chol_upper;
  if (!arma::chol(chol_upper, sigma)) {
    Rcpp::stop("draw_sign_restricted_q: the error covariance of a draw is not positive definite.");
  }
  const arma::mat P = arma::trans(chol_upper);

  const int h_max = static_cast<int>(restrictions.col(3).max());

  // Theta_h = Phi_h P, stacked by horizon. The rotation enters only as
  // Theta_h Q, so this is the whole of the draw that the attempts below need.
  arma::mat theta = arma::zeros<arma::mat>((h_max + 1) * k, k);
  arma::mat phi_store = arma::zeros<arma::mat>((h_max + 1) * k, k);
  arma::mat phi = arma::eye<arma::mat>(k, k);
  phi_store.rows(0, k - 1) = phi;
  theta.rows(0, k - 1) = P;

  if (h_max > 0) {
    // Lags beyond the furthest restricted horizon never enter the recursion,
    // the economy .ir and .vardecomp make as well.
    int n_use = coef.n_cols;
    const int lags = coef.n_cols / k;
    if (h_max < lags) {
      n_use = h_max * k;
    }

    arma::mat coef_use = arma::zeros<arma::mat>(k, h_max * k);
    coef_use.cols(0, n_use - 1) = coef.cols(0, n_use - 1);

    for (int i = 1; i <= h_max; i++) {
      phi.zeros();
      for (int j = 1; j <= i; j++) {
        phi += phi_store.rows((i - j) * k, (i - j + 1) * k - 1) *
               coef_use.cols((j - 1) * k, j * k - 1);
      }
      phi_store.rows(i * k, (i + 1) * k - 1) = phi;
      theta.rows(i * k, (i + 1) * k - 1) = phi * P;
    }
  }

  // The columns of Q the restrictions speak about. One that carries none is
  // left as drawn: nothing in the restrictions distinguishes one rotation of
  // an unrestricted shock from another, and pretending otherwise would narrow
  // the identified set by an assumption nobody made.
  std::vector<int> shocks;
  for (int r = 0; r < n_restrictions; r++) {
    const int j = static_cast<int>(restrictions(r, 0)) - 1;
    if (std::find(shocks.begin(), shocks.end(), j) == shocks.end()) {
      shocks.push_back(j);
    }
  }

  arma::mat draw(k, k), q_factor, r_factor, q_candidate;

  for (int attempt = 0; attempt < max_tries; attempt++) {

    // Haar uniform over the orthogonal group: the Q of a QR decomposition of a
    // standard normal matrix, with each column signed by the diagonal of R so
    // that the decomposition is unique and the draw therefore uniform.
    draw.randn();
    if (!arma::qr(q_factor, r_factor, draw)) {
      continue;
    }
    arma::vec diag_sign = arma::sign(r_factor.diag());
    diag_sign.replace(0.0, 1.0);
    q_candidate = q_factor * arma::diagmat(diag_sign);

    bool accepted = true;

    for (size_t s = 0; s < shocks.size() && accepted; s++) {
      const int j = shocks[s];

      // A sign restriction cannot tell a shock from its own negative, so a
      // column that fails is retried with its sign flipped before the whole
      // rotation is discarded. Negating one column leaves Q orthogonal, and
      // without this step rotations that satisfy the restrictions under the
      // other labelling of the shock would be rejected -- which both wastes
      // attempts and skews the set that survives.
      bool as_drawn = true;
      bool flipped = true;

      for (int r = 0; r < n_restrictions; r++) {
        if (static_cast<int>(restrictions(r, 0)) - 1 != j) {
          continue;
        }

        const int response = static_cast<int>(restrictions(r, 1)) - 1;
        const int horizon = static_cast<int>(restrictions(r, 3));
        const double value = arma::dot(theta.row(horizon * k + response),
                                       q_candidate.col(j));

        if (restrictions(r, 2) > 0) {
          if (value <= 0) as_drawn = false;
          if (value >= 0) flipped = false;
        } else {
          if (value >= 0) as_drawn = false;
          if (value <= 0) flipped = false;
        }
      }

      if (!as_drawn) {
        if (flipped) {
          q_candidate.col(j) *= -1;
        } else {
          accepted = false;
        }
      }
    }

    if (accepted) {
      return q_candidate;
    }
  }

  return arma::mat();
}
