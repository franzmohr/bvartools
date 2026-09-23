// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include <cmath>
#include <limits>
#include <vector>

// Inference in an SVAR identified with sign and zero restrictions, after
// Arias, Rubio-Ramirez and Waggoner (2018).
//
// The rejection sampler in sign_restrictions.cpp cannot impose a zero: the set
// of rotations that satisfies one has probability zero, so no number of
// attempts finds it. What this file implements instead is their Algorithm 2,
// which draws a rotation from the set that satisfies the zeros by construction,
// together with the importance weight of their Algorithm 3, which corrects the
// distribution that construction induces back to the posterior the researcher
// asked for.
//
// The whole thing rests on one observation. The restrictions are imposed on
// F(A_0, A_+), the responses stacked over the horizons they speak about, and
// that function has the property F(A_0 Q, A_+ Q) = F(A_0, A_+) Q. So a zero
// restriction Z_j F(A_0 Q, A_+ Q) e_j = 0 is a *linear* restriction on the jth
// column of Q, once the reduced form is fixed. Columns can therefore be built
// one at a time, each drawn uniformly from the sphere in the subspace left by
// the zeros and by the columns already placed.
//
// That construction does not draw from the posterior conditional on the zeros:
// it draws from a distribution whose density differs from it by a volume
// element, the Jacobian factor of a map between two manifolds of different
// dimension. Algorithm 3 computes that volume element numerically and divides
// it out. Section 7 of the paper is worth reading before touching the cost of
// anything here: in their timings the volume element is 4795 of the 4806
// seconds a run takes, and everything else together is eleven.

namespace {

const double kNegInf = -std::numeric_limits<double>::infinity();

// Everything about the identification that does not change from draw to draw.
//
// `w` is the one piece that is neither data nor restriction. Algorithm 2 needs
// an orthonormal basis for the null space of the constraints on each column,
// and any basis will do for the draw itself -- but the volume element
// differentiates the map that produces it, so the choice has to be a fixed,
// almost-everywhere differentiable function rather than whatever a
// decomposition happens to return. Appendix A.3 makes it one by appending a
// fixed random matrix to the constraints and taking a QR decomposition of the
// result. That matrix is drawn once per run, in R, and travels here.
struct ArwSetup {
  int k = 0;                      // endogenous variables
  int m = 0;                      // regressors per equation, deterministic terms included
  int lags = 0;
  std::vector<double> horizons;   // the horizons F stacks, R_PosInf for the long run
  std::vector<arma::mat> z;       // k blocks, z_j x (nh * k): the zero restrictions
  std::vector<arma::mat> s;       // k blocks, s_j x (nh * k): the sign restrictions
  std::vector<arma::mat> w;       // k blocks, (k - (j - 1 + z_j)) x k
  int dim = 0;                    // the w blocks' rows together: where the spheres live
  int nzeros = 0;
};

// F(A_0, A_+): the responses the restrictions are written on, one k x k block
// per horizon, each transposed so that its jth column is the response to the
// jth shock. Horizon zero is the impact response inv(A_0)' and R_PosInf is the
// long run inv(A_0 - A_1 - ... - A_p)'.
arma::mat stacked_responses(const arma::mat &a0, const arma::mat &aplus,
                            const ArwSetup &setup) {

  const int k = setup.k;
  const int nh = static_cast<int>(setup.horizons.size());

  arma::mat out(nh * k, k, arma::fill::zeros);
  if (nh == 0) {
    return out;
  }

  int max_finite = 0;
  bool long_run = false;
  for (int i = 0; i < nh; i++) {
    if (std::isinf(setup.horizons[i])) {
      long_run = true;
    } else if (static_cast<int>(setup.horizons[i]) > max_finite) {
      max_finite = static_cast<int>(setup.horizons[i]);
    }
  }

  // The lag blocks of A_+, which both recursions below read. Rows beyond
  // lags * k are the deterministic terms and never enter a response.
  std::vector<arma::mat> a(setup.lags);
  for (int i = 0; i < setup.lags; i++) {
    a[i] = aplus.rows(i * k, (i + 1) * k - 1);
  }

  // R(i) = (R(0) A_i + R(1) A_{i-1} + ... + R(i-1) A_1) inv(A_0), with the sum
  // running over the lags that exist. Only horizons that are actually asked for
  // are reached: a restriction at horizon zero alone computes nothing here.
  std::vector<arma::mat> r(max_finite + 1);
  const int first = std::min(max_finite, setup.lags);
  for (int i = 1; i <= first; i++) {
    arma::mat x = arma::solve(a0, a[i - 1]);
    for (int j = 1; j < i; j++) {
      x += r[j] * a[i - j - 1];
    }
    r[i] = arma::solve(a0.t(), x.t()).t();
  }
  for (int i = setup.lags + 1; i <= max_finite; i++) {
    arma::mat x(k, k, arma::fill::zeros);
    if (setup.lags > 0) {
      x = r[i - setup.lags] * a[setup.lags - 1];
    }
    for (int j = 1; j < setup.lags; j++) {
      x += r[i - j] * a[j - 1];
    }
    r[i] = arma::solve(a0.t(), x.t()).t();
  }

  arma::mat r_long_run;
  if (long_run) {
    arma::mat total = a0;
    for (int j = 0; j < setup.lags; j++) {
      total -= a[j];
    }
    r_long_run = arma::inv(total);
  }

  for (int i = 0; i < nh; i++) {
    if (std::isinf(setup.horizons[i])) {
      out.rows(i * k, (i + 1) * k - 1) = r_long_run.t();
    } else if (setup.horizons[i] == 0.0) {
      out.rows(i * k, (i + 1) * k - 1) = arma::inv(a0).t();
    } else {
      out.rows(i * k, (i + 1) * k - 1) = r[static_cast<int>(setup.horizons[i])].t();
    }
  }

  return out;
}

// Z_j F, one block per shock. The zero restrictions on column j of Q read
// Z_j F(B, Sigma, I) Q e_j = 0, so these blocks are what makes them linear.
std::vector<arma::mat> zf_blocks(const arma::mat &f, const ArwSetup &setup) {
  std::vector<arma::mat> out(setup.k);
  for (int j = 0; j < setup.k; j++) {
    out[j] = setup.z[j] * f;
  }
  return out;
}

// An orthonormal basis for the null space of the constraints on column j,
// signed so that it is a function of them rather than of a decomposition's
// conventions. `constraints` is the (j - 1 + z_j) x k matrix of columns already
// placed and zero restrictions; `fixed` is the w block that fills it out to a
// square matrix whose QR decomposition therefore has the null space in its
// trailing columns.
bool null_basis(const arma::mat &constraints, const arma::mat &fixed,
                arma::mat &basis) {

  const int k = static_cast<int>(fixed.n_cols);
  const int s = static_cast<int>(fixed.n_rows);

  arma::mat stacked(k, k);
  if (constraints.n_rows > 0) {
    stacked.rows(0, constraints.n_rows - 1) = constraints;
  }
  stacked.rows(k - s, k - 1) = fixed;

  arma::mat q_factor, r_factor;
  if (!arma::qr(q_factor, r_factor, stacked.t())) {
    return false;
  }

  for (int i = k - s; i < k; i++) {
    if (r_factor(i, i) < 0.0) {
      q_factor.col(i) = -q_factor.col(i);
    }
  }

  basis = q_factor.cols(k - s, k - 1);
  return true;
}

// The constraints on column j: the columns already placed, then the zeros.
arma::mat column_constraints(const arma::mat &q, int j, const arma::mat &zf) {
  arma::mat out(j + zf.n_rows, q.n_rows);
  if (j > 0) {
    out.rows(0, j - 1) = q.cols(0, j - 1).t();
  }
  if (zf.n_rows > 0) {
    out.rows(j, j + zf.n_rows - 1) = zf;
  }
  return out;
}

// Steps 2 and 3 of Algorithm 2: the unit vectors in `w` become the columns of
// an orthogonal matrix that satisfies the zero restrictions.
bool spheres_to_q(const arma::vec &w, const std::vector<arma::mat> &zf,
                  const ArwSetup &setup, arma::mat &q) {

  const int k = setup.k;
  q.zeros(k, k);

  int at = 0;
  for (int j = 0; j < k; j++) {
    const int s = static_cast<int>(setup.w[j].n_rows);
    arma::mat basis;
    if (!null_basis(column_constraints(q, j, zf[j]), setup.w[j], basis)) {
      return false;
    }
    q.col(j) = basis * w.subvec(at, at + s - 1);
    at += s;
  }

  return true;
}

// Its inverse. The same basis is rebuilt column by column and the column is
// resolved in it, which is what makes ff_h below invertible and so
// differentiable in the sense the volume element needs.
bool q_to_spheres(const arma::mat &q, const std::vector<arma::mat> &zf,
                  const ArwSetup &setup, arma::vec &w) {

  const int k = setup.k;
  w.zeros(setup.dim);

  int at = 0;
  for (int j = 0; j < k; j++) {
    const int s = static_cast<int>(setup.w[j].n_rows);
    arma::mat basis;
    if (!null_basis(column_constraints(q, j, zf[j]), setup.w[j], basis)) {
      return false;
    }
    w.subvec(at, at + s - 1) = basis.t() * q.col(j);
    at += s;
  }

  return true;
}

// The Choleski factor the paper writes h(Sigma): upper triangular with
// h(Sigma)' h(Sigma) = Sigma, symmetrised first because the numerical
// derivative below perturbs Sigma one element at a time and so walks off the
// symmetric matrices.
bool chol_tilde(const arma::mat &sigma, arma::mat &out) {
  return arma::chol(out, 0.5 * (sigma + sigma.t()));
}

// f_h^{-1}: the orthogonal reduced form (B, Sigma, Q) to the structural
// parameters (A_0, A_+).
bool orthogonal_to_structural(const arma::mat &b, const arma::mat &sigma,
                              const arma::mat &q, arma::mat &a0, arma::mat &aplus) {
  arma::mat h;
  if (!chol_tilde(sigma, h)) {
    return false;
  }
  a0 = arma::solve(arma::trimatu(h), q);
  aplus = b * a0;
  return true;
}

// f_h, the other way.
bool structural_to_orthogonal(const arma::mat &a0, const arma::mat &aplus,
                              arma::mat &b, arma::mat &sigma, arma::mat &q) {
  b = arma::solve(a0.t(), aplus.t()).t();
  arma::mat gram = a0 * a0.t();
  if (!arma::inv_sympd(sigma, 0.5 * (gram + gram.t()))) {
    return false;
  }
  arma::mat h;
  if (!chol_tilde(sigma, h)) {
    return false;
  }
  q = h * a0;
  return true;
}

// The zero restrictions as a function of the structural parameters, which is
// the map whose null space the volume element is taken over. Evaluated at the
// draw it returns zeros; what the derivative of it describes is the manifold
// the draw lies on.
bool zero_restrictions(const arma::vec &x, const ArwSetup &setup, arma::vec &out) {

  const int k = setup.k;
  const arma::mat a0 = arma::reshape(x.subvec(0, k * k - 1), k, k);
  const arma::mat aplus = arma::reshape(x.subvec(k * k, x.n_elem - 1), setup.m, k);

  const arma::mat f = stacked_responses(a0, aplus, setup);

  out.zeros(setup.nzeros);
  int at = 0;
  for (int j = 0; j < k; j++) {
    const int s = static_cast<int>(setup.z[j].n_rows);
    if (s == 0) {
      continue;
    }
    out.subvec(at, at + s - 1) = setup.z[j] * f.col(j);
    at += s;
  }

  return true;
}

// ff_h: the structural parameters to (B, Sigma, w). This is the map whose
// volume element, restricted to the zero manifold, Algorithm 3 divides by --
// it is Algorithm 2 read backwards, so differentiating it is differentiating
// the construction the draw came out of.
bool structural_to_spheres(const arma::vec &x, const ArwSetup &setup, arma::vec &out) {

  const int k = setup.k;
  const int m = setup.m;

  const arma::mat a0 = arma::reshape(x.subvec(0, k * k - 1), k, k);
  const arma::mat aplus = arma::reshape(x.subvec(k * k, x.n_elem - 1), m, k);

  arma::mat b, sigma, q;
  if (!structural_to_orthogonal(a0, aplus, b, sigma, q)) {
    return false;
  }

  // The blocks are taken at Q = I, where the restrictions become linear in the
  // columns of Q -- the same point Algorithm 2 builds them at.
  arma::mat a0_identity, aplus_identity;
  if (!orthogonal_to_structural(b, sigma, arma::eye(k, k), a0_identity, aplus_identity)) {
    return false;
  }
  const std::vector<arma::mat> zf =
    zf_blocks(stacked_responses(a0_identity, aplus_identity, setup), setup);

  arma::vec w;
  if (!q_to_spheres(q, zf, setup, w)) {
    return false;
  }

  out.set_size(m * k + k * k + setup.dim);
  out.subvec(0, m * k - 1) = arma::vectorise(b);
  out.subvec(m * k, m * k + k * k - 1) = arma::vectorise(sigma);
  if (setup.dim > 0) {
    out.subvec(m * k + k * k, out.n_elem - 1) = w;
  }

  return true;
}

// Two-sided unless `one_sided`, which the paper reports as forty percent
// faster and is the caller's choice because it is also less accurate.
template <typename Fn>
bool numerical_derivative(Fn fn, const arma::vec &x, double epsilon,
                          bool one_sided, arma::mat &out) {

  const int n = static_cast<int>(x.n_elem);
  arma::vec z = x;
  arma::vec base;

  if (one_sided && !fn(x, base)) {
    return false;
  }

  for (int j = 0; j < n; j++) {
    arma::vec high, low;

    z(j) = x(j) + epsilon;
    if (!fn(z, high)) {
      return false;
    }
    if (one_sided) {
      low = base;
    } else {
      z(j) = x(j) - epsilon;
      if (!fn(z, low)) {
        return false;
      }
    }
    z(j) = x(j);

    if (j == 0) {
      out.set_size(high.n_elem, n);
    }
    out.col(j) = (high - low) / (one_sided ? epsilon : 2.0 * epsilon);
  }

  return true;
}

// The log of the volume element of `fn` restricted to the manifold the zero
// restrictions define, which is Theorem 3 of the paper evaluated numerically.
// With no zero restrictions the manifold is the whole space and this is the
// unrestricted volume element of Theorem 2.
double log_volume_element(const ArwSetup &setup, const arma::vec &x,
                          double epsilon, bool one_sided) {

  arma::mat d_fn;
  if (!numerical_derivative(
        [&setup](const arma::vec &v, arma::vec &y) {
          return structural_to_spheres(v, setup, y);
        },
        x, epsilon, one_sided, d_fn)) {
    return NA_REAL;
  }

  arma::mat n_mat;
  if (setup.nzeros > 0) {
    arma::mat d_r;
    if (!numerical_derivative(
          [&setup](const arma::vec &v, arma::vec &y) {
            return zero_restrictions(v, setup, y);
          },
          x, epsilon, one_sided, d_r)) {
      return NA_REAL;
    }
    const arma::mat basis = arma::null(d_r);
    if (basis.n_cols == 0) {
      return NA_REAL;
    }
    n_mat = d_fn * basis;
  } else {
    n_mat = d_fn;
  }

  const arma::mat gram = n_mat.t() * n_mat;
  double value = 0.0;
  if (!arma::log_det_sympd(value, 0.5 * (gram + gram.t()))) {
    return NA_REAL;
  }

  return 0.5 * value;
}

// The restriction blocks of one R list, checked into the shape the rest of this
// file assumes.
std::vector<arma::mat> restriction_blocks(const Rcpp::List &blocks, int k,
                                          int cols, const char *name) {

  if (static_cast<int>(blocks.size()) != k) {
    Rcpp::stop("arw_draw_q: '%s' must hold one matrix per shock, so %d of them, not %d.",
               name, k, static_cast<int>(blocks.size()));
  }

  std::vector<arma::mat> out(k);
  for (int j = 0; j < k; j++) {
    if (Rf_isNull(blocks[j])) {
      out[j] = arma::mat(0, cols);
      continue;
    }
    out[j] = Rcpp::as<arma::mat>(blocks[j]);
    if (out[j].n_rows > 0 && static_cast<int>(out[j].n_cols) != cols) {
      Rcpp::stop("arw_draw_q: block %d of '%s' has %d columns, but the responses it "
                 "multiplies have %d rows.",
                 j + 1, name, static_cast<int>(out[j].n_cols), cols);
    }
    if (out[j].n_rows == 0) {
      out[j] = arma::mat(0, cols);
    }
  }

  return out;
}

ArwSetup read_setup(const Rcpp::List &setup_list) {

  ArwSetup setup;
  setup.k = Rcpp::as<int>(setup_list["k"]);
  setup.m = Rcpp::as<int>(setup_list["m"]);
  setup.lags = Rcpp::as<int>(setup_list["lags"]);

  const arma::vec horizons = Rcpp::as<arma::vec>(setup_list["horizons"]);
  setup.horizons.assign(horizons.begin(), horizons.end());
  const int cols = static_cast<int>(setup.horizons.size()) * setup.k;

  setup.z = restriction_blocks(Rcpp::List(setup_list["z"]), setup.k, cols, "z");
  setup.s = restriction_blocks(Rcpp::List(setup_list["s"]), setup.k, cols, "s");

  setup.nzeros = 0;
  for (int j = 0; j < setup.k; j++) {
    setup.nzeros += static_cast<int>(setup.z[j].n_rows);
  }

  const Rcpp::List w = Rcpp::List(setup_list["w"]);
  if (static_cast<int>(w.size()) != setup.k) {
    Rcpp::stop("arw_draw_q: 'w' must hold one matrix per shock, so %d of them, not %d.",
               setup.k, static_cast<int>(w.size()));
  }
  setup.w.resize(setup.k);
  setup.dim = 0;
  for (int j = 0; j < setup.k; j++) {
    setup.w[j] = Rcpp::as<arma::mat>(w[j]);
    const int expected = setup.k - (j + static_cast<int>(setup.z[j].n_rows));
    if (expected < 1) {
      Rcpp::stop("arw_draw_q: shock %d carries %d zero restrictions, which with the %d "
                 "shocks before it leaves nothing of its column to draw. Algorithm 2 "
                 "needs fewer than %d zero restrictions on it.",
                 j + 1, static_cast<int>(setup.z[j].n_rows), j, setup.k - j);
    }
    if (static_cast<int>(setup.w[j].n_rows) != expected ||
        static_cast<int>(setup.w[j].n_cols) != setup.k) {
      Rcpp::stop("arw_draw_q: block %d of 'w' is %dx%d, but the zero restrictions on that "
                 "shock leave a null space of dimension %d.",
                 j + 1, static_cast<int>(setup.w[j].n_rows),
                 static_cast<int>(setup.w[j].n_cols), expected);
    }
    setup.dim += expected;
  }

  return setup;
}

} // namespace

//' One draw of the sign and zero restricted rotation, with its importance weight
//'
//' Draws a rotation that satisfies the zero restrictions by construction,
//' following Algorithm 2 of Arias, Rubio-Ramirez and Waggoner (2018), and
//' returns the log of the unnormalised importance weight of their Algorithm 3,
//' which corrects the distribution that construction induces back to the
//' posterior conditional on the sign and zero restrictions.
//'
//' Unlike the rejection sampler in \code{.draw_sign_restricted_q}, there is one
//' rotation per posterior draw and no second attempt: a draw whose rotation
//' fails the sign restrictions gets weight zero rather than another try, and a
//' shock is not retried with its sign flipped, since the sphere the column is
//' drawn from already covers both signs of it.
//'
//' @param A a list with elements \code{A}, the k x m coefficients of one draw
//'   with the deterministic terms included, and \code{Sigma}, its k x k error
//'   covariance.
//' @param setup a list with the elements \code{k}, \code{m}, \code{lags},
//'   \code{horizons}, the zero and sign restriction blocks \code{z} and
//'   \code{s}, and \code{w}, the fixed matrices that complete each column's
//'   constraints to a square system. See \code{.arw_setup}.
//' @param weight logical. Should the importance weight be computed? It is the
//'   expensive part of the algorithm by an order of magnitude and is needed only
//'   for a draw that satisfies the sign restrictions.
//' @param epsilon the step of the numerical derivative.
//' @param one_sided logical. Should the numerical derivative be one sided?
//'
//' @return A list with the accepted k x k rotation in \code{q}, or a 0 x 0
//'   matrix when the draw fails the sign restrictions, and the log of the
//'   unnormalised importance weight in \code{log_weight}, which is \code{NA}
//'   when it was not asked for or could not be computed.
//'
//' @noRd
// [[Rcpp::export(.arw_draw_q)]]
Rcpp::List arw_draw_q(Rcpp::List A, Rcpp::List setup_list, bool weight = true,
                      double epsilon = 1e-6, bool one_sided = false) {

  const ArwSetup setup = read_setup(setup_list);
  const int k = setup.k;
  const int m = setup.m;

  const arma::mat coef = Rcpp::as<arma::mat>(A["A"]);
  const arma::mat sigma = Rcpp::as<arma::mat>(A["Sigma"]);

  if (static_cast<int>(coef.n_rows) != k || static_cast<int>(coef.n_cols) != m) {
    Rcpp::stop("arw_draw_q: the coefficients of a draw are %dx%d, but the model has %d "
               "endogenous variables and %d regressors.",
               static_cast<int>(coef.n_rows), static_cast<int>(coef.n_cols), k, m);
  }

  // The paper's B is regressors by equations, which is the transpose of the way
  // the rest of this package carries a draw.
  const arma::mat b = coef.t();

  Rcpp::List failed = Rcpp::List::create(Rcpp::Named("q") = arma::mat(0, 0),
                                         Rcpp::Named("log_weight") = NA_REAL);

  // Step 2 of Algorithm 2 evaluated at Q = I, where the zero restrictions are
  // linear in the columns of the rotation.
  arma::mat a0_identity, aplus_identity;
  if (!orthogonal_to_structural(b, sigma, arma::eye(k, k), a0_identity, aplus_identity)) {
    return failed;
  }
  const std::vector<arma::mat> zf =
    zf_blocks(stacked_responses(a0_identity, aplus_identity, setup), setup);

  // Uniform on the product of spheres, one per shock, each of the dimension the
  // zeros and the columns already placed leave over.
  arma::vec w(setup.dim);
  int at = 0;
  for (int j = 0; j < k; j++) {
    const int s = static_cast<int>(setup.w[j].n_rows);
    arma::vec wj(s, arma::fill::randn);
    const double norm = arma::norm(wj);
    if (!(norm > 0.0)) {
      return failed;
    }
    w.subvec(at, at + s - 1) = wj / norm;
    at += s;
  }

  arma::mat q;
  if (!spheres_to_q(w, zf, setup, q)) {
    return failed;
  }

  arma::mat a0, aplus;
  if (!orthogonal_to_structural(b, sigma, q, a0, aplus)) {
    return failed;
  }

  // Step 2 of Algorithm 3. A draw that fails is kept out of the sample by its
  // weight rather than by being dropped here, so that the effective sample size
  // the caller computes counts the draws it should.
  const arma::mat f = stacked_responses(a0, aplus, setup);
  for (int j = 0; j < k; j++) {
    if (setup.s[j].n_rows == 0) {
      continue;
    }
    const arma::vec signs = setup.s[j] * f.col(j);
    if (signs.min() <= 0.0) {
      return Rcpp::List::create(Rcpp::Named("q") = arma::mat(0, 0),
                                Rcpp::Named("log_weight") = kNegInf);
    }
  }

  double log_weight = NA_REAL;
  if (weight) {
    arma::vec x(k * k + m * k);
    x.subvec(0, k * k - 1) = arma::vectorise(a0);
    x.subvec(k * k, x.n_elem - 1) = arma::vectorise(aplus);

    // The volume element of f_h, which Proposition 1 gives in closed form,
    // against the volume element of Algorithm 2's construction, which does not
    // have one. The normal-inverse-Wishart density appears in both the target
    // and the proposal and cancels, which is why nothing about the prior
    // reaches this function.
    double log_abs_det = 0.0, det_sign = 0.0;
    if (!arma::log_det(log_abs_det, det_sign, a0) || !std::isfinite(log_abs_det)) {
      return failed;
    }
    const double log_ve_fh =
      0.5 * k * (k + 1) * std::log(2.0) - (2.0 * k + m + 1.0) * log_abs_det;

    const double log_ve_gfh = log_volume_element(setup, x, epsilon, one_sided);
    if (std::isfinite(log_ve_gfh)) {
      log_weight = log_ve_fh - log_ve_gfh;
    }
  }

  return Rcpp::List::create(Rcpp::Named("q") = q,
                            Rcpp::Named("log_weight") = log_weight);
}
