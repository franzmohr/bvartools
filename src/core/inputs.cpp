// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/inputs.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{

std::string dims(const arma::mat &m)
{
    return std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols);
}

/// Three significant digits in whichever notation suits: std::to_string() prints
/// a fixed six decimals, which turns every asymmetry worth reporting into 0.000000.
std::string number(double x)
{
    char buffer[32];
    std::snprintf(buffer, sizeof buffer, "%.3g", x);
    return buffer;
}

void require_square(const arma::mat &m, arma::uword side, const char *what)
{
    if (m.n_rows != side || m.n_cols != side)
    {
        throw std::invalid_argument(std::string(what) + " must be " + std::to_string(side) +
                                    "x" + std::to_string(side) + ", got " + dims(m));
    }
}

/// How far a matrix may be from its transpose and still count as symmetric: the
/// largest |m(i,j) - m(j,i)| against the largest |m(i,j)|.
///
/// Relative, because rounding scales with the entries. A prior precision of
/// 1e4 I and one of 1e-4 I pick up the same error in proportion, and an absolute
/// bound would be loose for the one and tight for the other.
///
/// And loose, because what rounding leaves is either nothing or a matter of
/// conditioning. An outer product, a Kronecker product and (P + P') / 2 -- how
/// bvartools symmetrises its "ml" prior -- are symmetric to the bit. A matrix
/// inverted through a general LU rather than a Cholesky is off by about its
/// condition number times machine epsilon, so 1e-8 passes such a matrix up to a
/// condition number near 1e7, past the point where its inverse had eight digits
/// to give. It also passes everything bvartools accepts: add_priors() tests with
/// isSymmetric(), whose mean relative difference of 100 * .Machine$double.eps
/// bounds the largest one by n^2 times that, under 1e-8 for any n below 670.
///
/// What it has to stop is a matrix that is not the one meant -- a triangle left
/// empty, a product taken in the wrong order, an element mistyped -- and those
/// are off by many orders of magnitude more.
constexpr double kSymmetryTolerance = 1e-8;

/// Refuses a matrix further from symmetric than rounding explains. For a matrix
/// already known to be square, which is what makes the transpose comparable.
///
/// Refused rather than symmetrised on the way in. The samplers read these
/// matrices partly through one triangle -- draw_normal_precision() factorises
/// the upper -- and partly whole, in products with a mean or a draw, so an
/// asymmetric one runs and stands for a prior that is neither the matrix nor its
/// transpose. Averaging it would pick one of the readings for the file, which
/// is the file's to say.
void require_symmetric(const arma::mat &m, const char *what)
{
    if (m.is_empty())
    {
        return;
    }

    const double largest = arma::abs(m).max();
    const double asymmetry = arma::abs(m - arma::trans(m)).max();

    if (!(asymmetry <= kSymmetryTolerance * largest))
    {
        throw std::invalid_argument(
            std::string(what) + " must be symmetric, but differs from its transpose by up to " +
            number(asymmetry) + " against a largest absolute element of " + number(largest) +
            "; the tolerance is " + number(kSymmetryTolerance) + " times that element");
    }
}

/// Refuses a symmetric matrix that is not positive definite, read off its
/// eigenvalues rather than a Cholesky attempt, which Armadillo may report on the
/// console when it fails.
///
/// For the prior precision of a random walk's state before the sample. The
/// samplers integrate that state out of the prior of the first period, which
/// takes the precision's inverse (see initial_state_variance()), so a zero or
/// singular one -- a flat prior on the state -- has no variance to hand the
/// smoother.
void require_positive_definite(const arma::mat &m, const char *what)
{
    arma::vec eigenvalues;
    if (m.is_empty() || !arma::eig_sym(eigenvalues, m) || !(eigenvalues.min() > 0.0))
    {
        throw std::invalid_argument(
            std::string(what) + " must be positive definite, since the state it describes is "
            "integrated out of the prior of the first period, which takes its inverse" +
            (eigenvalues.n_elem > 0 ? "; its smallest eigenvalue is " + number(eigenvalues.min())
                                    : std::string()));
    }
}

void require_length(const arma::vec &v, arma::uword n, const char *what)
{
    if (v.n_elem != n)
    {
        throw std::invalid_argument(std::string(what) + " must have " + std::to_string(n) +
                                    " elements, got " + std::to_string(v.n_elem));
    }
}

void require_shape(const arma::mat &m, arma::uword rows, arma::uword cols, const char *what)
{
    if (m.n_rows != rows || m.n_cols != cols)
    {
        throw std::invalid_argument(std::string(what) + " must be " + std::to_string(rows) + "x" +
                                    std::to_string(cols) + ", got " + dims(m));
    }
}

/// True if the value is a finite number, read off the exponent bits.
///
/// Not std::isfinite(), for the reason stochvol_mixture.h sets out: a host that
/// compiles these sources with -ffast-math is licensed to fold it to true, and
/// the checks below exist to catch exactly the NaN it would fold away.
bool is_finite(const double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "expected IEEE-754 binary64");

    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
}

/// Refuses an element that is not finite, or not above `floor` -- strictly
/// above it where `strict`. The one test behind every value check below, so the
/// messages say the same thing the same way.
///
/// Shapes were all these files were ever checked for, and a value no prior can
/// have ran regardless: a negative gamma rate went through `bayests check` and a
/// chain with exit code 0 and plausible numbers, since the posterior rate stays
/// positive while the data outweigh it.
void require_above(const arma::vec &v, const double floor, const bool strict,
                   const std::string &what)
{
    for (arma::uword i = 0; i < v.n_elem; i++)
    {
        const double x = v[i];
        if (!is_finite(x) || (strict ? !(x > floor) : !(x >= floor)))
        {
            throw std::invalid_argument(what + " must be finite and " +
                                        (strict ? "greater than " : "at least ") + number(floor) +
                                        ", but element " + std::to_string(i + 1) + " is " +
                                        number(x));
        }
    }
}

/// Refuses a probability outside [0, 1]. A selection block takes log(p) and
/// log(1 - p) of it, and outside that interval one of them is a NaN -- which the
/// indicator draw reads as exclusion, so an inclusion probability of 1.5 ran a
/// whole chain with every coefficient it covered held at zero.
void require_probabilities(const arma::vec &v, const std::string &what)
{
    require_above(v, 0.0, false, what);
    for (arma::uword i = 0; i < v.n_elem; i++)
    {
        if (v[i] > 1.0)
        {
            throw std::invalid_argument(what + " must lie in [0, 1], but element " +
                                        std::to_string(i + 1) + " is " + number(v[i]));
        }
    }
}

/// Refuses a matrix with anything off its diagonal. For a starting value the
/// sampler redraws only the diagonal of, where an off-diagonal element would
/// otherwise stay in the chain from the first draw to the last.
void require_diagonal(const arma::mat &m, const std::string &what)
{
    for (arma::uword j = 0; j < m.n_cols; j++)
    {
        for (arma::uword i = 0; i < m.n_rows; i++)
        {
            if (i != j && m(i, j) != 0.0)
            {
                throw std::invalid_argument(
                    what + " must be diagonal, since only its diagonal is ever redrawn, but "
                    "element (" + std::to_string(i + 1) + ", " + std::to_string(j + 1) + ") is " +
                    number(m(i, j)));
            }
        }
    }
}

} // namespace

namespace bayests
{
namespace
{

/// Number of periods, after checking that the response divides into them.
/// Every validator starts here: a zero or ragged `y` turns the first reshape
/// into a division by zero or a silently misaligned sample.
arma::uword checked_periods(const VarSpec &spec, const TrainData &train)
{
    spec.validate();

    const arma::uword k = static_cast<arma::uword>(spec.k);

    if (train.y.n_elem == 0)
    {
        throw std::invalid_argument("no training observations");
    }
    if (train.y.n_elem % k != 0)
    {
        throw std::invalid_argument("training observations (" + std::to_string(train.y.n_elem) +
                                    ") are not a multiple of k (" + std::to_string(k) + ")");
    }
    return train.y.n_elem / k;
}

void require_stacked_regressors(const TrainData &train, arma::uword tt, arma::uword k)
{
    if (train.z.n_rows != tt * k)
    {
        throw std::invalid_argument("z must have " + std::to_string(tt * k) +
                                    " rows to match the stacked response, got " +
                                    std::to_string(train.z.n_rows));
    }
}

/// The values of a gamma prior. Zero is allowed for either: an improper prior
/// that the sample makes proper. A negative or non-finite one is not a gamma
/// prior at all.
void require_gamma_values(const GammaPrior &prior, const std::string &what)
{
    require_above(prior.shape, 0.0, false, "gamma prior shape of " + what);
    require_above(prior.rate, 0.0, false, "gamma prior rate of " + what);
}

/// The checks a selection block needs whichever coefficient vector it applies
/// to. `n` is the length of that vector, and the labels name it so the message
/// says whether it was the coefficients or the covariance block that was wrong.
void validate_varsel(const VarSelPrior &prior, const arma::vec &initial_lambda,
                     arma::uword n, VarSelection scheme, const char *block)
{
    const std::string what(block);

    require_length(prior.inprior, n, (what + " prior inclusion probabilities").c_str());
    require_length(initial_lambda, n, (what + " initial inclusion indicators").c_str());
    require_probabilities(prior.inprior, what + " prior inclusion probabilities");

    if (prior.include.n_elem == 0)
    {
        throw std::invalid_argument("variable selection is enabled for " + what +
                                    " but no positions were marked for selection");
    }
    if (prior.include.max() >= n)
    {
        throw std::invalid_argument(what + " variable selection position " +
                                    std::to_string(prior.include.max() + 1) +
                                    " is out of range for " + std::to_string(n) + " elements");
    }
    if (scheme == VarSelection::ssvs)
    {
        require_length(prior.ssvs.tau0, n, (what + " SSVS tau0").c_str());
        require_length(prior.ssvs.tau1, n, (what + " SSVS tau1").c_str());
        require_above(prior.ssvs.tau0, 0.0, true, what + " SSVS tau0");
        require_above(prior.ssvs.tau1, 0.0, true, what + " SSVS tau1");

        // The spike has to be the narrower component. Swapped, the sampler runs
        // exactly as before and every indicator it writes means the opposite of
        // what the file says it means: lambda = 1 would be the tight component.
        // Only the selected positions are read, so only those are held to it.
        for (const arma::uword pos : prior.include)
        {
            if (!(prior.ssvs.tau0(pos) < prior.ssvs.tau1(pos)))
            {
                throw std::invalid_argument(
                    what + " SSVS tau0 must be smaller than tau1 at every selected position, "
                    "the spike being the excluded component and the slab the included one, "
                    "but at position " + std::to_string(pos + 1) + " tau0 is " +
                    number(prior.ssvs.tau0(pos)) + " and tau1 " + number(prior.ssvs.tau1(pos)));
            }
        }
    }
}

/// What SSVS assumes of the normal prior it replaces at the selected positions,
/// and which the sweep cannot check for itself. George, Sun and Ni (2008,
/// eq. 12) centre both mixture components at zero with R = I, and that is what
/// ssvs_sweep() scores against: N(0, tau0^2) against N(0, tau1^2), one position
/// at a time. The coefficient draw, though, reads the prior as the file gives
/// it -- `mu`, and every off-diagonal of `v_inv` -- so where the file departs
/// from the paper the two halves of the Gibbs step are the conditionals of two
/// different models and the chain targets neither. Nothing fails; the
/// inclusion probabilities are simply not the posterior of anything.
///
/// Called after validate_normal_block(), which has checked the sizes.
void validate_ssvs_normal_prior(const VarSelPrior &prior, const NormalPrior &normal,
                                const char *block)
{
    const std::string what(block);

    for (const arma::uword pos : prior.include)
    {
        if (normal.mu(pos) != 0.0)
        {
            throw std::invalid_argument(
                "SSVS centres both components of its prior at zero, so the prior mean of " +
                what + " must be zero at every selected position, but at position " +
                std::to_string(pos + 1) + " it is " + number(normal.mu(pos)) +
                ". A non-zero mean is a prior to shrink towards, not a coefficient to select: "
                "leave that position out of `include`, or use bvs");
        }

        for (arma::uword other = 0; other < normal.v_inv.n_cols; other++)
        {
            if (other != pos && normal.v_inv(pos, other) != 0.0)
            {
                throw std::invalid_argument(
                    "SSVS draws each inclusion indicator from its own coefficient alone, which "
                    "holds only when the prior makes the selected coefficients independent of "
                    "everything else, but the prior precision of " + what + " couples position " +
                    std::to_string(pos + 1) + " to position " + std::to_string(other + 1) +
                    " with " + number(normal.v_inv(pos, other)) +
                    ". Make the rows and columns of the selected positions zero off the diagonal");
            }
        }
    }
}

/// The normal prior plus starting value that both the coefficient block and the
/// covariance block of the constant-coefficient models carry.
void validate_normal_block(const NormalPrior &prior, const arma::vec &initial,
                           arma::uword n, const char *block)
{
    const std::string what(block);
    require_length(prior.mu, n, ("prior mean of " + what).c_str());
    require_square(prior.v_inv, n, ("prior precision of " + what).c_str());
    require_symmetric(prior.v_inv, ("prior precision of " + what).c_str());
    require_length(initial, n, ("initial value of " + what).c_str());
}

/// The random walk state equation plus starting values that every time-varying
/// block carries: a path, the precision of its innovations, the state before
/// the sample, and the prior on both halves.
///
/// `noun` names the thing that drifts and `name` names the vector it is stored
/// in -- ("coefficient", "a") and ("psi", "psi") are the two in use. Two labels
/// rather than one because the messages read better that way and because these
/// are the exact strings the models have always produced.
/// The prior on how far a random walk moves, under whichever parameterisation
/// the file chose: an inverse gamma on the variance, or -- with `omega_v` set --
/// a normal on the signed standard deviation. Both at once is refused rather
/// than resolved, since either reading would ignore half of what the file says.
void validate_state_variance_prior(const RandomWalkPrior &prior, arma::uword n,
                                   const std::string &thing)
{
    if (prior.noncentred())
    {
        if (prior.sigma.shape.n_elem > 0 || prior.sigma.rate.n_elem > 0)
        {
            throw std::invalid_argument(
                "the " + thing + " innovations have both an inverse gamma prior on their "
                "variance (shape, rate) and a normal prior on their standard deviation (omega_v); "
                "give one: omega_v for the non-centred parameterisation, shape and rate for the "
                "centred one");
        }
        require_length(prior.omega_v, n,
                       ("prior variance of the standard deviation of the " + thing + " innovations")
                           .c_str());
        require_above(prior.omega_v, 0.0, true,
                      "prior variance of the standard deviation of the " + thing + " innovations");
        return;
    }

    require_length(prior.sigma.shape, n, ("prior shape of the " + thing + " innovations").c_str());
    require_length(prior.sigma.rate, n, ("prior rate of the " + thing + " innovations").c_str());
    require_gamma_values(prior.sigma, "the " + thing + " innovations");
}

void validate_tvp_block(const RandomWalkPrior &prior, const arma::mat &path,
                        const arma::mat &sigma_inv, const arma::vec &init, arma::uword n,
                        arma::uword tt, const char *noun, const char *name)
{
    const std::string thing(noun);
    const std::string vec(name);

    require_shape(path, n, tt, ("initial " + thing + " path").c_str());
    require_square(sigma_inv, n, ("initial precision of the " + thing + " innovations").c_str());
    require_diagonal(sigma_inv, "initial precision of the " + thing + " innovations");
    require_length(init, n, ("initial value of " + vec + " before the sample").c_str());

    // The non-centred chain starts its standard deviations at the square root
    // of the variance this inverts, so a zero would start it at infinity.
    if (prior.noncentred())
    {
        require_above(arma::vec(sigma_inv.diag()), 0.0, true,
                      "initial precision of the " + thing + " innovations");
    }

    validate_state_variance_prior(prior, n, thing);
    require_length(prior.initial_state.mu, n, ("prior mean of " + vec + " before the sample").c_str());
    require_square(prior.initial_state.v_inv, n,
                   ("prior precision of " + vec + " before the sample").c_str());
    require_symmetric(prior.initial_state.v_inv,
                      ("prior precision of " + vec + " before the sample").c_str());
    require_positive_definite(prior.initial_state.v_inv,
                              ("prior precision of " + vec + " before the sample").c_str());
}

/// A structural model's contemporaneous coefficients are identified only
/// against a diagonal error covariance.
///
/// The data determine the reduced form and nothing else: the coefficients
/// A_0^-1 A_i, and the reduced-form error covariance
/// Omega = A_0^-1 Sigma A_0^-T, which has k(k+1)/2 free elements. A_0 is unit
/// lower triangular and so contributes k(k-1)/2 of its own. Leave Sigma
/// unrestricted -- another k(k+1)/2 -- and the structural side carries k^2
/// parameters mapping onto k(k+1)/2, so a k(k-1)/2-dimensional set of
/// (A_0, Sigma) pairs produces exactly the same Omega and the likelihood is flat
/// along it. Make Sigma diagonal and the count is k(k-1)/2 + k = k(k+1)/2
/// exactly: A_0 and diag(Sigma) are then the LDL factor of Omega, which is
/// unique, and the model is the recursive SVAR.
///
/// Two things in these files leave Sigma unrestricted, and the second is the one
/// that is easy to miss. A Wishart prior does, obviously. So does a covariance
/// block: Sigma^-1 = Psi' Omega^-1 Psi with Psi unit lower triangular and Omega
/// diagonal is a full Sigma, and Psi is then a second contemporaneous matrix
/// doing the same job as A_0 -- only their composition is pinned down.
///
/// Rejected rather than warned about. Inference would still be coherent under a
/// proper prior, and everything these models report is a function of the reduced
/// form alone -- forecasts and the pointwise log likelihood are invariant to
/// where the chain sits on the ridge, and would be correct. But the reason to
/// set the flag at all is to read A_0, and a draw of it here is the prior plus
/// whatever the sampler last wandered onto. That is the failure this codebase
/// refuses elsewhere: output that looks like output.
///
/// `sigma_is_unrestricted` is the model's own answer, and `sigma_source` names
/// what makes it so.
void require_identified_structural(const VarSpec &spec, bool sigma_is_unrestricted,
                                   const char *sigma_source)
{
    if (spec.n_structural() == 0 || !sigma_is_unrestricted)
    {
        return;
    }

    const int free_a0 = spec.n_structural();
    const int free_sigma = spec.k * (spec.k + 1) / 2;

    throw std::invalid_argument(
        "a structural model is not identified alongside " + std::string(sigma_source) +
        ": A_0 contributes " + std::to_string(free_a0) +
        " free elements and an unrestricted error covariance another " +
        std::to_string(free_sigma) + ", against the " + std::to_string(free_sigma) +
        " the reduced-form error covariance determines -- " + std::to_string(free_a0) +
        " more than the data can separate. Estimate the contemporaneous coefficients against a "
        "diagonal covariance instead (the gamma or sv error specification without a covariance "
        "block), or drop them");
}

/// A VEC's `z` has to have the columns the model dimensions describe.
///
/// Counted rather than trusted, because the two arrive from different places in
/// a file -- z from /data/train, the dimensions from /model -- and a spec that
/// disagrees with its data fails nowhere on its own: the sampler sizes
/// everything off z and runs to completion on a model that is not the one the
/// dimensions name. That is the worst kind of wrong, because the output looks
/// like output. It matters twice over for a VEC, whose loading columns are
/// addressed by position and rewritten in place on every draw.
void validate_vec_columns(const VarSpec &spec, arma::uword n_a)
{
    const arma::uword expected = static_cast<arma::uword>(spec.nparams_per_period_vec());
    if (n_a != expected)
    {
        throw std::invalid_argument(
            "z has " + std::to_string(n_a) + " columns but the model dimensions describe " +
            std::to_string(expected) + " coefficients (k*rank + k*(k*(p-1) + m*s + n) with "
            "k=" + std::to_string(spec.k) + ", p=" + std::to_string(spec.p) +
            ", m=" + std::to_string(spec.m) + ", s=" + std::to_string(spec.s) +
            ", n=" + std::to_string(spec.n) + ", rank=" + std::to_string(spec.rank) + ")");
    }
}

/// Selection may not reach a VEC's loadings, in either scheme.
///
/// Excluding one is a change in the rank of Pi, which nothing downstream models.
/// The constant VECs have a second reason: the loadings' prior precision is
/// rebuilt from the cointegration space prior at the top of every draw -- the
/// same matrix SSVS moves between spike and slab, and the same block whose
/// regressors BVS masks -- so whichever writes last wins and neither gets what
/// it meant. In the time-varying VECs the clash is with the beta block, which
/// rewrites those regressors instead. bvartools' .bvectvpalg does apply BVS to
/// the whole of `a`; none of these do.
void validate_vec_varsel_scope(const VarSpec &spec, const VarSelPrior &prior, bool use_beta)
{
    if (use_beta && prior.include.n_elem > 0 &&
        prior.include.min() < static_cast<arma::uword>(spec.n_alpha()))
    {
        throw std::invalid_argument(
            "variable selection cannot be applied to the " + std::to_string(spec.n_alpha()) +
            " loading coefficients at the front of a VEC's a; restrict the selected positions to "
            "the coefficients after them");
    }
}

/// The error correction regressors, which the samplers transpose to k_beta x tt
/// and read one column of per period -- so a `w` of the wrong width fails inside
/// a Kronecker product rather than here.
void validate_vec_w(const VarSpec &spec, const TrainData &train, arma::uword tt)
{
    require_shape(train.w, tt, static_cast<arma::uword>(spec.k_beta),
                  "error correction regressors w");
}

/// A cointegration relation without regressors is not a model any of these
/// samplers can express, and left alone it is silent rather than wrong: the beta
/// block is drawn inside the coefficient block -- beta's own regressors are
/// built from the loadings -- so with no `z` it never runs. Rank r means k*r
/// loading columns, so this cannot happen to a well-formed file; it is a spec
/// and a `z` that disagree, which is what validate_vec_columns() catches when
/// there is a `z` to count.
void require_vec_regressors(const VarSpec &spec, bool use_a)
{
    if (!use_a)
    {
        throw std::invalid_argument(
            "the model has a cointegration relation of rank " + std::to_string(spec.rank) +
            " but no regressors; a VEC of positive rank carries " +
            std::to_string(spec.n_alpha()) + " loading columns in z");
    }
}

/// The quantile a quantile regression model estimates, and the shared checks
/// every asymmetric Laplace model makes before it looks at its own blocks.
///
/// Three refusals rather than one. The quantile has to be a proper one: at zero
/// or one the loss has no minimiser and theta is infinite. A covariance block
/// has to be absent, because Psi rotates the equations into each other and the
/// q-th quantile of a combination is not the combination of q-th quantiles --
/// silently ignoring it would leave a file whose name says quantile and whose
/// numbers do not. And the horizon has to be zero, because iterating a one step
/// quantile does not give an h step one; refusing here rather than in forecast()
/// means the file is rejected before a chain is spent on it.
void validate_ald_spec(const VarSpec &spec)
{
    if (!(spec.quantile > 0.0 && spec.quantile < 1.0))
    {
        throw std::invalid_argument(
            "the quantile of an asymmetric Laplace model must lie in (0, 1), got " +
            std::to_string(spec.quantile));
    }

    if (spec.covar)
    {
        throw std::invalid_argument(
            "a covariance block is not available for a quantile regression model: rotating the "
            "equations into each other leaves a residual whose quantile is not the one asked for");
    }

    if (spec.h != 0)
    {
        throw std::invalid_argument(
            "a quantile regression model does not forecast, so its horizon must be zero, got " +
            std::to_string(spec.h) +
            "; the h step quantile is not the quantile of the iterated one step quantiles");
    }

    if (spec.varsel == VarSelection::ssvs)
    {
        throw std::invalid_argument("SSVS is not implemented for a quantile regression model; "
                                    "expected one of none, bvs");
    }
}

/// The latent scales and the scale of the asymmetric Laplace, shared by both
/// quantile models: the same two blocks at the same widths whatever the
/// coefficients do.
void validate_ald_errors(const GammaPrior &u_scale_prior, const arma::mat &w,
                         const arma::vec &u_scale, arma::uword k, arma::uword tt)
{
    require_length(u_scale_prior.shape, k, "prior shape of the asymmetric Laplace scale");
    require_length(u_scale_prior.rate, k, "prior rate of the asymmetric Laplace scale");
    require_gamma_values(u_scale_prior, "the asymmetric Laplace scale");

    require_shape(w, tt, k, "initial latent scales");
    require_length(u_scale, k, "initial scale of the asymmetric Laplace");

    if (w.min() <= 0.0)
    {
        throw std::invalid_argument(
            "every initial latent scale must be positive: it is the variance the first draw of the "
            "coefficients is weighted by, and a zero divides by nothing");
    }

    if (u_scale.min() <= 0.0)
    {
        throw std::invalid_argument("every initial scale of the asymmetric Laplace must be positive");
    }
}

/// rho scales the cointegration state path itself, so a value outside (0, 1]
/// either reverses the sign of the relation from period to period or lets it
/// grow without bound. One is the random walk bvartools' .bvectvpalg uses.
///
/// Where the file puts a prior on rho rather than fixing it, the support has to
/// be an interval of the same (0, 1], and the value the chain starts at has to
/// be inside it: a starting value its own prior gives no weight to is a file
/// that means two different things at once, and the first draw would move it
/// without saying so.
void validate_tvp_coint_rho(const TvpCointSpacePrior &prior)
{
    if (!(prior.rho > 0.0 && prior.rho <= 1.0))
    {
        throw std::invalid_argument(
            "the autoregression of the cointegration state equation (rho) must lie in (0, 1], "
            "got " + std::to_string(prior.rho));
    }

    if (!prior.rho_prior.draw)
    {
        return;
    }

    const double min = prior.rho_prior.min;
    const double max = prior.rho_prior.max;

    if (!(min > 0.0 && max <= 1.0 && min < max))
    {
        throw std::invalid_argument(
            "the prior support of the autoregression of the cointegration state equation must be "
            "an interval within (0, 1], got [" + std::to_string(min) + ", " + std::to_string(max) +
            "]");
    }

    if (prior.rho < min || prior.rho > max)
    {
        throw std::invalid_argument(
            "the starting value of the autoregression of the cointegration state equation (rho) "
            "lies outside its own prior support: got " + std::to_string(prior.rho) + " for [" +
            std::to_string(min) + ", " + std::to_string(max) + "]");
    }
}

} // namespace

/// The selected positions a BVS sweep would be scoring against a prior draw
/// rather than against the data. See the declaration in priors.h for why the
/// diagonal is what to read, and why this is a report and not a refusal.
FlatSelectionPrior flat_selection_prior(const VarSelPrior &prior, const arma::mat &v_inv,
                                        const double variance_threshold)
{
    FlatSelectionPrior report;
    report.selected = prior.include.n_elem;

    const arma::uword side = std::min(v_inv.n_rows, v_inv.n_cols);
    for (const arma::uword pos : prior.include)
    {
        // validate_varsel() refuses a position past the end of the block, but
        // this is a diagnostic and may be called before or instead of it.
        if (pos >= side)
        {
            continue;
        }

        const double precision = v_inv(pos, pos);
        const double variance = precision > 0.0 ? 1.0 / precision
                                                : std::numeric_limits<double>::infinity();
        if (variance < variance_threshold)
        {
            continue;
        }

        report.flat++;
        if (report.flat == 1 || variance > report.worst_variance)
        {
            report.worst_variance = variance;
            report.worst_position = pos;
        }
    }

    return report;
}

std::string flat_selection_message(const FlatSelectionPrior &report, const std::string &block)
{
    const std::string variance = std::isinf(report.worst_variance)
                                     ? std::string("infinity -- no prior precision at all")
                                     : number(report.worst_variance);

    return "bvs is selecting over " + std::to_string(report.selected) + " position(s) of " +
           block + ", and " + std::to_string(report.flat) +
           " of them have a prior too flat to select against: /priors/" + block +
           "/v_inv leaves position " + std::to_string(report.worst_position + 1) +
           " -- one-based, as `include` counts -- a conditional prior variance of " + variance +
           ". An excluded coefficient is drawn from that prior before it is scored against the "
           "data, so the wider that prior, the harder it is for anything to get back in once it "
           "is out. Expect inclusion probabilities pinned near zero that say more about the "
           "prior than about the data. Korobilis (2013) puts the usable range at a prior "
           "variance of 0.25 to 25";
}

} // namespace bayests

namespace bayests
{

void VarNormalWishartInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword nparams = train.nparams();

    require_identified_structural(spec, true, "a Wishart prior on the error precision");

    if (use_a())
    {
        require_stacked_regressors(train, tt, k);
        validate_normal_block(a_prior, initial.a, nparams, "a");

        if (spec.uses_varsel())
        {
            validate_varsel(varsel_prior, initial.a_lambda, nparams, spec.varsel, "a");
        }
        if (spec.varsel == VarSelection::ssvs)
        {
            validate_ssvs_normal_prior(varsel_prior, a_prior, "a");
        }
    }

    if (u_sigma_prior.df <= 0)
    {
        throw std::invalid_argument("Wishart prior degrees of freedom must be positive");
    }
    require_square(u_sigma_prior.scale, k, "Wishart prior scale");
    require_symmetric(u_sigma_prior.scale, "Wishart prior scale");
    require_square(initial.u_sigma_inv, k, "initial error precision");
}

void VarNormalGammaInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword nparams = train.nparams();

    require_identified_structural(spec, use_psi(), "a covariance block");

    if (use_a())
    {
        require_stacked_regressors(train, tt, k);
        validate_normal_block(a_prior, initial.a, nparams, "a");

        if (spec.uses_varsel())
        {
            validate_varsel(a_varsel_prior, initial.a_lambda, nparams, spec.varsel, "a");
        }
        if (spec.varsel == VarSelection::ssvs)
        {
            validate_ssvs_normal_prior(a_varsel_prior, a_prior, "a");
        }
    }

    if (use_psi())
    {
        const arma::uword n_psi = static_cast<arma::uword>(spec.n_psi());
        validate_normal_block(psi_prior, initial.psi, n_psi, "psi");

        if (spec.uses_varsel())
        {
            validate_varsel(psi_varsel_prior, initial.psi_lambda, n_psi, spec.varsel, "psi");
        }
        if (spec.varsel == VarSelection::ssvs)
        {
            validate_ssvs_normal_prior(psi_varsel_prior, psi_prior, "psi");
        }
    }

    require_length(u_sigma_prior.shape, k, "gamma prior shape of the error precision");
    require_length(u_sigma_prior.rate, k, "gamma prior rate of the error precision");
    require_gamma_values(u_sigma_prior, "the error precision");
    require_square(initial.u_sigma_inv, k, "initial error precision");
    require_diagonal(initial.u_sigma_inv, "initial error precision");
}

void VarNormalStochvolInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword nparams = train.nparams();

    require_identified_structural(spec, use_psi(), "a covariance block");

    if (spec.varsel == VarSelection::ssvs)
    {
        throw std::invalid_argument("SSVS is not implemented for a stochastic volatility model; "
                                    "expected one of none, bvs");
    }

    if (use_a())
    {
        require_stacked_regressors(train, tt, k);
        validate_normal_block(a_prior, initial.a, nparams, "a");

        if (spec.uses_varsel())
        {
            validate_varsel(a_varsel_prior, initial.a_lambda, nparams, spec.varsel, "a");
        }
    }

    if (use_psi())
    {
        const arma::uword n_psi = static_cast<arma::uword>(spec.n_psi());
        validate_normal_block(psi_prior, initial.psi, n_psi, "psi");

        if (spec.uses_varsel())
        {
            validate_varsel(psi_varsel_prior, initial.psi_lambda, n_psi, spec.varsel, "psi");
        }
    }

    require_length(u_sigma_prior.offset, k, "log-volatility offset");
    require_length(u_sigma_prior.state.sigma.shape, k, "prior shape of the log-volatility variance");
    require_length(u_sigma_prior.state.sigma.rate, k, "prior rate of the log-volatility variance");
    require_length(u_sigma_prior.state.initial_state.mu, k, "prior mean of the initial log-volatility");
    require_square(u_sigma_prior.state.initial_state.v_inv, k,
                   "prior precision of the initial log-volatility");

    require_above(u_sigma_prior.offset, 0.0, true, "log-volatility offset");
    require_gamma_values(u_sigma_prior.state.sigma, "the log-volatility variance");
    require_symmetric(u_sigma_prior.state.initial_state.v_inv,
                      "prior precision of the initial log-volatility");

    require_shape(initial.h, tt, k, "initial log-volatility");
    require_length(initial.h_init, k, "initial value of the log-volatility before the sample");
    require_length(initial.h_sigma, k, "initial variance of the log-volatility innovations");
    require_above(initial.h_sigma, 0.0, true, "initial variance of the log-volatility innovations");

    // The random walk differences h against its own lag, so a single period
    // leaves nothing to difference.
    if (tt < 2)
    {
        throw std::invalid_argument("a stochastic volatility model needs at least two periods");
    }
}

void VarNormalAldInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword nparams = train.nparams();

    validate_ald_spec(spec);

    // Sigma is diagonal here and always will be, so a structural model is
    // identified by that alone -- which is what this call says with `false`.
    require_identified_structural(spec, false, "a covariance block");

    if (use_a())
    {
        require_stacked_regressors(train, tt, k);
        validate_normal_block(a_prior, initial.a, nparams, "a");

        if (spec.uses_varsel())
        {
            validate_varsel(a_varsel_prior, initial.a_lambda, nparams, spec.varsel, "a");
        }
    }

    validate_ald_errors(u_scale_prior, initial.w, initial.u_scale, k, tt);
}

void VarTvpAldInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword nparams = train.nparams();

    validate_ald_spec(spec);
    require_identified_structural(spec, false, "a covariance block");

    if (use_a())
    {
        require_stacked_regressors(train, tt, k);
        validate_tvp_block(a_prior, initial.a, initial.a_sigma_inv, initial.a_init, nparams, tt,
                           "coefficient", "a");

        if (spec.uses_varsel())
        {
            validate_varsel(a_varsel_prior, initial.a_lambda, nparams, spec.varsel, "a");
        }
    }

    validate_ald_errors(u_scale_prior, initial.w, initial.u_scale, k, tt);

    // The random walk differences a against its own lag, so a single period
    // leaves nothing to difference.
    if (tt < 2)
    {
        throw std::invalid_argument("a time-varying parameter model needs at least two periods");
    }
}

void VarTvpGammaInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword nparams = train.nparams();

    require_identified_structural(spec, use_psi(), "a covariance block");

    if (spec.varsel == VarSelection::ssvs || psi_varsel == VarSelection::ssvs)
    {
        throw std::invalid_argument("SSVS is not implemented for a time-varying parameter model; "
                                    "expected one of none, bvs");
    }

    // The state equation lags the path against itself, and the smoother is
    // handed columns 1..tt of a t+1 wide result.
    if (tt < 2)
    {
        throw std::invalid_argument("a time-varying parameter model needs at least two periods");
    }

    if (use_a())
    {
        require_stacked_regressors(train, tt, k);
        validate_tvp_block(a_prior, initial.a, initial.a_sigma_inv, initial.a_init, nparams, tt,
                           "coefficient", "a");

        if (spec.uses_varsel())
        {
            validate_varsel(a_varsel_prior, initial.a_lambda, nparams, spec.varsel, "a");
        }
    }

    if (use_psi())
    {
        const arma::uword n_psi = static_cast<arma::uword>(spec.n_psi());
        validate_tvp_block(psi_prior, initial.psi, initial.psi_sigma_inv, initial.psi_init, n_psi,
                           tt, "psi", "psi");

        if (uses_psi_varsel())
        {
            validate_varsel(psi_varsel_prior, initial.psi_lambda, n_psi, psi_varsel, "psi");
        }
    }

    require_length(u_sigma_prior.shape, k, "gamma prior shape of the error precision");
    require_length(u_sigma_prior.rate, k, "gamma prior rate of the error precision");
    require_gamma_values(u_sigma_prior, "the error precision");
    require_square(initial.u_omega_inv, k, "initial error precision");
    require_diagonal(initial.u_omega_inv, "initial error precision");
}

void VarTvpWishartInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword nparams = train.nparams();

    require_identified_structural(spec, true, "a Wishart prior on the error precision");

    if (spec.varsel == VarSelection::ssvs)
    {
        throw std::invalid_argument("SSVS is not implemented for a time-varying parameter model; "
                                    "expected one of none, bvs");
    }

    // The state equation lags the path against itself, and the smoother is
    // handed columns 1..tt of a t+1 wide result.
    if (tt < 2)
    {
        throw std::invalid_argument("a time-varying parameter model needs at least two periods");
    }

    if (use_a())
    {
        require_stacked_regressors(train, tt, k);
        validate_tvp_block(a_prior, initial.a, initial.a_sigma_inv, initial.a_init, nparams, tt,
                           "coefficient", "a");

        if (spec.uses_varsel())
        {
            validate_varsel(a_varsel_prior, initial.a_lambda, nparams, spec.varsel, "a");
        }
    }

    // Unlike VarTvpGamma this model has no psi block: the covariance is carried
    // by the Wishart precision alone, so there is nothing here to match against
    // spec.n_psi().
    if (u_sigma_prior.df <= 0)
    {
        throw std::invalid_argument("Wishart prior degrees of freedom must be positive");
    }
    require_square(u_sigma_prior.scale, k, "Wishart prior scale");
    require_symmetric(u_sigma_prior.scale, "Wishart prior scale");
    require_square(initial.u_sigma_inv, k, "initial error precision");
}


void VarTvpStochvolInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword nparams = train.nparams();

    require_identified_structural(spec, use_psi(), "a covariance block");

    if (spec.varsel == VarSelection::ssvs || psi_varsel == VarSelection::ssvs)
    {
        throw std::invalid_argument("SSVS is not implemented for a time-varying parameter model "
                                    "with stochastic volatility; expected one of none, bvs");
    }

    // The state equation lags the path against itself, and the smoother is
    // handed columns 1..tt of a t+1 wide result.
    if (tt < 2)
    {
        throw std::invalid_argument("a time-varying parameter model needs at least two periods");
    }

    if (use_a())
    {
        require_stacked_regressors(train, tt, k);
        validate_tvp_block(a_prior, initial.a, initial.a_sigma_inv, initial.a_init, nparams, tt,
                           "coefficient", "a");

        if (spec.uses_varsel())
        {
            validate_varsel(a_varsel_prior, initial.a_lambda, nparams, spec.varsel, "a");
        }
    }

    if (use_psi())
    {
        const arma::uword n_psi = static_cast<arma::uword>(spec.n_psi());
        validate_tvp_block(psi_prior, initial.psi, initial.psi_sigma_inv, initial.psi_init, n_psi,
                           tt, "psi", "psi");

        if (uses_psi_varsel())
        {
            validate_varsel(psi_varsel_prior, initial.psi_lambda, n_psi, psi_varsel, "psi");
        }
    }

    require_length(u_sigma_prior.offset, k, "offset of the log-volatility measurement equation");
    validate_state_variance_prior(u_sigma_prior.state, k, "log-volatility");
    require_length(u_sigma_prior.state.initial_state.mu, k, "prior mean of the log-volatility before the sample");
    require_square(u_sigma_prior.state.initial_state.v_inv, k, "prior precision of the log-volatility before the sample");

    require_above(u_sigma_prior.offset, 0.0, true,
                  "offset of the log-volatility measurement equation");
    require_symmetric(u_sigma_prior.state.initial_state.v_inv,
                      "prior precision of the log-volatility before the sample");

    require_length(initial.h_sigma, k, "initial variance of the log-volatility innovations");
    require_above(initial.h_sigma, 0.0, true, "initial variance of the log-volatility innovations");
    require_length(initial.h_init, k, "initial log-volatility before the sample");
    require_shape(initial.h, tt, k, "initial log-volatility path");
}

namespace
{

/// The half of a VEC's validation that does not depend on how its coefficients
/// or its errors are modelled: the shape of `z` and the scope of any selection.
void validate_vec_coefficients(const VarSpec &spec, const TrainData &train,
                               const VarSelPrior &varsel_prior, const arma::vec &initial_lambda,
                               arma::uword tt, arma::uword k, arma::uword n_a, bool use_beta)
{
    require_stacked_regressors(train, tt, k);
    validate_vec_columns(spec, n_a);

    if (spec.uses_varsel())
    {
        validate_varsel(varsel_prior, initial_lambda, n_a, spec.varsel, "a");
        validate_vec_varsel_scope(spec, varsel_prior, use_beta);
    }
}

/// The constant cointegration block: a starting value, the error correction
/// regressors, and the central location of the cointegration space.
void validate_constant_coint_block(const VarSpec &spec, const TrainData &train,
                                   const ConstantCointSpacePrior &prior, const arma::vec &initial,
                                   arma::uword tt, bool use_a)
{
    // Checked here as well as in the time-varying block, and for a sharper
    // reason: the constant VECs read the loadings straight out of the front of
    // `a`, so without regressors that subvector is taken out of an empty
    // vector rather than merely left undrawn.
    require_vec_regressors(spec, use_a);

    require_length(initial, static_cast<arma::uword>(spec.n_beta()), "initial beta");
    validate_vec_w(spec, train, tt);

    // k_beta square, not n_beta square: P_tau^-1 is the prior's central location
    // for the cointegration *space*, so it is indexed by the rows of beta and
    // carries no rank dimension. It enters the draws only through
    // kron(., P_tau^-1) and beta' P_tau^-1 beta, both of which want k_beta.
    // Demanding n_beta here happened to pass at rank one, where the two
    // coincide, and rejected every well-formed file above it.
    require_square(prior.p_tau_inv, static_cast<arma::uword>(spec.k_beta),
                   "prior precision of the cointegration space");

    // And symmetric, which both of those assume without checking.
    // kron(., P_tau^-1) goes into the posterior precision of beta, which the draw
    // reads through its upper triangle alone. beta' P_tau^-1 beta goes into the
    // prior precision of the loadings, whose diagonal sees only the symmetric
    // part of P_tau^-1 and whose off-diagonal sees the whole of it. So an
    // asymmetric one was not even the same prior in the two blocks.
    require_symmetric(prior.p_tau_inv, "prior precision of the cointegration space");

    // v scales the prior precision of the loadings, so a negative one is a
    // negative precision -- a prior that is not a density, which the draws
    // nonetheless run on. Zero is the flat prior on alpha and the uniform one on
    // the space, and is fine.
    if (!std::isfinite(prior.v_inv) || prior.v_inv < 0.0)
    {
        throw std::invalid_argument("the shrinkage v of the cointegration space prior, "
                                    "/priors/beta/v_inv, must be finite and at least zero, got " +
                                    number(prior.v_inv));
    }

    // P_tau^-1 has to be the inverse of something: the matrix angular central
    // Gaussian is defined for a positive definite P_tau, and tau = 0 -- the
    // dogmatic prior that puts the space exactly on sp(H) -- has no inverse to
    // give. Only read when v > 0; at v = 0 the space is uniform whatever it is.
    if (prior.v_inv > 0.0)
    {
        require_positive_definite(prior.p_tau_inv, "prior precision of the cointegration space");
    }
}

/// What Koop, Leon-Gonzalez and Strachan's (2010) collapsed sampler assumes of
/// the loadings' prior and cannot check for itself. Their Proposition 1 -- that
/// B is normal given A -- and the term the prior adds to Sigma's posterior in
/// their eq. (8) both rest on alpha | beta ~ N(0, v^-1 (beta' P^-1 beta)^-1
/// kron G), independent of every other coefficient. The samplers rebuild that
/// prior's precision block every draw, so the file's values inside it are never
/// read; what the file does get to say is the mean of alpha and its prior
/// correlation with the rest of `a`, and either one non-zero makes the three
/// Gibbs blocks the conditionals of three different priors. Nothing fails. The
/// chain just targets none of them.
void validate_coint_loadings_prior(const VarSpec &spec, const NormalPrior &a_prior)
{
    const arma::uword n_alpha = static_cast<arma::uword>(spec.n_alpha());
    const arma::uword n_a = a_prior.mu.n_elem;

    for (arma::uword i = 0; i < n_alpha; i++)
    {
        if (a_prior.mu(i) != 0.0)
        {
            throw std::invalid_argument(
                "the cointegration space prior centres the loadings at zero, so the first " +
                std::to_string(n_alpha) + " elements of /priors/a/mu -- the loadings alpha -- "
                "must be zero, but element " + std::to_string(i + 1) + " is " +
                number(a_prior.mu(i)));
        }
        for (arma::uword j = n_alpha; j < n_a; j++)
        {
            if (a_prior.v_inv(i, j) != 0.0)
            {
                throw std::invalid_argument(
                    "the cointegration space prior makes the loadings independent of the other "
                    "coefficients, so /priors/a/v_inv must be zero between the first " +
                    std::to_string(n_alpha) + " positions and the rest, but it couples position " +
                    std::to_string(i + 1) + " to position " + std::to_string(j + 1) + " with " +
                    number(a_prior.v_inv(i, j)));
            }
        }
    }
}

/// G, the matrix the loadings' prior is scaled by, is the error covariance in
/// every constant VEC but VecNormalStochvol, where there is a different one in
/// every period and the file supplies G instead. A G given to one of the others
/// would be ignored, so it is refused.
void refuse_coint_g(const ConstantCointSpacePrior &prior)
{
    if (!prior.g_inv.is_empty())
    {
        throw std::invalid_argument(
            "/priors/beta/g_inv is read by VecNormalStochvol alone: every other constant VEC "
            "scales the loadings' prior by its error covariance, which is what Koop, "
            "Leon-Gonzalez and Strachan's eq. (8) needs, so a G given here would never be used");
    }
}

/// The time-varying cointegration block: a path, where it starts, the error
/// correction regressors, and the state equation's autoregression. No state
/// variance among them -- see TvpCointSpacePrior for why it is fixed.
void validate_tvp_coint_block(const VarSpec &spec, const TrainData &train,
                              const TvpCointSpacePrior &prior, const arma::mat &initial_path,
                              const arma::vec &initial_state, arma::uword tt, bool use_a)
{
    const arma::uword n_beta = static_cast<arma::uword>(spec.n_beta());

    require_vec_regressors(spec, use_a);
    require_shape(initial_path, n_beta, tt, "initial cointegration path");
    require_length(initial_state, n_beta, "initial value of beta before the sample");
    validate_vec_w(spec, train, tt);

    require_length(prior.initial_state.mu, n_beta, "prior mean of beta before the sample");
    require_square(prior.initial_state.v_inv, n_beta, "prior precision of beta before the sample");

    // Symmetric for the reason P_tau^-1 is in the constant block: the draw of the
    // state before the sample factorises v_inv + rho^2 P'P through its upper
    // triangle, and multiplies the prior mean by the whole of v_inv.
    require_symmetric(prior.initial_state.v_inv, "prior precision of beta before the sample");

    validate_tvp_coint_rho(prior);

    // P_tau is the transition with rho taken out, and it has to leave the state
    // equation what Koop, Leon-Gonzalez and Strachan's is: a pull towards sp(H)
    // that never pushes away from it and never flips a direction's sign from one
    // period to the next. Symmetric with eigenvalues in [0, 1] is exactly that --
    // one along H, below one off it, a tau per direction.
    if (!prior.p_tau.is_empty())
    {
        require_square(prior.p_tau, static_cast<arma::uword>(spec.k_beta),
                       "transition P_tau of the cointegration state equation");
        require_symmetric(prior.p_tau, "transition P_tau of the cointegration state equation");

        const arma::vec eigval = arma::eig_sym(arma::symmatu(prior.p_tau));
        if (eigval.min() < -1e-10 || eigval.max() > 1.0 + 1e-10)
        {
            throw std::invalid_argument(
                "the eigenvalues of the transition P_tau of the cointegration state equation must "
                "lie in [0, 1], got [" + std::to_string(eigval.min()) + ", " +
                std::to_string(eigval.max()) + "]");
        }
    }
}

/// The two periods a random walk needs to be differenced against itself, and the
/// SSVS rejection every time-varying model shares.
void validate_tvp_preconditions(VarSelection varsel, VarSelection psi_varsel, arma::uword tt)
{
    if (varsel == VarSelection::ssvs || psi_varsel == VarSelection::ssvs)
    {
        throw std::invalid_argument("SSVS is not implemented for a time-varying parameter model; "
                                    "expected one of none, bvs");
    }

    // The state equation lags the path against itself, and the smoother is
    // handed columns 1..tt of a t+1 wide result.
    if (tt < 2)
    {
        throw std::invalid_argument("a time-varying parameter model needs at least two periods");
    }
}

/// The stochastic volatility block's priors and starting values.
void validate_stochvol_block(const StochvolPrior &prior, const arma::vec &h_sigma,
                             const arma::vec &h_init, const arma::mat &h, arma::uword k,
                             arma::uword tt)
{
    require_length(prior.offset, k, "offset of the log-volatility measurement equation");
    validate_state_variance_prior(prior.state, k, "log-volatility");
    require_length(prior.state.initial_state.mu, k,
                   "prior mean of the log-volatility before the sample");
    require_square(prior.state.initial_state.v_inv, k,
                   "prior precision of the log-volatility before the sample");

    require_above(prior.offset, 0.0, true, "offset of the log-volatility measurement equation");
    require_symmetric(prior.state.initial_state.v_inv,
                      "prior precision of the log-volatility before the sample");

    require_length(h_sigma, k, "initial variance of the log-volatility innovations");
    require_above(h_sigma, 0.0, true, "initial variance of the log-volatility innovations");
    require_length(h_init, k, "initial log-volatility before the sample");
    require_shape(h, tt, k, "initial log-volatility path");
}

/// The Wishart prior on the error precision, and the value the chain starts it
/// at.
void validate_wishart_block(const WishartPrior &prior, const arma::mat &initial, arma::uword k)
{
    if (prior.df <= 0)
    {
        throw std::invalid_argument("Wishart prior degrees of freedom must be positive");
    }
    require_square(prior.scale, k, "Wishart prior scale");
    require_symmetric(prior.scale, "Wishart prior scale");
    require_square(initial, k, "initial error precision");
}

} // namespace

void VecNormalWishartInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n_a = train.nparams();

    require_identified_structural(spec, true, "a Wishart prior on the error precision");

    if (use_a())
    {
        validate_vec_coefficients(spec, train, varsel_prior, initial.a_lambda, tt, k, n_a,
                                  use_beta());
        validate_normal_block(a_prior, initial.a, n_a, "a");
        if (spec.varsel == VarSelection::ssvs)
        {
            validate_ssvs_normal_prior(varsel_prior, a_prior, "a");
        }
    }

    if (use_beta())
    {
        validate_constant_coint_block(spec, train, beta_prior, initial.beta, tt, use_a());
        validate_coint_loadings_prior(spec, a_prior);
        refuse_coint_g(beta_prior);
    }

    validate_wishart_block(u_sigma_prior, initial.u_sigma_inv, k);
}

namespace
{

/// The checks both dynamic factor models share: what makes the dimensions a
/// factor model at all, and the two coefficient blocks. Only the error
/// specification differs between them -- gamma priors on two precisions against
/// two stochastic volatility blocks -- so only that is left to the callers.
void validate_dfm_dimensions(const VarSpec &spec, arma::uword tt)
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword n = static_cast<arma::uword>(spec.n_factors);

    if (spec.n_factors <= 0)
    {
        throw std::invalid_argument("a dynamic factor model must have at least one factor "
                                    "(n_factors), got " + std::to_string(spec.n_factors));
    }

    // More factors than series and the loading matrix has no identifying block
    // to build: the rotation that Lambda's unit lower triangle pins down needs
    // one row per factor to pin it with.
    if (n > k)
    {
        throw std::invalid_argument(
            "a dynamic factor model cannot have more factors (" + std::to_string(spec.n_factors) +
            ") than observed series (" + std::to_string(spec.k) + ")");
    }

    // The transition regresses on p lags of the factors, and the factors are the
    // sample: there has to be a sample left after the longest of them.
    if (static_cast<arma::uword>(spec.p) >= tt)
    {
        throw std::invalid_argument(
            "a factor transition of order " + std::to_string(spec.p) +
            " needs more than that many periods, and the sample has " + std::to_string(tt));
    }

    if (spec.uses_varsel())
    {
        throw std::invalid_argument("variable selection is not implemented for a dynamic factor "
                                    "model; expected varsel none");
    }

    if (spec.structural)
    {
        throw std::invalid_argument("a dynamic factor model has no contemporaneous coefficients "
                                    "to identify; expected structural false");
    }
}

/// The two coefficient blocks of a dynamic factor model whose coefficients do
/// not move: a normal prior over the free loadings and one over the transition.
void validate_dfm_shape(const VarSpec &spec, arma::uword tt, const NormalPrior &lambda_prior,
                        const arma::vec &initial_lambda, const NormalPrior &a_prior,
                        const arma::vec &initial_a, bool use_a)
{
    validate_dfm_dimensions(spec, tt);

    // The free loadings, in the row-major order the sampler draws them.
    const arma::uword n_lambda = static_cast<arma::uword>(spec.n_lambda());
    validate_normal_block(lambda_prior, initial_lambda, n_lambda, "lambda");

    if (use_a)
    {
        validate_normal_block(a_prior, initial_a, static_cast<arma::uword>(spec.n_factor_a()), "a");
    }
}

/// The two gamma-distributed error precisions a dynamic factor model carries.
/// Both are diagonal, so both arrive as the diagonal rather than as a matrix --
/// and a starting value that is not positive is not one, since the first factor
/// draw inverts it.
void validate_dfm_gamma_errors(const GammaPrior &u_prior, const arma::vec &initial_u,
                               const GammaPrior &v_prior, const arma::vec &initial_v,
                               arma::uword k, arma::uword n)
{
    require_length(u_prior.shape, k, "gamma prior shape of the idiosyncratic precision");
    require_length(u_prior.rate, k, "gamma prior rate of the idiosyncratic precision");
    require_length(initial_u, k, "initial idiosyncratic precision");

    require_length(v_prior.shape, n, "gamma prior shape of the factor innovation precision");
    require_length(v_prior.rate, n, "gamma prior rate of the factor innovation precision");
    require_length(initial_v, n, "initial factor innovation precision");

    require_gamma_values(u_prior, "the idiosyncratic precision");
    require_gamma_values(v_prior, "the factor innovation precision");

    if (initial_u.min() <= 0.0)
    {
        throw std::invalid_argument("every element of the initial idiosyncratic precision must be "
                                    "positive");
    }
    if (initial_v.min() <= 0.0)
    {
        throw std::invalid_argument("every element of the initial factor innovation precision must "
                                    "be positive");
    }
}

/// One of the two stochastic volatility blocks a dynamic factor model carries.
/// They differ in their width and in what a message calls them, and in nothing
/// else, so reading both through one function is what keeps the two from drifting
/// apart. `errors` and `block` are the only two words that differ.
void validate_dfm_stochvol_block(const StochvolPrior &prior, const arma::mat &h,
                                 const arma::vec &h_init, const arma::vec &h_sigma,
                                 arma::uword width, arma::uword tt, const char *errors,
                                 const char *block)
{
    const std::string of_errors(errors);
    const std::string of_block(block);

    require_length(prior.offset, width, ("log-volatility offset of the " + of_errors).c_str());
    require_length(prior.state.sigma.shape, width,
                   ("prior shape of the " + of_block + " log-volatility variance").c_str());
    require_length(prior.state.sigma.rate, width,
                   ("prior rate of the " + of_block + " log-volatility variance").c_str());
    require_length(prior.state.initial_state.mu, width,
                   ("prior mean of the initial " + of_block + " log-volatility").c_str());
    require_square(prior.state.initial_state.v_inv, width,
                   ("prior precision of the initial " + of_block + " log-volatility").c_str());

    require_above(prior.offset, 0.0, true, "log-volatility offset of the " + of_errors);
    require_gamma_values(prior.state.sigma, "the " + of_block + " log-volatility variance");
    require_symmetric(prior.state.initial_state.v_inv,
                      ("prior precision of the initial " + of_block + " log-volatility").c_str());

    require_shape(h, tt, width, ("initial " + of_block + " log-volatility").c_str());
    require_length(h_init, width,
                   ("initial " + of_block + " log-volatility before the sample").c_str());
    require_length(h_sigma, width,
                   ("initial variance of the " + of_block + " log-volatility innovations").c_str());

    // The variance is divided by, once per period, inside the banded draw of the
    // log-volatility path. A zero there is an infinity that only surfaces as a
    // non-finite draw several steps later.
    if (h_sigma.min() <= 0.0)
    {
        throw std::invalid_argument("every element of the initial variance of the " + of_block +
                                    " log-volatility innovations must be positive");
    }
}

} // namespace

void DfmNormalGammaInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n = static_cast<arma::uword>(spec.n_factors);

    validate_dfm_shape(spec, tt, lambda_prior, initial.lambda, a_prior, initial.a, use_a());

    validate_dfm_gamma_errors(u_sigma_prior, initial.u_sigma_inv, v_sigma_prior,
                              initial.v_sigma_inv, k, n);
}

void DfmTvpGammaInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n = static_cast<arma::uword>(spec.n_factors);

    validate_dfm_dimensions(spec, tt);

    // Both state equations difference their path against its own lag, so a
    // single period leaves nothing to difference -- the same floor every
    // time-varying parameter model here has.
    if (tt < 2)
    {
        throw std::invalid_argument("a time-varying parameter model needs at least two periods");
    }

    // A single observed series leaves no free loading to give a state equation
    // to: the whole of Lambda is then the identifying block. Skipped rather than
    // demanded as a set of empty datasets, which is what a file would otherwise
    // have to carry to describe a model with nothing in the block.
    const arma::uword n_lambda = static_cast<arma::uword>(spec.n_lambda());
    if (n_lambda > 0)
    {
        validate_tvp_block(lambda_prior, initial.lambda, initial.lambda_sigma_inv,
                           initial.lambda_init, n_lambda, tt, "loading", "lambda");
    }

    if (use_a())
    {
        validate_tvp_block(a_prior, initial.a, initial.a_sigma_inv, initial.a_init,
                           static_cast<arma::uword>(spec.n_factor_a()), tt, "transition", "a");
    }

    validate_dfm_gamma_errors(u_sigma_prior, initial.u_sigma_inv, v_sigma_prior,
                              initial.v_sigma_inv, k, n);
}

void DfmTvpStochvolInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n = static_cast<arma::uword>(spec.n_factors);

    validate_dfm_dimensions(spec, tt);

    // Four random walks here, and every one of them differences its path against
    // its own lag, so a single period leaves nothing to difference.
    if (tt < 2)
    {
        throw std::invalid_argument("a time-varying parameter model needs at least two periods");
    }

    // A single observed series leaves no free loading to give a state equation
    // to; see DfmTvpGammaInput::validate().
    const arma::uword n_lambda = static_cast<arma::uword>(spec.n_lambda());
    if (n_lambda > 0)
    {
        validate_tvp_block(lambda_prior, initial.lambda, initial.lambda_sigma_inv,
                           initial.lambda_init, n_lambda, tt, "loading", "lambda");
    }

    if (use_a())
    {
        validate_tvp_block(a_prior, initial.a, initial.a_sigma_inv, initial.a_init,
                           static_cast<arma::uword>(spec.n_factor_a()), tt, "transition", "a");
    }

    validate_dfm_stochvol_block(u_sigma_prior, initial.u_h, initial.u_h_init, initial.u_h_sigma,
                                k, tt, "idiosyncratic errors", "idiosyncratic");
    validate_dfm_stochvol_block(v_sigma_prior, initial.v_h, initial.v_h_init, initial.v_h_sigma,
                                n, tt, "factor innovations", "factor-innovation");
}

void DfmNormalStochvolInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n = static_cast<arma::uword>(spec.n_factors);

    validate_dfm_shape(spec, tt, lambda_prior, initial.lambda, a_prior, initial.a, use_a());

    // Both random walks difference their path against its own lag, so a single
    // period leaves nothing to difference -- the same floor every stochastic
    // volatility model here has.
    if (tt < 2)
    {
        throw std::invalid_argument("a stochastic volatility model needs at least two periods");
    }

    // One block per error term, at its own width. That pair of lengths is what a
    // file is most likely to get wrong -- k and n_factors are the two counts a
    // factor model carries -- so the two calls name which is which.
    validate_dfm_stochvol_block(u_sigma_prior, initial.u_h, initial.u_h_init, initial.u_h_sigma,
                                k, tt, "idiosyncratic errors", "idiosyncratic");
    validate_dfm_stochvol_block(v_sigma_prior, initial.v_h, initial.v_h_init, initial.v_h_sigma,
                                n, tt, "factor innovations", "factor-innovation");
}

void FavarNormalWishartInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n_obs = static_cast<arma::uword>(spec.n_obs_factors);
    const arma::uword n_state = static_cast<arma::uword>(spec.n_state());

    // The unobserved half is a factor model and is checked as one: at least one
    // factor, no more of them than series, a sample longer than the transition,
    // and neither variable selection nor a structural block.
    validate_dfm_dimensions(spec, tt);

    if (spec.uses_covar())
    {
        throw std::invalid_argument("a factor augmented VAR has no psi block; its state "
                                    "innovation covariance is the Wishart precision alone");
    }

    // The observed half is data, so this is the one dimension check that is
    // about the sample rather than about the model. A file that leaves it out
    // has described a dynamic factor model and should say so in the algorithm.
    if (n_obs == 0)
    {
        throw std::invalid_argument("a factor augmented VAR must have at least one observed "
                                    "factor (n_obs_factors); a model with none is a dynamic "
                                    "factor model, and DfmNormalGamma estimates it");
    }
    require_shape(train.f_obs, tt, n_obs, "the observed factors (f_obs)");

    // The free loadings, in the row-major order the sampler draws them -- the
    // FAVAR count, never n_lambda(). The two agree at more than one dimension
    // (see VarSpec::n_favar_lambda()), so this length passing is not on its own
    // evidence that the file was written to a FAVAR's convention.
    validate_normal_block(lambda_prior, initial.lambda,
                          static_cast<arma::uword>(spec.n_favar_lambda()), "lambda");

    if (use_a())
    {
        validate_normal_block(a_prior, initial.a, static_cast<arma::uword>(spec.n_favar_a()), "a");
    }

    require_length(u_sigma_prior.shape, k, "gamma prior shape of the idiosyncratic precision");
    require_length(u_sigma_prior.rate, k, "gamma prior rate of the idiosyncratic precision");
    require_gamma_values(u_sigma_prior, "the idiosyncratic precision");
    require_length(initial.u_sigma_inv, k, "initial idiosyncratic precision");
    if (initial.u_sigma_inv.min() <= 0.0)
    {
        throw std::invalid_argument("every element of the initial idiosyncratic precision must be "
                                    "positive");
    }

    // The state innovation precision is a matrix here, not a diagonal, which is
    // what separates this model from every DFM. Both of these are therefore
    // n_state square rather than n_state long, and that is the shape a file
    // written against a DFM gets wrong.
    require_square(v_sigma_prior.scale, n_state,
                   "Wishart prior scale of the state innovation precision");
    require_symmetric(v_sigma_prior.scale, "Wishart prior scale of the state innovation precision");
    require_square(initial.v_sigma_inv, n_state, "initial state innovation precision");

    if (static_cast<arma::uword>(v_sigma_prior.df) < n_state)
    {
        throw std::invalid_argument(
            "the Wishart prior on the state innovation precision needs at least " +
            std::to_string(n_state) + " degrees of freedom, one per state element, got " +
            std::to_string(v_sigma_prior.df));
    }
}

void VecKlgs2010Input::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);

    require_identified_structural(spec, true, "a Wishart prior on the error precision");

    // Both schemes act on the columns of the SUR design matrix -- SSVS by
    // moving a column's prior between spike and slab, BVS by masking the column
    // itself -- and this sampler exists precisely because it never builds one.
    // Refused rather than silently ignored: a file that asks for selection and
    // gets none back is output that looks like output.
    if (spec.uses_varsel())
    {
        throw std::invalid_argument(
            "variable selection is not implemented for the non-SUR Koop, Leon-Gonzalez and "
            "Strachan (2010) sampler, which draws the coefficients without forming the design "
            "matrix a selection scheme would act on; VecNormalWishart is the same model in SUR "
            "form and carries both schemes");
    }

    if (use_a())
    {
        // One column per regressor, not k -- that is the whole difference from
        // the other VECs, so it is worth a check of its own rather than
        // require_stacked_regressors(). A model whose only regressor is the
        // error correction term has no `x` at all, and an empty one is right.
        const arma::uword n_x = static_cast<arma::uword>(spec.n_x_vec());
        if (n_x > 0)
        {
            require_shape(train.x, tt, n_x, "compact regressors x");
        }

        validate_normal_block(a_prior, initial.a, static_cast<arma::uword>(n_a()), "a");
    }

    if (use_beta())
    {
        validate_constant_coint_block(spec, train, beta_prior, initial.beta, tt, use_a());
        validate_coint_loadings_prior(spec, a_prior);
        refuse_coint_g(beta_prior);
    }

    validate_wishart_block(u_sigma_prior, initial.u_sigma_inv, k);
}

void VecNormalGammaInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n_a = train.nparams();

    require_identified_structural(spec, use_psi(), "a covariance block");

    if (use_a())
    {
        validate_vec_coefficients(spec, train, varsel_prior, initial.a_lambda, tt, k, n_a,
                                  use_beta());
        validate_normal_block(a_prior, initial.a, n_a, "a");
        if (spec.varsel == VarSelection::ssvs)
        {
            validate_ssvs_normal_prior(varsel_prior, a_prior, "a");
        }
    }

    if (use_beta())
    {
        validate_constant_coint_block(spec, train, beta_prior, initial.beta, tt, use_a());
        validate_coint_loadings_prior(spec, a_prior);
        refuse_coint_g(beta_prior);
    }

    if (use_psi())
    {
        const arma::uword n_psi = static_cast<arma::uword>(spec.n_psi());
        validate_normal_block(psi_prior, initial.psi, n_psi, "psi");

        if (spec.uses_varsel())
        {
            validate_varsel(psi_varsel_prior, initial.psi_lambda, n_psi, spec.varsel, "psi");
        }
        if (spec.varsel == VarSelection::ssvs)
        {
            validate_ssvs_normal_prior(psi_varsel_prior, psi_prior, "psi");
        }
    }

    require_length(u_sigma_prior.shape, k, "gamma prior shape of the error precision");
    require_length(u_sigma_prior.rate, k, "gamma prior rate of the error precision");
    require_gamma_values(u_sigma_prior, "the error precision");
    require_square(initial.u_sigma_inv, k, "initial error precision");
    require_diagonal(initial.u_sigma_inv, "initial error precision");
}

void VecNormalStochvolInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n_a = train.nparams();

    require_identified_structural(spec, use_psi(), "a covariance block");

    if (spec.varsel == VarSelection::ssvs)
    {
        throw std::invalid_argument("SSVS is not implemented for a stochastic volatility model; "
                                    "expected one of none, bvs");
    }

    if (use_a())
    {
        validate_vec_coefficients(spec, train, varsel_prior, initial.a_lambda, tt, k, n_a,
                                  use_beta());
        validate_normal_block(a_prior, initial.a, n_a, "a");
    }

    if (use_beta())
    {
        validate_constant_coint_block(spec, train, beta_prior, initial.beta, tt, use_a());
        validate_coint_loadings_prior(spec, a_prior);

        // The G the loadings' prior is scaled by, when the file fixes it. See
        // ConstantCointSpacePrior::g_inv.
        if (!beta_prior.g_inv.is_empty())
        {
            require_square(beta_prior.g_inv, k, "G^-1 of the cointegration space prior");
            require_symmetric(beta_prior.g_inv, "G^-1 of the cointegration space prior");
            require_positive_definite(beta_prior.g_inv, "G^-1 of the cointegration space prior");
        }
    }

    if (use_psi())
    {
        const arma::uword n_psi = static_cast<arma::uword>(spec.n_psi());
        validate_normal_block(psi_prior, initial.psi, n_psi, "psi");

        if (spec.uses_varsel())
        {
            validate_varsel(psi_varsel_prior, initial.psi_lambda, n_psi, spec.varsel, "psi");
        }
    }

    validate_stochvol_block(u_sigma_prior, initial.h_sigma, initial.h_init, initial.h, k, tt);

    // The random walk differences h against its own lag, so a single period
    // leaves nothing to difference.
    if (tt < 2)
    {
        throw std::invalid_argument("a stochastic volatility model needs at least two periods");
    }
}

void VecTvpWishartInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n_a = train.nparams();

    require_identified_structural(spec, true, "a Wishart prior on the error precision");

    validate_tvp_preconditions(spec.varsel, VarSelection::none, tt);

    if (use_a())
    {
        validate_vec_coefficients(spec, train, a_varsel_prior, initial.a_lambda, tt, k, n_a,
                                  use_beta());
        validate_tvp_block(a_prior, initial.a, initial.a_sigma_inv, initial.a_init, n_a, tt,
                           "coefficient", "a");
    }

    if (use_beta())
    {
        validate_tvp_coint_block(spec, train, beta_prior, initial.beta, initial.beta_init, tt,
                                 use_a());
    }

    // Unlike VecTvpGamma this model has no psi block: the covariance is carried
    // by the Wishart precision alone, so there is nothing here to match against
    // spec.n_psi().
    validate_wishart_block(u_sigma_prior, initial.u_sigma_inv, k);
}

void VecTvpGammaInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n_a = train.nparams();

    require_identified_structural(spec, use_psi(), "a covariance block");

    validate_tvp_preconditions(spec.varsel, psi_varsel, tt);

    if (use_a())
    {
        validate_vec_coefficients(spec, train, a_varsel_prior, initial.a_lambda, tt, k, n_a,
                                  use_beta());
        validate_tvp_block(a_prior, initial.a, initial.a_sigma_inv, initial.a_init, n_a, tt,
                           "coefficient", "a");
    }

    if (use_beta())
    {
        validate_tvp_coint_block(spec, train, beta_prior, initial.beta, initial.beta_init, tt,
                                 use_a());
    }

    if (use_psi())
    {
        const arma::uword n_psi = static_cast<arma::uword>(spec.n_psi());
        validate_tvp_block(psi_prior, initial.psi, initial.psi_sigma_inv, initial.psi_init, n_psi,
                           tt, "psi", "psi");

        if (uses_psi_varsel())
        {
            validate_varsel(psi_varsel_prior, initial.psi_lambda, n_psi, psi_varsel, "psi");
        }
    }

    require_length(u_sigma_prior.shape, k, "gamma prior shape of the error precision");
    require_length(u_sigma_prior.rate, k, "gamma prior rate of the error precision");
    require_gamma_values(u_sigma_prior, "the error precision");
    require_square(initial.u_omega_inv, k, "initial error precision");
    require_diagonal(initial.u_omega_inv, "initial error precision");
}

void VecTvpStochvolInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
    const arma::uword n_a = train.nparams();

    require_identified_structural(spec, use_psi(), "a covariance block");

    validate_tvp_preconditions(spec.varsel, psi_varsel, tt);

    if (use_a())
    {
        validate_vec_coefficients(spec, train, a_varsel_prior, initial.a_lambda, tt, k, n_a,
                                  use_beta());
        validate_tvp_block(a_prior, initial.a, initial.a_sigma_inv, initial.a_init, n_a, tt,
                           "coefficient", "a");
    }

    if (use_beta())
    {
        validate_tvp_coint_block(spec, train, beta_prior, initial.beta, initial.beta_init, tt,
                                 use_a());
    }

    if (use_psi())
    {
        const arma::uword n_psi = static_cast<arma::uword>(spec.n_psi());
        validate_tvp_block(psi_prior, initial.psi, initial.psi_sigma_inv, initial.psi_init, n_psi,
                           tt, "psi", "psi");

        if (uses_psi_varsel())
        {
            validate_varsel(psi_varsel_prior, initial.psi_lambda, n_psi, psi_varsel, "psi");
        }
    }

    validate_stochvol_block(u_sigma_prior, initial.h_sigma, initial.h_init, initial.h, k, tt);
}

} // namespace bayests
