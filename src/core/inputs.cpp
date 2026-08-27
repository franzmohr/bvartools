// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "bayests/inputs.h"

#include <stdexcept>
#include <string>

namespace
{

std::string dims(const arma::mat &m)
{
    return std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols);
}

void require_square(const arma::mat &m, arma::uword side, const char *what)
{
    if (m.n_rows != side || m.n_cols != side)
    {
        throw std::invalid_argument(std::string(what) + " must be " + std::to_string(side) +
                                    "x" + std::to_string(side) + ", got " + dims(m));
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

/// The checks a selection block needs whichever coefficient vector it applies
/// to. `n` is the length of that vector, and the labels name it so the message
/// says whether it was the coefficients or the covariance block that was wrong.
void validate_varsel(const VarSelPrior &prior, const arma::vec &initial_lambda,
                     arma::uword n, VarSelection scheme, const char *block)
{
    const std::string what(block);

    require_length(prior.inprior, n, (what + " prior inclusion probabilities").c_str());
    require_length(initial_lambda, n, (what + " initial inclusion indicators").c_str());

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
void validate_tvp_block(const RandomWalkPrior &prior, const arma::mat &path,
                        const arma::mat &sigma_inv, const arma::vec &init, arma::uword n,
                        arma::uword tt, const char *noun, const char *name)
{
    const std::string thing(noun);
    const std::string vec(name);

    require_shape(path, n, tt, ("initial " + thing + " path").c_str());
    require_square(sigma_inv, n, ("initial precision of the " + thing + " innovations").c_str());
    require_length(init, n, ("initial value of " + vec + " before the sample").c_str());

    require_length(prior.sigma.shape, n, ("prior shape of the " + thing + " innovations").c_str());
    require_length(prior.sigma.rate, n, ("prior rate of the " + thing + " innovations").c_str());
    require_length(prior.initial_state.mu, n, ("prior mean of " + vec + " before the sample").c_str());
    require_square(prior.initial_state.v_inv, n,
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

/// rho scales the cointegration state path itself, so a value outside (0, 1]
/// either reverses the sign of the relation from period to period or lets it
/// grow without bound. One is the random walk bvartools' .bvectvpalg uses.
void validate_tvp_coint_rho(double rho)
{
    if (!(rho > 0.0 && rho <= 1.0))
    {
        throw std::invalid_argument(
            "the autoregression of the cointegration state equation (rho) must lie in (0, 1], "
            "got " + std::to_string(rho));
    }
}

} // namespace
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
    }

    if (u_sigma_prior.df <= 0)
    {
        throw std::invalid_argument("Wishart prior degrees of freedom must be positive");
    }
    require_square(u_sigma_prior.scale, k, "Wishart prior scale");
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

    require_length(u_sigma_prior.shape, k, "gamma prior shape of the error precision");
    require_length(u_sigma_prior.rate, k, "gamma prior rate of the error precision");
    require_square(initial.u_sigma_inv, k, "initial error precision");
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

    require_shape(initial.h, tt, k, "initial log-volatility");
    require_length(initial.h_init, k, "initial value of the log-volatility before the sample");
    require_length(initial.h_sigma, k, "initial variance of the log-volatility innovations");

    // The random walk differences h against its own lag, so a single period
    // leaves nothing to difference.
    if (tt < 2)
    {
        throw std::invalid_argument("a stochastic volatility model needs at least two periods");
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
    require_square(initial.u_omega_inv, k, "initial error precision");
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
    require_length(u_sigma_prior.state.sigma.shape, k, "prior shape of the log-volatility innovations");
    require_length(u_sigma_prior.state.sigma.rate, k, "prior rate of the log-volatility innovations");
    require_length(u_sigma_prior.state.initial_state.mu, k, "prior mean of the log-volatility before the sample");
    require_square(u_sigma_prior.state.initial_state.v_inv, k, "prior precision of the log-volatility before the sample");

    require_length(initial.h_sigma, k, "initial variance of the log-volatility innovations");
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

    validate_tvp_coint_rho(prior.rho);
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
    require_length(prior.state.sigma.shape, k, "prior shape of the log-volatility innovations");
    require_length(prior.state.sigma.rate, k, "prior rate of the log-volatility innovations");
    require_length(prior.state.initial_state.mu, k,
                   "prior mean of the log-volatility before the sample");
    require_square(prior.state.initial_state.v_inv, k,
                   "prior precision of the log-volatility before the sample");

    require_length(h_sigma, k, "initial variance of the log-volatility innovations");
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
    }

    if (use_beta())
    {
        validate_constant_coint_block(spec, train, beta_prior, initial.beta, tt, use_a());
    }

    validate_wishart_block(u_sigma_prior, initial.u_sigma_inv, k);
}

void DfmNormalGammaInput::validate() const
{
    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword tt = checked_periods(spec, train);
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

    // The free loadings, in the row-major order the sampler draws them.
    const arma::uword n_lambda = static_cast<arma::uword>(spec.n_lambda());
    validate_normal_block(lambda_prior, initial.lambda, n_lambda, "lambda");

    if (use_a())
    {
        validate_normal_block(a_prior, initial.a, static_cast<arma::uword>(spec.n_factor_a()), "a");
    }

    // Both precisions are diagonal, so both arrive as the diagonal rather than
    // as a matrix -- and a starting value that is not positive is not one, since
    // the first factor draw inverts it.
    require_length(u_sigma_prior.shape, k, "gamma prior shape of the idiosyncratic precision");
    require_length(u_sigma_prior.rate, k, "gamma prior rate of the idiosyncratic precision");
    require_length(initial.u_sigma_inv, k, "initial idiosyncratic precision");

    require_length(v_sigma_prior.shape, n, "gamma prior shape of the factor innovation precision");
    require_length(v_sigma_prior.rate, n, "gamma prior rate of the factor innovation precision");
    require_length(initial.v_sigma_inv, n, "initial factor innovation precision");

    if (initial.u_sigma_inv.min() <= 0.0)
    {
        throw std::invalid_argument("every element of the initial idiosyncratic precision must be "
                                    "positive");
    }
    if (initial.v_sigma_inv.min() <= 0.0)
    {
        throw std::invalid_argument("every element of the initial factor innovation precision must "
                                    "be positive");
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
    }

    if (use_beta())
    {
        validate_constant_coint_block(spec, train, beta_prior, initial.beta, tt, use_a());
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

    require_length(u_sigma_prior.shape, k, "gamma prior shape of the error precision");
    require_length(u_sigma_prior.rate, k, "gamma prior rate of the error precision");
    require_square(initial.u_sigma_inv, k, "initial error precision");
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
    require_square(initial.u_omega_inv, k, "initial error precision");
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
