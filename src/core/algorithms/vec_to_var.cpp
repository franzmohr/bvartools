// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "core/algorithms/vec_to_var.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace bayests::core
{

namespace
{

void require_rows(const arma::mat &values, const arma::uword rows, const char *what)
{
    if (values.n_rows != rows)
    {
        throw std::invalid_argument(std::string("posterior draws of ") + what + " must have " +
                                    std::to_string(rows) + " rows, got " +
                                    std::to_string(values.n_rows));
    }
}

} // namespace

VarSpec vec_to_var_spec(const VarSpec &spec)
{
    VarSpec out = spec;

    out.p = std::max(spec.p, 1);
    out.s = spec.m > 0 ? std::max(spec.s, 1) : 0;
    out.n = spec.n + (spec.uses_coint() ? spec.n_restricted : 0);

    out.rank = 0;
    out.n_restricted = 0;
    out.varsel = VarSelection::none;

    return out;
}

VarNormalWishartDraws vec_to_var_coefficients(const VarSpec &spec,
                                              const VecNormalWishartDraws &draws)
{
    // Deliberately not spec.validate(): that also insists on a positive
    // iteration count, which is chain bookkeeping a forecast-only input does not
    // carry and this transformation never reads -- the number of draws comes
    // from the posterior itself. What does matter is that the dimensions can be
    // cast to uword and indexed with.
    if (spec.k <= 0)
    {
        throw std::invalid_argument("model must have at least one endogenous variable (k)");
    }
    if (spec.p < 0 || spec.m < 0 || spec.s < 0 || spec.n < 0 || spec.rank < 0 ||
        spec.n_restricted < 0)
    {
        throw std::invalid_argument("model dimensions (p, m, s, n, rank, n_restricted) cannot be "
                                    "negative");
    }

    const VarSpec var_spec = vec_to_var_spec(spec);

    const arma::uword k = static_cast<arma::uword>(spec.k);
    const arma::uword m = static_cast<arma::uword>(spec.m);
    const arma::uword rank = static_cast<arma::uword>(spec.rank);
    const arma::uword k_ect = static_cast<arma::uword>(spec.k_ect());
    const arma::uword n_alpha = static_cast<arma::uword>(spec.n_alpha());
    const arma::uword n_structural = static_cast<arma::uword>(spec.n_structural());

    // Blocks of the VEC's coefficient matrix, and of the level VAR's. `p` and
    // `s` count level lags, so the VEC has one block fewer of each than they
    // say; the level counts come from vec_to_var_spec() rather than from spec
    // directly, so that its floor of one endogenous lag applies here too.
    const arma::uword n_gamma = static_cast<arma::uword>(spec.n_gamma());
    const arma::uword n_upsilon = static_cast<arma::uword>(spec.n_upsilon());
    const arma::uword n_a = static_cast<arma::uword>(var_spec.p);
    const arma::uword n_b = m > 0 ? static_cast<arma::uword>(var_spec.s) + 1 : 0;

    // Columns of each coefficient matrix, VEC and VAR.
    const arma::uword ncols_vec = k * n_gamma + m * n_upsilon + static_cast<arma::uword>(spec.n);
    const arma::uword ncols_var = k * n_a + m * n_b + static_cast<arma::uword>(var_spec.n);

    if (draws.u_sigma_inv.n_elem == 0)
    {
        throw std::invalid_argument("posterior draws of u_sigma_inv are missing");
    }
    require_rows(draws.u_sigma_inv, k * k, "u_sigma_inv");

    const arma::uword iterations = draws.iterations();
    const arma::uword nparams_vec = static_cast<arma::uword>(spec.nparams_per_period_vec());

    if (nparams_vec > 0)
    {
        if (!draws.has_a())
        {
            throw std::invalid_argument("the VEC has " + std::to_string(nparams_vec) +
                                        " coefficients but posterior draws of a are missing");
        }
        require_rows(draws.a, nparams_vec, "a");
        if (draws.a.n_cols != iterations)
        {
            throw std::invalid_argument("posterior draws of a and u_sigma_inv must come from the "
                                        "same chain, got " + std::to_string(draws.a.n_cols) +
                                        " and " + std::to_string(iterations) + " draws");
        }
    }

    const bool use_coint = spec.uses_coint();
    if (use_coint)
    {
        if (!draws.has_beta())
        {
            throw std::invalid_argument("the VEC has a cointegration relation of rank " +
                                        std::to_string(spec.rank) +
                                        " but posterior draws of beta are missing");
        }
        require_rows(draws.beta, static_cast<arma::uword>(spec.n_beta()), "beta");
        if (draws.beta.n_cols != iterations)
        {
            throw std::invalid_argument("posterior draws of beta and u_sigma_inv must come from "
                                        "the same chain, got " + std::to_string(draws.beta.n_cols) +
                                        " and " + std::to_string(iterations) + " draws");
        }
    }

    VarNormalWishartDraws out;
    out.u_sigma_inv = draws.u_sigma_inv;
    out.a = arma::mat(static_cast<arma::uword>(var_spec.nparams_per_period()), iterations);

    // Where the deterministic terms sit in each coefficient matrix.
    const arma::uword det_vec = k * n_gamma + m * n_upsilon;
    const arma::uword det_var = k * n_a + m * n_b;
    const arma::uword n_unrestricted = static_cast<arma::uword>(spec.n);
    const arma::uword n_restricted = use_coint ? static_cast<arma::uword>(spec.n_restricted) : 0;

    const arma::mat diag_k = arma::eye<arma::mat>(k, k);
    arma::mat vec_coef, var_coef(k, ncols_var), pi;

    for (arma::uword draw = 0; draw < iterations; draw++)
    {
        var_coef.zeros();

        // Pi = alpha beta', whose columns are the loadings on the levels of the
        // endogenous and unmodelled variables and on the restricted
        // deterministic terms, in that order.
        if (use_coint)
        {
            pi = arma::reshape(draws.a.submat(0, draw, n_alpha - 1, draw), k, rank) *
                 arma::trans(arma::reshape(draws.beta.col(draw), k_ect, rank));
        }

        // [Gamma_1 .. Gamma_p, Upsilon_0 .. Upsilon_s, C], recovered from the
        // column-major vec the sampler stored.
        if (ncols_vec > 0)
        {
            vec_coef = arma::reshape(draws.a.submat(n_alpha, draw, n_alpha + k * ncols_vec - 1, draw),
                                     k, ncols_vec);
        }

        // A_i = Gamma_i - Gamma_{i-1}, with the identity and Pi_y joining the
        // first block and Gamma_{p-1} appearing in A_p with a negative sign
        // only. Blocks the model does not have simply drop out of the sum, which
        // is what makes p <= 1 -- a VEC with no lagged differences, whose level
        // VAR is the single block A_1 = I + Pi_y -- fall out of the same loop.
        for (arma::uword i = 0; i < n_a; i++)
        {
            const arma::uword first = i * k;
            const arma::uword last = first + k - 1;

            if (i < n_gamma)
            {
                var_coef.cols(first, last) += vec_coef.cols(i * k, i * k + k - 1);
            }
            if (i > 0)
            {
                var_coef.cols(first, last) -= vec_coef.cols((i - 1) * k, i * k - 1);
            }
            else
            {
                var_coef.cols(first, last) += diag_k;
                if (use_coint)
                {
                    var_coef.cols(first, last) += pi.cols(0, k - 1);
                }
            }
        }

        // B_j = Upsilon_j - Upsilon_{j-1}, the same telescoping sum one block
        // later: the VEC regresses on the current difference dx_t, so the levels
        // it implies run from x_t to x_{t-s} and Pi_x lands on x_{t-1}.
        for (arma::uword j = 0; j < n_b; j++)
        {
            const arma::uword first = k * n_a + j * m;
            const arma::uword last = first + m - 1;
            const arma::uword ups = k * n_gamma + j * m;

            if (j < n_upsilon)
            {
                var_coef.cols(first, last) += vec_coef.cols(ups, ups + m - 1);
            }
            // Guarded on both sides, unlike the endogenous loop: with s = 0 the
            // VEC has no Upsilon at all, yet Pi_x still has to reach x_{t-1}, so
            // vec_to_var_spec() floors the level order at one and this block runs
            // for a j whose Upsilon_{j-1} does not exist.
            if (j > 0 && j - 1 < n_upsilon)
            {
                var_coef.cols(first, last) -= vec_coef.cols(ups - m, ups - 1);
            }
            if (j == 1 && use_coint)
            {
                var_coef.cols(first, last) += pi.cols(k, k + m - 1);
            }
        }

        // Deterministic terms are not differenced, so the unrestricted ones
        // carry over unchanged; the restricted ones arrive as the last columns
        // of Pi and become ordinary regressors of the level VAR.
        if (n_unrestricted > 0)
        {
            var_coef.cols(det_var, det_var + n_unrestricted - 1) =
                vec_coef.cols(det_vec, det_vec + n_unrestricted - 1);
        }
        if (n_restricted > 0)
        {
            var_coef.cols(det_var + n_unrestricted, ncols_var - 1) = pi.cols(k + m, k_ect - 1);
        }

        out.a.submat(0, draw, k * ncols_var - 1, draw) = arma::vectorise(var_coef);

        // Contemporaneous coefficients are relations among the errors, which the
        // transformation leaves alone. They stay last, as the VAR expects.
        if (n_structural > 0)
        {
            out.a.submat(k * ncols_var, draw, k * ncols_var + n_structural - 1, draw) =
                draws.a.submat(nparams_vec - n_structural, draw, nparams_vec - 1, draw);
        }
    }

    return out;
}

} // namespace bayests::core
