
# bvartools 0.3.0

* Migrated `bvarpost` into `add_posterior_coefficients.bvarmodel`
* Added the core layer of BayesTS -- the samplers' declarations in `inst/include/bayests` and their numerics in `src/core` -- as vendored source. See `src/core/VENDORED.md`.
* All six `bvarmodel` algorithms -- `VarNormalWishart`, `VarNormalGamma`, `VarNormalStochvol`, `VarTvpGamma`, `VarTvpWishart` and `VarTvpStochvol` -- are now the samplers from the vendored core layer, and each model's three entry points are consolidated into one source file, `src/<Model>.cpp`. No posterior simulation is implemented in this package any more. Draws agree with the previous implementations to a rounding error (worst relative difference 5.6e-14 across eighteen model configurations; all inclusion indicators identical).
* BVS selection of the contemporaneous error coefficients was not selecting in `VarTvpGamma` and `VarTvpStochvol`. Those coefficients are a path, one column per period, but the candidate the likelihood ratio scored zeroed the coefficient in the first period alone rather than across the sample -- a single linear index into the matrix instead of a whole row. The two candidates being compared therefore differed in one period out of `T`, the ratio between them was correspondingly close to zero, and inclusion won essentially every time: at `T = 24` all three free coefficients stayed in for all of a recorded 80 draws, against 26 of 240 once the candidate spans the sample. The mask applied to the draws that are reported was never affected, so the coefficients themselves were consistent with the indicators -- it was the indicators that carried no information. Fixed in BayesTS. Draws of `VarTvpGamma` and `VarTvpStochvol` with both a covariance block and variable selection change; no other model or configuration is affected.
* The two variable selection schemes of the vendored core now live in `src/core/algorithms/ssvs.cpp` and `src/core/algorithms/bvs.cpp`, one file per scheme, in place of thirteen near-identical copies spread across the six samplers. Draws from a given seed are unchanged for every model and configuration.
* `VarTvpWishart` draws change beyond a rounding error. Its sampler built the block diagonal error precision once, before the first iteration, and then kept using it although the precision is redrawn every iteration, so the coefficient smoother conditioned on the initial error covariance for the whole chain. Fixed in BayesTS.
* The `VarNormalWishart` algorithm for objects of class `bvarmodel` is now the sampler from the vendored core layer. `bvarpost`, `add_posterior_forecasts` and `add_posterior_loglik` are unchanged and still dispatch to `.VarNormalWishartCoefficients`, `.VarNormalWishartForecasts` and `.VarNormalWishartLogLik`, which have become bindings over `bayests::VarNormalWishartSampler`; the R model object and the posterior it returns keep their shape. Draws agree with the previous implementation to a rounding error (relative differences below 1e-13, inclusion indicators identical), which comes from the core reusing one Cholesky factor for the posterior mean and the draw and inverting the Wishart scale with `inv_sympd`. Inconsistent input is now reported by name -- "initial error precision must be 3x3, got 2x2" -- rather than as a matrix dimension mismatch.
* The simulation smoother of Durbin and Koopman (2002) now comes from that core layer, and the TVP samplers call it directly instead of through the package's own C++ interface. The interface wrapped every call in an `Rcpp::RNGScope`, which rewound the random number stream and handed the same numbers out twice per iteration, so posterior draws of TVP models change from the second iteration on. The C++ entry point `bvartools::kalman_durbin_koopman_2002` no longer exists; it was never available from R.
* Fixed the same random number stream rewind in `post_bvs`, which called `bvartools::sur_const_to_tvp` for time varying parameters, and in the VEC sampler, which called `bvartools::stochvol_ocsn2007_internal` for stochastic volatility. Both now call the function directly. Draws of the affected models change. The `bvartools::sur_const_to_tvp` wrapper itself is unchanged and remains available to other packages; `bvartools::stochvol_ocsn2007_internal` has since been removed, see the next entry.
* The stochastic volatility draw of Omori, Chib, Shephard and Nakajima (2007) now also comes from the vendored core layer, as `stochvol_ocsn_2007()`, and replaces `src/stochvol_ocsn2007_internal.cpp`. It validates its arguments instead of indexing them on trust, and draws the mixture indicator in logs, so an observation far out in the tails of all ten mixture components no longer normalises to `NaN` and selects a component that does not exist. The C++ entry point `bvartools::stochvol_ocsn2007_internal` no longer exists; it was never available from R. Draws of the VEC sampler with stochastic volatility change by a rounding error, because the posterior mean is now formed from the same Cholesky factor as the draw rather than by a separate LU solve.
* The exported R function `stochvol_ocsn2007` is a separate implementation and received the same three fixes. It validates `y`, `h`, `sigma`, `h_init` and `constant` rather than reporting a bad argument as a failed matrix factorisation or, for `h_init`, not at all; it draws the mixture indicator in logs, which stops an observation far out in the tails of all ten components from returning a matrix of `NA`; and it reuses one Cholesky factor for the posterior mean and the draw. Draws from a given seed are unchanged up to a rounding error. The documentation of the mixture, previously described as having seven components, now says ten.
* Fixed missing transformations for structural modesl in `irf.bvar` and `fevd.bvar`.
* Added function `choose_best_model`.
* Added function `create_first_difference_matrix`.
* Added function `create_second_difference_matrix`.
* Added functions `create_bvarmodel`, `create_bvecmodel` and `create_dfmodel` to replace `gen_var`, `gen_vec` and `gen_dfm`, respectively, in the future.
* Added `generate_lower_block_diagonal` for faster simulation of autocorrelated TVP coefficients.
* Updated functions for dynamic factor models to the revised design pattern.
* Added a warning message that the functionality for dynamic factor models will be exported to a separate package with package in the future.
* Add a generic function `mean_absolute_forecast_error` for MAFE calculation.
* Add a generic function `prediction_matrix` which helps with forecast generation.
* Add a generic function `selection_criteria` to calculate information criteria for model selection.
* Add a generic function `forecast_errors` to generate forecast errors.
* Add a generic function `add_initial_values` to separate initial value generation from prior specification.
* Added `gen_artificial_vec` to generate artificial data sets for algorithm and model testing.
* Added `gen_artificial_var` to generate artificial data sets for algorithm and model testing.
* Added `coint_kls2010_reparameterise_two` for more convenient data transformation.
* Added `coint_prepare_sur_data` for more convenient input data preparation for cointegration simulation.
* `plot.bvar` and `plot.bvec` allow to specify whether a horizontal line should be added or not.
* Added convenience function `covar_vector_to_matrix`.
* Added convenience function `sur_const_to_tvp`.
* General updates in documentations including update of `Rcpp` dependency in DESCRIPTION file to version 1.0.12.
* Added `post_gamma_state_variance` for posterior simulation of constant error variances of the state equation.
* Added `post_gamma_measurement_variance` for posterior simulation of constant error variances of the measurement equation.
* Renamed `.prep_covar_data` to `covar_prepare_data` and made it visible in R and also callable from C++.

# bvartools 0.2.4

* Using an updated version of `Rcpp` to address an issue with `Rcpp::stop`.
* `stochvol_ocsn2007` can handle multi-column input.
* `stochvol_ksc1998` can handle multi-column input.
* Added `post_normal_covar_tvp` for posterior simulation of time varying, lower triangular covariance matrices.
* Added `post_normal_covar_const` for posterior simulation of constant, lower triangular covariance matrices.

# bvartools 0.2.3

* Fixed alias issue resulting from use of `roxygen2`.
* Made `kalman_dk` callable from C++.
* Stochastic volatility algorithms allow to set the offsetting constant manually.
* Changed `stoch_vol` to a wrapper for `stochvol_ksc1998`.
* Added stochastic volatility algorithm of Kim et al. (1998) in a separate function `stochvol_ksc1998`.
* Added stochastic volatility algorithm of Omori et al. (2007) in function `stochvol_ocsn2007`.
* Fixed bug with detection of deterministic terms in `bvar`.
* Implemented recursive iterations for forecasts in C++.
* Replaced erroneous `|` in C++ sampling functions by `||`.

# bvartools 0.2.2

* Addressed CRAN NOTE on CITATION file
* Addressed the CRAN NOTE "Specified C++11: please drop specification unless essential" by dropping the specification from "src/Makevars"
* Improved the treatment of `bvar` and `bvec` objects if Gibbs sampler fails.
* Fix erroneous SUR-matrix generation for VEC models with r = 0 in `.bvecalg`.
* Fix bug in `.bvecalg` and `.bvectvpalg` with the storing of posterior draws of beta.
* Fix bug of `predict.bvar`, which could not handle only VARX models with contemporaneous exogenous variables only.
* Model plot functions support boxplots.
* Fix typos in documentation.

# bvartools 0.2.1

* Added functionality for the simulation of models with time varying parameters, both for VAR and VEC models.
* Added functionality for the simulation of models with stochastic volatility, both for VAR and VEC models.
* Added a plot function for classes `bvar` and `bvec` for visual inspection of posterior draws.
* Changed the generation of the output object in the Gibbs sampler functions `bvaralg` and `bvecalg` to make them more stable for especially large output.
* Changed `draw_posterior` to a generic function and added the corresponding methods for BVAR, BVEC and DFM input.
* Changed `irf` and `fevd` to generic functions.
* Corrected typos in documentation.
* `thin_posterior` methods were renamed to `thin` and are now methods of `coda::thin`.
* Function `irf` allows to specify the size of a shock.
* Fixed a bug in `ssvs_prior` concerning BVEC models.
* Fixed a bug with the prior in the BVEC algorithm.

# bvartools 0.2.0

* Changed `thin_posterior` to a generic function and added methods for BVAR, BVEC and dynamic factor model input.
* Changed `add_prior` to a generic function and added methods for BVAR, BVEC and dynamic factor model input.
* Added funcionality to estimate dynamic factor models (DFM).
* `predict` requires to specify an object of class `ts` as input for argument `exogen`.
* Additioal argument checks for `add_priors` methods.
* Updated documentation in `minnesota_prior` and for `add_prior` methods.
* Using \doi instead of \url in documentation

# bvartools 0.1.0

* Omitted package `Matrix` from "Imports"" in DESCRIPTION, which caused a note in version 0.0.3.
* Added function `bvarpost` for posterior simulation of BVAR models.
* Added function `bvecpost` for posterior simulation of BVEC models.
* Added function `draw_posterior` for estimation of multiple models.
* Fixed erroneous calculation of structural forecast error variance decompositions.
* More specification checks and increased robustness against erroneous model specificaions.
* Function `fevd` calculates FEVDs based on means of posterior draws of FEVDs and not based on the means of the coefficient draws.
* Function `bvar` and `summary.bvar` can deal with inclusion parameters.
* Added funtion `add_priors` for easier construction of prior matrices for multiple models.
* `gen_var` and `gen_vec` can produce multiple models.
* Changed all argument names of `predict.bvar` to lower cases.

# bvartools 0.0.3

* Changed all argument names of `post_normal`, `post_normal_sur`, `post_coint_kls` and `post_coint_kls_sur` to lower case letters.
* Replaced output element in function `ssvs` from `V_i` to `v_i`.
* Refined function `minnesota_prior` and added additional functionaliy.
* Fixed error message when creating seasonal dummies with `gen_var` and `gen_vec`.
* New data set `us_macrodata`.
* Added additional checks in `gen_vec`.
* Added functions `inclusion_prior` for the calculation of inclusion probability priors as used by `bvs` and `ssvs`.
* Added `summary` functions.
* Fixed conversion and collection of exogenous regressors in `bvec_to_bvar`.
* Fixed detection of deterministic terms in `bvec_to_bvar`.
* Updated documentation in `kalman_dk`.
* `irf` contains a new argument `keep_draws`.
* Additional checks in `post_normal`, `post_normal_sur`, `post_coint_kls` and `post_coint_kls_sur`.
* Adapt vignette `bvec`.
* Added `loglik_normal` for the calculation of a multivariate normal log-likelihood.

# bvartools 0.0.2

* Updated vignette `ssvs` after the introduction of function `ssvs_prior`.
* Added `ssvs_prior` for the calculation of prior matrices for the SSVS algorithm.
* Added `minnesota_prior` for the calculation of the Minnesota prior.
* Use unsigned integers for indices in Cpp code to address warnings during installation.
* Better error handling in `irf`.
* In `post_coint_kls_sur` the prior matrix `g_i` can be time varying.
* `bvar` and `predict` also work only with deterministic terms, i.e. p can be zero.
* Use SVD to obtain a draw of beta in `post_coint_kls` and `post_coint_kls_sur`.
* `predict` allows for p = 1.
* Add legend to `plot.bvarfevd`.

# bvartools 0.0.1

* Initial release
