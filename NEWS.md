# bvartools 0.3.0

This is a transition release. It is the functionality of the previous CRAN
version with the fixes listed below, and it is the last release before 1.0.0
reorganises the package around a different set of functions. Nothing here stops
working.

* **Functions that bvartools 1.0.0 does not have any more announce themselves.**
  The first time in a session that such a function is used it emits a message
  naming its successor, or saying that it has no replacement. The message is
  shown once per function per session and can be switched off with
  `options(bvartools.transition.messages = FALSE)`. It is a message rather than
  a warning, so it cannot become an error under `options(warn = 2)`.

    Renamed: `gen_var` to `create_bvarmodel`, `gen_vec` to `create_bvecmodel`,
  `bvec_to_bvar` to `vec_to_var`, `kalman_dk` to `kalman_durbin_koopman_2002`,
  `stochvol_ksc1998` to `stochvol_ksc_1998`, `stochvol_ocsn2007` to
  `stochvol_ocsn_2007`, `stoch_vol` moved into `stochvol_ksc_1998`, and `bvs` to
  `post_bvs`. Replaced by a different workflow: `draw_posterior`, `bvarpost`
  and `bvecpost`, which become `add_posterior_coefficients` alongside
  `add_posterior_forecasts` and `add_posterior_loglik`. Removed with no
  successor: `post_normal_covar_const`, `post_normal_covar_tvp`, and the whole
  dynamic factor model branch -- `dfm`, `dfmpost`, `gen_dfm`, their methods, and
  the data set `bem_dfmdata`, which cannot announce itself.

    Methods on the renamed classes `bvar`, `bvec` and `bvarlist` are silent,
  because the generic that dispatches them is unchanged and only the class is
  renamed, to `bvarmodel`, `bvecmodel` and `modellist`. So are the functions
  that keep their name in 1.0.0 but take the reorganised model object:
  `add_priors`, `bvar`, `bvec`, `irf`, `fevd`, `inclusion_prior`,
  `minnesota_prior` and `ssvs_prior`. The new vignette, *Moving from bvartools
  0.3.0 to 1.0.0*, lists all of it.

* **Fixed: `stochvol_ksc1998` and `stochvol_ocsn2007` failed on an observation
  far out in the tails of every mixture component.** Both sampled the mixture
  indicator from weights formed as densities and normalised by their sum. Where
  the log of the squared observation lies far enough below the log-volatility --
  about 112 for the seven components of Kim, Shephard and Chib, about 158 for
  the ten of Omori, Chib, Shephard and Nakajima -- every density underflows to
  zero, the row sums to zero, the weights become `NaN`, and the sampled
  indicator runs one past the last component, ending the call with
  `Mat::elem(): index out of bounds`. The weights are now formed in logs and
  shifted by their row maximum before they are exponentiated, and the indicator
  is clamped to the components that exist.

    This is algebraically the same calculation, and draws are unchanged where
  they were being produced at all: verified bit for bit against the previous
  implementation over the `us_macrodata` series from a fixed seed. `stoch_vol`
  is a wrapper for the first of the two and inherits the fix.

* **Fixed: neither function checked the size of `sigma`, `h_init` or
  `constant`.** They were indexed on trust, so a vector of the wrong length was
  reported as `Mat::elem(): index out of bounds` instead of as a statement about
  the argument. Each is now checked against the number of columns of `y`.

* **Fixed: `irf` and `fevd` produced reduced form quantities from a structural
  model.** A structural model keeps its contemporaneous block separately, so its
  coefficient draws are the structural `A_i` and its covariance draws the
  covariance of the structural errors. The forecast error, orthogonalised and
  generalised recursions want the reduced form. Given the structural quantities
  they returned numbers that belong to no model at all, and the same numbers for
  all three types, since the structural error covariance makes the
  orthogonalisation degenerate. `irf` with `type` of `"feir"`, `"oir"` or
  `"gir"`, and `fevd` with `"oir"` or `"gir"`, now stop on a structural model
  and say why. The structural types `"sir"` and `"sgir"` are unchanged, as is
  every reduced form model.

* **Fixed: the structural variance decomposition ignored the variances of the
  structural shocks.** `fevd` with `type = "sir"` used `A_0^-1` as the impulse
  matrix and `A_0^-1 A_0^-1'` as the forecast error covariance, leaving the
  covariance of the structural errors out of both, so every structural shock was
  decomposed as if it had unit variance. Weight moved to whichever shock loads
  most heavily in `A_0`, and since the shares are normalised they still summed
  to one, so the output gave nothing away. On a three variable example whose
  shock variances stand at 4, 1 and 0.25 the reported shares were 0.21, 0.13 and
  0.65 where they should be 0.74, 0.12 and 0.14. The impulse matrix is now
  `A_0^-1 chol(Sigma)'` and the forecast error covariance `A_0^-1 Sigma A_0^-1'`.
  `"sgir"` already carried `Sigma` and is unchanged, as is every reduced form
  type.

* **Fixed: `gen_vec` stopped on seasonal terms for data of frequency one.** It
  warned that no seasonal dummies are generated and then added them anyway,
  failing with `object 'seas' not found` because the dummies had never been
  built. It now does what the warning says. `gen_var` was never affected.

* Added a test suite, which the package did not have before. It covers the
  announcements and carries one regression test per fix above.


* Carried over from the development version that was never released as 0.2.5.
  All six survive into 1.0.0 and so announce nothing.

    * Added function `covar_vector_to_matrix`.
    * Added function `sur_const_to_tvp`.
    * Updated `Rcpp` dependency in DESCRIPTION file to version 1.0.12.
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
