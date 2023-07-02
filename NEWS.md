# 0.2.3

* Implement recursive iterations for forecasts in C++.
* Replace erroneous `|` in C++ sampling functions by `||`.

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