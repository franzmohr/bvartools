# Package index

## Model set-up

- [`bvartools_model`](https://franzmohr.github.io/bvartools/reference/bvartools_model.md)
  : The Structure of a Model Object
- [`create_bvarmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvarmodel.md)
  : Create a Vector Autoregressive Model
- [`create_bvecmodel()`](https://franzmohr.github.io/bvartools/reference/create_bvecmodel.md)
  : Create Vector Error Correction Models
- [`is_discount_model()`](https://franzmohr.github.io/bvartools/reference/discount_models.md)
  [`check_discount_specification()`](https://franzmohr.github.io/bvartools/reference/discount_models.md)
  : The Discounted Models
- [`use_expanding_window()`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.md)
  : Expanding Window Estimation
- [`use_expanding_window(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvarmodel.md)
  : Expanding Window Estimation
- [`use_expanding_window(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.bvecmodel.md)
  : Expanding Window Estimation
- [`use_expanding_window(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/use_expanding_window.modellist.md)
  : Expanding Window Estimation
- [`combine_models()`](https://franzmohr.github.io/bvartools/reference/combine_models.md)
  : Combine Models
- [`align_model_obs()`](https://franzmohr.github.io/bvartools/reference/align_model_obs.md)
  : Align Observations Across Models
- [`align_model_obs(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/align_model_obs.modellist.md)
  : Align Observations Across Models
- [`transform_variables()`](https://franzmohr.github.io/bvartools/reference/transform_variables.md)
  : Apply Transformations
- [`window(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/window.bvarmodel.md)
  : Time Series Windows
- [`window(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/window.bvecmodel.md)
  : Time Series Windows
- [`window(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/window.modellist.md)
  : Time Series Windows
- [`generate_artificial_var()`](https://franzmohr.github.io/bvartools/reference/generate_artificial_var.md)
  : Generate Artificial VAR Data
- [`generate_artificial_vec()`](https://franzmohr.github.io/bvartools/reference/generate_artificial_vec.md)
  : Generate Artificial VEC Data

## Priors and initial values

- [`add_priors()`](https://franzmohr.github.io/bvartools/reference/add_priors.md)
  : Add Priors to Bayesian Models
- [`add_priors(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_priors.bvarmodel.md)
  : Add Priors to Bayesian Models
- [`add_priors(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_priors.bvecmodel.md)
  : Add Priors to Bayesian Models
- [`add_priors(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/add_priors.expandingwindow.md)
  : Add Priors to Bayesian Models
- [`add_priors(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/add_priors.modellist.md)
  : Add Priors to Bayesian Models
- [`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_priors(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_initial_values(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_seed(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_coefficients(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_loglik(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_forecast_input(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_forecasts(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`minnesota_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`ssvs_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`inclusion_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`thin(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`print(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  : Objects for Externally Produced Forecasts
- [`add_priors(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_initial_values(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_seed(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_posterior_coefficients(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_posterior_loglik(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`thin(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`selection_criteria(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  : Estimation Steps on a Folder of Models
- [`minnesota_prior()`](https://franzmohr.github.io/bvartools/reference/minnesota_prior.md)
  : Minnesota Prior
- [`minnesota_prior(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/minnesota_prior.bvarmodel.md)
  : Minnesota Prior
- [`minnesota_prior(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/minnesota_prior.bvecmodel.md)
  : Minnesota Prior
- [`minnesota_prior(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/minnesota_prior.modellist.md)
  : Minnesota Prior
- [`ssvs_prior()`](https://franzmohr.github.io/bvartools/reference/ssvs_prior.md)
  : Stochastic Search Variable Selection Prior
- [`ssvs_prior(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/ssvs_prior.bvarmodel.md)
  : Stochastic Search Variable Selection Prior
- [`ssvs_prior(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/ssvs_prior.bvecmodel.md)
  : Stochastic Search Variable Selection Prior
- [`inclusion_prior()`](https://franzmohr.github.io/bvartools/reference/inclusion_prior.md)
  : Prior Inclusion Probabilities
- [`inclusion_prior(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/inclusion_prior.bvarmodel.md)
  : Prior Inclusion Probabilities
- [`inclusion_prior(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/inclusion_prior.bvecmodel.md)
  : Prior Inclusion Probabilities
- [`cointspace_prior()`](https://franzmohr.github.io/bvartools/reference/cointspace_prior.md)
  : Build the Prior on the Cointegration Space
- [`add_sign_restrictions()`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.md)
  : Sign Restrictions
- [`add_sign_restrictions(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md)
  [`add_sign_restrictions(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md)
  [`add_sign_restrictions(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/add_sign_restrictions.bvarmodel.md)
  : Sign Restrictions
- [`add_initial_values()`](https://franzmohr.github.io/bvartools/reference/add_initial_values.md)
  : Add Initial Values of an MCMC Chain
- [`add_initial_values(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvarmodel.md)
  : Add Initial Values of an MCMC Chain
- [`add_initial_values(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_initial_values.bvecmodel.md)
  : Add Initial Values of an MCMC Chain
- [`add_initial_values(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/add_initial_values.expandingwindow.md)
  : Add Initial Values of an MCMC Chain
- [`add_initial_values(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/add_initial_values.modellist.md)
  : Add Initial Values of an MCMC Chain

## Posterior simulation

- [`add_seed()`](https://franzmohr.github.io/bvartools/reference/add_seed.md)
  : Seed of the Posterior Simulation
- [`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_priors(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_initial_values(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_seed(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_coefficients(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_loglik(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_forecast_input(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_forecasts(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`minnesota_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`ssvs_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`inclusion_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`thin(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`print(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  : Objects for Externally Produced Forecasts
- [`add_priors(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_initial_values(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_seed(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_posterior_coefficients(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_posterior_loglik(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`thin(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`selection_criteria(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  : Estimation Steps on a Folder of Models
- [`add_posterior_coefficients()`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.md)
  : Posterior Simulation of Model Coefficients
- [`add_posterior_coefficients(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvarmodel.md)
  : Posterior Simulation of Model Coefficients
- [`add_posterior_coefficients(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.bvecmodel.md)
  : Posterior Simulation of Model Coefficients
- [`add_posterior_coefficients(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.expandingwindow.md)
  : Posterior Simulation of Model Coefficients
- [`add_posterior_coefficients(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_coefficients.modellist.md)
  : Posterior Simulation of Model Coefficients
- [`add_posterior_loglik()`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.md)
  : Add Log-Likelihood
- [`add_posterior_loglik(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvarmodel.md)
  : Add Log-Likelihood
- [`add_posterior_loglik(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.bvecmodel.md)
  : Add Log-Likelihood
- [`add_posterior_loglik(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.expandingwindow.md)
  : Add Log-Likelihood
- [`add_posterior_loglik(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_loglik.modellist.md)
  : Add Log-Likelihood
- [`chain_diagnostics()`](https://franzmohr.github.io/bvartools/reference/chain_diagnostics.md)
  : Convergence Diagnostics of Several Chains
- [`bayests_posterior()`](https://franzmohr.github.io/bvartools/reference/bayests_posterior.md)
  : Posterior Simulation with the BayesTS Executable
- [`bayests_files()`](https://franzmohr.github.io/bvartools/reference/bayests_files.md)
  : Run BayesTS on Stored Models
- [`thin(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/thin.bvarmodel.md)
  : Thinning Posterior Draws
- [`thin(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/thin.bvecmodel.md)
  : Thinning Posterior Draws
- [`thin(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/thin.expandingwindow.md)
  : Thinning Posterior Draws
- [`thin(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/thin.modellist.md)
  : Thinning Posterior Draws

## Forecasting and model evaluation

- [`add_forecast_input()`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.md)
  : Add Forecast Input Data
- [`add_forecast_input(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvarmodel.md)
  : Add Forecast Input Data
- [`add_forecast_input(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.bvecmodel.md)
  : Add Forecast Input Data
- [`add_forecast_input(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.expandingwindow.md)
  : Add Forecast Input Data
- [`add_forecast_input(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/add_forecast_input.modellist.md)
  : Add Forecast Input Data
- [`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_priors(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_initial_values(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_seed(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_coefficients(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_loglik(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_forecast_input(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_forecasts(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`minnesota_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`ssvs_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`inclusion_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`thin(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`print(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  : Objects for Externally Produced Forecasts
- [`prepare_forecast_input()`](https://franzmohr.github.io/bvartools/reference/prepare_forecast_input.md)
  : Prepare Forecast Input
- [`prepare_forecast_input(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/prepare_forecast_input.bvarmodel.md)
  : Prepare Forecast Input
- [`prepare_forecast_input(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/prepare_forecast_input.bvecmodel.md)
  : Prepare Forecast Input
- [`add_posterior_forecasts()`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.md)
  : Add Forecasts
- [`add_posterior_forecasts(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvarmodel.md)
  : Add Forecasts
- [`add_posterior_forecasts(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.bvecmodel.md)
  : Add Forecasts
- [`add_posterior_forecasts(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.expandingwindow.md)
  : Add Forecasts
- [`add_posterior_forecasts(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/add_posterior_forecasts.modellist.md)
  : Add Forecasts
- [`predict(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/predict.bvarmodel.md)
  : Predict Method for Objects of Class bvar
- [`predict(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/predict.bvecmodel.md)
  : Predict Method for Objects of Class bvecmodel
- [`predict(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/predict.expandingwindow.md)
  : Predict Method for Objects of Class expandingwindow
- [`predict(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/predict.modellist.md)
  : Predict Method for Objects of Class modellist
- [`add_forecast_errors()`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.md)
  : Add Forecast Errors
- [`add_forecast_errors(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvarmodel.md)
  : Add Forecast Errors
- [`add_forecast_errors(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.bvecmodel.md)
  : Add Forecast Errors
- [`add_forecast_errors(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.expandingwindow.md)
  : Add Forecast Errors
- [`add_forecast_errors(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/add_forecast_errors.modellist.md)
  : Add Forecast Errors
- [`add_predictive_loglik()`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.md)
  : Add Predictive Log-Likelihood
- [`add_predictive_loglik(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvarmodel.md)
  : Add the Log Predictive Density of a Forecast
- [`add_predictive_loglik(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/add_predictive_loglik.bvecmodel.md)
  : Add the Log Predictive Density of a Forecast
- [`get_forecast_errors()`](https://franzmohr.github.io/bvartools/reference/get_forecast_errors.md)
  : Get Forecast Errors
- [`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md)
  : Plotting Forecast Errors per Period
- [`aggregate_forecasts()`](https://franzmohr.github.io/bvartools/reference/aggregate_forecasts.md)
  : Aggregate Forecasts to Annual Figures
- [`add_priors(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_initial_values(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_seed(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_posterior_coefficients(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_posterior_loglik(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`thin(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`selection_criteria(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  : Estimation Steps on a Folder of Models
- [`selection_criteria()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  [`print(`*`<selcrit>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  [`print(`*`<selcritlist>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  : Model Selection Criteria
- [`selection_criteria(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvarmodel.md)
  : Model Selection Criteria
- [`selection_criteria(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.bvecmodel.md)
  : Model Selection Criteria
- [`selection_criteria(`*`<default>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.default.md)
  : Selection Criteria
- [`selection_criteria(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.expandingwindow.md)
  : Model Selection Criteria
- [`selection_criteria(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.externalforecast.md)
  : Model Selection Criteria
- [`selection_criteria(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.modellist.md)
  : Model Selection Criteria
- [`time_variation_test()`](https://franzmohr.github.io/bvartools/reference/time_variation_test.md)
  : Test for Time Variation
- [`time_variation_test(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/time_variation_test.bvarmodel.md)
  [`time_variation_test(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/time_variation_test.bvarmodel.md)
  : Test for Time Variation in a VAR or VEC Model
- [`choose_best_model()`](https://franzmohr.github.io/bvartools/reference/choose_best_model.md)
  : Choose Best Model
- [`choose_best_model(`*`<selcritlist>`*`)`](https://franzmohr.github.io/bvartools/reference/choose_best_model.selcritlist.md)
  : Choose Best Model
- [`get_model_specifications()`](https://franzmohr.github.io/bvartools/reference/get_model_specifications.md)
  : Get Model Specifications

## Structural analysis

- [`irf(`*`<bvarfile>`*`)`](https://franzmohr.github.io/bvartools/reference/analysis_of_stored_models.md)
  [`fevd(`*`<bvarfile>`*`)`](https://franzmohr.github.io/bvartools/reference/analysis_of_stored_models.md)
  : Impulse Responses and Variance Decompositions of a Stored Model
- [`irf()`](https://franzmohr.github.io/bvartools/reference/irf.md) :
  Impulse Response Function
- [`irf(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/irf.bvarmodel.md)
  : Impulse Response Function
- [`irf(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/irf.bvecmodel.md)
  : Impulse Response Function
- [`fevd()`](https://franzmohr.github.io/bvartools/reference/fevd.md) :
  Forecast Error Variance Decomposition
- [`fevd(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/fevd.bvarmodel.md)
  : Forecast Error Variance Decomposition
- [`fevd(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/fevd.bvecmodel.md)
  : Forecast Error Variance Decomposition
- [`spillover()`](https://franzmohr.github.io/bvartools/reference/spillover.md)
  : Spillover Index
- [`spillover(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/spillover.bvarmodel.md)
  : Spillover Index
- [`spillover(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/spillover.bvecmodel.md)
  : Spillover Index
- [`spillover(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/spillover.expandingwindow.md)
  : Rolling Spillover Index
- [`spillover(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/spillover.modellist.md)
  : Spillover Index
- [`multipliers()`](https://franzmohr.github.io/bvartools/reference/multipliers.md)
  : Dynamic Multipliers
- [`multipliers(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/multipliers.bvarmodel.md)
  : Dynamic Multipliers of a VAR Model with Exogenous Variables
- [`multipliers(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/multipliers.bvecmodel.md)
  : Dynamic Multipliers of a VEC Model with Exogenous Variables
- [`vec_to_var()`](https://franzmohr.github.io/bvartools/reference/vec_to_var.md)
  : Transform a VEC Model to a VAR in Levels
- [`vec_to_var(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/vec_to_var.bvecmodel.md)
  : Transform a VEC Model to a VAR in Levels
- [`vec_to_var(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/vec_to_var.expandingwindow.md)
  : Transform VEC Models to VARs in Levels
- [`vec_to_var(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/vec_to_var.modellist.md)
  : Transform a VEC Model to a VAR in Levels
- [`scale_error_correction()`](https://franzmohr.github.io/bvartools/reference/scale_error_correction.md)
  : Scale Error Correction
- [`scale_error_correction(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/scale_error_correction.bvecmodel.md)
  : Scale Error Correction
- [`rescale_error_correction()`](https://franzmohr.github.io/bvartools/reference/rescale_error_correction.md)
  : Rescale Error Correction
- [`rescale_error_correction(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/rescale_error_correction.bvecmodel.md)
  : Rescale Error Correction

## Summaries and plots

- [`summary(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/summary.bvarmodel.md)
  [`print(`*`<summary.bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/summary.bvarmodel.md)
  : Summarising Bayesian VAR Coefficients
- [`summary(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/summary.bvecmodel.md)
  [`print(`*`<summary.bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/summary.bvecmodel.md)
  : Summarising Bayesian VEC Coefficients
- [`summary(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/summary.expandingwindow.md)
  : Summarising Bayesian Models
- [`summary(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/summary.modellist.md)
  : Summarising Bayesian Models
- [`create_external_forecast()`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_priors(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_initial_values(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_seed(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_coefficients(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_loglik(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_forecast_input(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`add_posterior_forecasts(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`minnesota_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`ssvs_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`inclusion_prior(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`thin(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  [`print(`*`<externalforecast>`*`)`](https://franzmohr.github.io/bvartools/reference/create_external_forecast.md)
  : Objects for Externally Produced Forecasts
- [`open_model()`](https://franzmohr.github.io/bvartools/reference/open_model.md)
  [`print(`*`<bvarfile>`*`)`](https://franzmohr.github.io/bvartools/reference/open_model.md)
  : Open a Model Stored in an HDF5 File
- [`open_models()`](https://franzmohr.github.io/bvartools/reference/open_models.md)
  [`print(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/open_models.md)
  [`model_files()`](https://franzmohr.github.io/bvartools/reference/open_models.md)
  : Open a Folder of Stored Models
- [`print(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/print.bvarmodel.md)
  : Printing Model Information
- [`print(`*`<bvarspillover>`*`)`](https://franzmohr.github.io/bvartools/reference/print.bvarspillover.md)
  : Print a Spillover Index
- [`print(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/print.bvecmodel.md)
  : Printing Model Information
- [`print(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/print.modellist.md)
  : Printing Model Information
- [`selection_criteria()`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  [`print(`*`<selcrit>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  [`print(`*`<selcritlist>`*`)`](https://franzmohr.github.io/bvartools/reference/selection_criteria.md)
  : Model Selection Criteria
- [`plot(`*`<bvarfevd>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.bvarfevd.md)
  : Plotting Forecast Error Variance Decompositions of Bayesian Vector
  Autoregression
- [`plot(`*`<bvarirf>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.bvarirf.md)
  : Plotting Impulse Responses of Bayesian Vector Autoregression
- [`plot(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.bvarmodel.md)
  : Plotting Draws of a Bayesian VAR Model
- [`plot(`*`<bvarprd>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.bvarprd.md)
  : Plotting Forecasts of BVAR Models
- [`plot(`*`<bvarprdlist>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.bvarprdlist.md)
  : Plotting Forecasts of a List of Bayesian VAR Models
- [`plot(`*`<bvarspillover>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.bvarspillover.md)
  : Plot a Spillover Index
- [`plot(`*`<bvarspilloverts>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.bvarspilloverts.md)
  : Plot a Rolling Spillover Index
- [`plot(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.bvecmodel.md)
  : Plotting Draws of a Bayesian VEC Model
- [`plot(`*`<expandwindbvarprdlist>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.expandwindbvarprdlist.md)
  : Plotting Forecasts of a List of Bayesian VAR Models
- [`plot(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.modellist.md)
  : Plotting Draws of a Bayesian Time Series Models
- [`plot(`*`<selcritlist>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.selcritlist.md)
  [`plot(`*`<selcrit>`*`)`](https://franzmohr.github.io/bvartools/reference/plot.selcritlist.md)
  : Plotting Selection Criteria
- [`plot_forecast_errors_by_period()`](https://franzmohr.github.io/bvartools/reference/plot_forecast_errors_by_period.md)
  : Plotting Forecast Errors per Period

## Storage

- [`write_to_hdf5()`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.md)
  : Export to HDF5 File
- [`write_to_hdf5(`*`<bvarmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.bvarmodel.md)
  : Export to HDF5 File
- [`write_to_hdf5(`*`<bvecmodel>`*`)`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.bvecmodel.md)
  : Export to HDF5 File
- [`write_to_hdf5(`*`<expandingwindow>`*`)`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.expandingwindow.md)
  : Export to HDF5 File
- [`write_to_hdf5(`*`<modellist>`*`)`](https://franzmohr.github.io/bvartools/reference/write_to_hdf5.modellist.md)
  : Export to HDF5 File
- [`read_model_from_hdf5()`](https://franzmohr.github.io/bvartools/reference/read_model_from_hdf5.md)
  : Import Models from HDF5 Files
- [`read_models_from_folder()`](https://franzmohr.github.io/bvartools/reference/read_models_from_folder.md)
  : Import Models from a Folder of HDF5 Files
- [`read_expanding_window_model_from_folder()`](https://franzmohr.github.io/bvartools/reference/read_expanding_window_model_from_folder.md)
  : Import Models from HDF5 Files
- [`list_models_in_hdf5()`](https://franzmohr.github.io/bvartools/reference/list_models_in_hdf5.md)
  : Models in an HDF5 File
- [`open_model()`](https://franzmohr.github.io/bvartools/reference/open_model.md)
  [`print(`*`<bvarfile>`*`)`](https://franzmohr.github.io/bvartools/reference/open_model.md)
  : Open a Model Stored in an HDF5 File
- [`open_models()`](https://franzmohr.github.io/bvartools/reference/open_models.md)
  [`print(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/open_models.md)
  [`model_files()`](https://franzmohr.github.io/bvartools/reference/open_models.md)
  : Open a Folder of Stored Models
- [`add_priors(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_initial_values(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_seed(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_posterior_coefficients(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`add_posterior_loglik(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`thin(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  [`selection_criteria(`*`<bvarfolder>`*`)`](https://franzmohr.github.io/bvartools/reference/folder_steps.md)
  : Estimation Steps on a Folder of Models
- [`irf(`*`<bvarfile>`*`)`](https://franzmohr.github.io/bvartools/reference/analysis_of_stored_models.md)
  [`fevd(`*`<bvarfile>`*`)`](https://franzmohr.github.io/bvartools/reference/analysis_of_stored_models.md)
  : Impulse Responses and Variance Decompositions of a Stored Model
- [`map_models()`](https://franzmohr.github.io/bvartools/reference/map_models.md)
  : Apply a Function to the Models of a Folder
- [`map_draws()`](https://franzmohr.github.io/bvartools/reference/map_draws.md)
  : Apply a Function to the Draws of a Stored Model

## Building blocks for custom samplers

Posterior draws and data preparation steps for writing your own Gibbs
sampler, and the constructors that collect its output in a standard
object.

- [`bvar()`](https://franzmohr.github.io/bvartools/reference/bvar.md) :
  Bayesian Vector Autoregression Objects
- [`bvec()`](https://franzmohr.github.io/bvartools/reference/bvec.md) :
  Bayesian Vector Error Correction Objects
- [`post_bvs()`](https://franzmohr.github.io/bvartools/reference/post_bvs.md)
  : Bayesian Variable Selection
- [`post_coint_kls()`](https://franzmohr.github.io/bvartools/reference/post_coint_kls.md)
  : Posterior Draw for Cointegration Models
- [`post_coint_kls_sur()`](https://franzmohr.github.io/bvartools/reference/post_coint_kls_sur.md)
  : Posterior Draw for Cointegration Models
- [`post_gamma_measurement_variance()`](https://franzmohr.github.io/bvartools/reference/post_gamma_measurement_variance.md)
  : Posterior Draws of Error Variances
- [`post_gamma_state_variance()`](https://franzmohr.github.io/bvartools/reference/post_gamma_state_variance.md)
  : Posterior Draws of Error Variances
- [`post_normal()`](https://franzmohr.github.io/bvartools/reference/post_normal.md)
  : Posterior Draw from a Normal Distribution
- [`post_normal_sur()`](https://franzmohr.github.io/bvartools/reference/post_normal_sur.md)
  : Posterior Draw from a Normal Distribution
- [`ssvs()`](https://franzmohr.github.io/bvartools/reference/ssvs.md) :
  Stochastic Search Variable Selection
- [`stochvol_ksc_1998()`](https://franzmohr.github.io/bvartools/reference/stochvol_ksc_1998.md)
  : Stochastic Volatility
- [`stochvol_ocsn_2007()`](https://franzmohr.github.io/bvartools/reference/stochvol_ocsn_2007.md)
  : Stochastic Volatility
- [`kalman_durbin_koopman_2002()`](https://franzmohr.github.io/bvartools/reference/kalman_durbin_koopman_2002.md)
  : Durbin and Koopman Simulation Smoother
- [`covar_prepare_data()`](https://franzmohr.github.io/bvartools/reference/covar_prepare_data.md)
  : Covariance: Data Preparation
- [`covar_vector_to_matrix()`](https://franzmohr.github.io/bvartools/reference/covar_vector_to_matrix.md)
  : Covariance: Vector to Matrix
- [`coint_kls2010_reparameterise_two()`](https://franzmohr.github.io/bvartools/reference/coint_kls2010_reparameterise_two.md)
  : Cointegration Reparameterisation
- [`coint_prepare_sur_data()`](https://franzmohr.github.io/bvartools/reference/coint_prepare_sur_data.md)
  : Posterior Data Preparation
- [`generate_lower_block_diagonal()`](https://franzmohr.github.io/bvartools/reference/generate_lower_block_diagonal.md)
  : Posterior Data Preparation
- [`sur_const_to_tvp()`](https://franzmohr.github.io/bvartools/reference/sur_const_to_tvp.md)
  : SUR Matrix Transformation
- [`loglik_normal()`](https://franzmohr.github.io/bvartools/reference/loglik_normal.md)
  : Log-Likelihood of a Multivariate Normal Distribution

## Data

- [`at_macrodata`](https://franzmohr.github.io/bvartools/reference/at_macrodata.md)
  : Austrian sub-model of a global VAR
- [`e1`](https://franzmohr.github.io/bvartools/reference/e1.md) : West
  German economic time series data
- [`e6`](https://franzmohr.github.io/bvartools/reference/e6.md) : German
  interest and inflation rate data
- [`uk_macrodata`](https://franzmohr.github.io/bvartools/reference/uk_macrodata.md)
  : UK interest and inflation rate data
- [`us_macrodata`](https://franzmohr.github.io/bvartools/reference/us_macrodata.md)
  : US macroeconomic data
