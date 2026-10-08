# The generics of the workflow and the infrastructure every model package
# shares live in bayests. They are imported for this package's methods and
# re-exported, so that code written against bvartools finds them where it
# always did.

#' @importFrom bayests add_forecast_errors
#' @export
bayests::add_forecast_errors

#' @importFrom bayests add_forecast_input
#' @export
bayests::add_forecast_input

#' @importFrom bayests add_initial_values
#' @export
bayests::add_initial_values

#' @importFrom bayests add_posterior_coefficients
#' @export
bayests::add_posterior_coefficients

#' @importFrom bayests add_posterior_forecasts
#' @export
bayests::add_posterior_forecasts

#' @importFrom bayests add_posterior_loglik
#' @export
bayests::add_posterior_loglik

#' @importFrom bayests add_predictive_loglik
#' @export
bayests::add_predictive_loglik

#' @importFrom bayests add_priors
#' @export
bayests::add_priors

#' @importFrom bayests add_seed
#' @export
bayests::add_seed

#' @importFrom bayests add_sign_restrictions
#' @export
bayests::add_sign_restrictions

#' @importFrom bayests add_sign_zero_restrictions
#' @export
bayests::add_sign_zero_restrictions

#' @importFrom bayests align_model_obs
#' @export
bayests::align_model_obs

#' @importFrom bayests bayests_files
#' @export
bayests::bayests_files

#' @importFrom bayests bayests_posterior
#' @export
bayests::bayests_posterior

#' @importFrom bayests chain_diagnostics
#' @export
bayests::chain_diagnostics

#' @importFrom bayests choose_best_model
#' @export
bayests::choose_best_model

#' @importFrom bayests combine_models
#' @export
bayests::combine_models

#' @importFrom bayests expected_model_size
#' @export
bayests::expected_model_size

#' @importFrom bayests fevd
#' @export
bayests::fevd

#' @importFrom bayests from_bayests_tree
#' @export
bayests::from_bayests_tree

#' @importFrom bayests read_hdf5_model
#' @export
bayests::read_hdf5_model

#' @importFrom bayests get_forecast_errors
#' @export
bayests::get_forecast_errors

#' @importFrom bayests get_model_specifications
#' @export
bayests::get_model_specifications

#' @importFrom bayests irf
#' @export
bayests::irf

#' @importFrom bayests read_model_groups
#' @export
bayests::read_model_groups

#' @importFrom bayests map_draws
#' @export
bayests::map_draws

#' @importFrom bayests map_models
#' @export
bayests::map_models

#' @importFrom bayests model_files
#' @export
bayests::model_files

#' @importFrom bayests read_model_file
#' @export
bayests::read_model_file

#' @importFrom bayests read_model_folder
#' @export
bayests::read_model_folder

#' @importFrom bayests plot_predictive_loglik
#' @export
bayests::plot_predictive_loglik

#' @importFrom bayests prepare_forecast_input
#' @export
bayests::prepare_forecast_input

#' @importFrom bayests read_bayests_tree
#' @export
bayests::read_bayests_tree

#' @importFrom bayests read_expanding_window_model_from_folder
#' @export
bayests::read_expanding_window_model_from_folder

#' @importFrom bayests read_model_from_hdf5
#' @export
bayests::read_model_from_hdf5

#' @importFrom bayests read_models_from_folder
#' @export
bayests::read_models_from_folder

#' @importFrom bayests selection_criteria
#' @export
bayests::selection_criteria

#' @importFrom bayests time_variation_test
#' @export
bayests::time_variation_test

#' @importFrom bayests use_expanding_window
#' @export
bayests::use_expanding_window

#' @importFrom bayests write_bayests_tree
#' @export
bayests::write_bayests_tree

#' @importFrom bayests write_to_hdf5
#' @export
bayests::write_to_hdf5

# What this package's methods need of bayests beyond its generics: seeds,
# the HDF5 writers and readers, and the pieces of the selection criteria.
#' @importFrom bayests draw_error_bars forecast_plot_series .check_max_groups .limit_fevd_groups .check_seed .draw_model_seed .offset_seed .draws_matrix .hdf5_draws_attrs .hdf5_series_attrs .hdf5_read_matrix .hdf5_read_value .hdf5_restore_class .hdf5_write_model .hdf5_group .hdf5_write .hdf5_write_attr .log_mean_exp .nse_log_mean_exp .lpl_entry .loo_data_frame .waic_data_frame .read_posterior_block
NULL
