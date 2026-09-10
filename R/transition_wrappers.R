# Transition announcements for the functions that are implemented in C++.
#
# The message cannot be emitted from the C++ body, and the compiled routines
# must keep both their names and their registered symbols: they carry
# '[[Rcpp::interfaces(r, cpp)]]', so they are part of the C++ interface in
# 'inst/include' that other packages link to, and the posterior simulation
# algorithms in 'src' call them directly.
#
# The R interface is therefore wrapped instead. This file is collated after
# 'RcppExports.R', so each compiled function is first bound to an internal name
# and then shadowed by a wrapper that announces the change before forwarding to
# it. Callers from C++ are unaffected, since they reach the compiled routine
# through its registered symbol rather than through these bindings.

.kalman_dk_cpp <- kalman_dk
.stochvol_ksc1998_cpp <- stochvol_ksc1998
.stochvol_ocsn2007_cpp <- stochvol_ocsn2007
.bvs_cpp <- bvs
.stoch_vol_cpp <- stoch_vol

kalman_dk <- function(y, z, sigma_u, sigma_v, B, a_init, P_init) {
  .transition_message("kalman_dk", "kalman_durbin_koopman_2002")

  .kalman_dk_cpp(y, z, sigma_u, sigma_v, B, a_init, P_init)
}

stochvol_ksc1998 <- function(y, h, sigma, h_init, constant) {
  .transition_message("stochvol_ksc1998", "stochvol_ksc_1998")

  .stochvol_ksc1998_cpp(y, h, sigma, h_init, constant)
}

stochvol_ocsn2007 <- function(y, h, sigma, h_init, constant) {
  .transition_message("stochvol_ocsn2007", "stochvol_ocsn_2007")

  .stochvol_ocsn2007_cpp(y, h, sigma, h_init, constant)
}

bvs <- function(y, z, a, lambda, sigma_i, prob_prior, include = NULL) {
  .transition_message("bvs", "post_bvs")

  .bvs_cpp(y, z, a, lambda, sigma_i, prob_prior, include)
}

stoch_vol <- function(y, h, sigma, h_init, constant) {
  .transition_message("stoch_vol", "stochvol_ksc_1998",
                      note = paste("'stoch_vol()' is a wrapper for the algorithm of Kim, Shephard",
                                   "and Chib (1998), which 'stochvol_ksc_1998()' implements directly."))

  .stoch_vol_cpp(y, h, sigma, h_init, constant)
}
