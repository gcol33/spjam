// scatR: Joint Spatial Modelling Under Preferential Sampling
// C++ initialization and registration
// Note: Rcpp::compileAttributes() generates RcppExports.cpp with actual registrations.
// This file provides the R_init hook and any non-exported helpers.

#include <Rcpp.h>

// Placeholder for future HMC sampler
// [[Rcpp::export]]
Rcpp::List scatr_hmc_placeholder() {
    return Rcpp::List::create(
        Rcpp::Named("status") = "not_implemented",
        Rcpp::Named("message") = "HMC sampler not yet implemented. Use backend = 'laplace'."
    );
}
