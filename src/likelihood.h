// likelihood.h
// Log-likelihood functions with gradients and Hessians
// Adapted from numdenom v1.3.0 laplace_core.cpp lines 24-106

#ifndef SCATR_LIKELIHOOD_H
#define SCATR_LIKELIHOOD_H

#include "linalg_utils.h"
#include <cmath>

namespace scatr {

// =====================================================================
// Poisson: y ~ Poisson(mu = exp(eta))
// =====================================================================

inline double log_lik_poisson(int y, double eta) {
  return y * eta - scatr_linalg::safe_exp(eta) - R::lgammafn(y + 1.0);
}

inline double grad_log_lik_poisson(int y, double eta) {
  return y - scatr_linalg::safe_exp(eta);
}

inline double neg_hess_log_lik_poisson(int y, double eta) {
  return scatr_linalg::safe_exp(eta);
}

// =====================================================================
// Negative binomial: y ~ NegBin(mu = exp(eta), phi)
// Parameterization: E[Y] = mu, Var[Y] = mu + mu^2/phi
// =====================================================================

inline double log_lik_negbin(int y, double eta, double phi) {
  double mu = scatr_linalg::safe_exp(eta);
  return R::lgammafn(y + phi) - R::lgammafn(phi) - R::lgammafn(y + 1.0)
       + phi * std::log(phi / (mu + phi))
       + y * std::log(mu / (mu + phi));
}

inline double grad_log_lik_negbin(int y, double eta, double phi) {
  double mu = scatr_linalg::safe_exp(eta);
  double p = mu / (mu + phi);
  return y - (y + phi) * p;
}

inline double neg_hess_log_lik_negbin(int y, double eta, double phi) {
  double mu = scatr_linalg::safe_exp(eta);
  double denom = mu + phi;
  return (y + phi) * mu * phi / (denom * denom);
}

// =====================================================================
// Binomial: y ~ Binomial(n, logit^{-1}(eta))
// =====================================================================

inline double log_lik_binomial(int y, int n, double eta) {
  double log_p;
  if (eta > 0) {
    log_p = y * eta - n * eta - n * std::log(1.0 + std::exp(-eta));
  } else {
    log_p = y * eta - n * std::log(1.0 + std::exp(eta));
  }
  return log_p;
}

inline double grad_log_lik_binomial(int y, int n, double eta) {
  double p;
  if (eta > 0) {
    p = 1.0 / (1.0 + std::exp(-eta));
  } else {
    double exp_eta = std::exp(eta);
    p = exp_eta / (1.0 + exp_eta);
  }
  return y - n * p;
}

inline double neg_hess_log_lik_binomial(int y, int n, double eta) {
  double p;
  if (eta > 0) {
    p = 1.0 / (1.0 + std::exp(-eta));
  } else {
    double exp_eta = std::exp(eta);
    p = exp_eta / (1.0 + exp_eta);
  }
  return n * p * (1.0 - p);
}

// =====================================================================
// Dispatch by family string
// =====================================================================

inline double compute_grad(const std::string& family, int y, int n_trials,
                            double eta, double phi) {
  if (family == "poisson") return grad_log_lik_poisson(y, eta);
  if (family == "negbin") return grad_log_lik_negbin(y, eta, phi);
  if (family == "binomial") return grad_log_lik_binomial(y, n_trials, eta);
  return 0.0;
}

inline double compute_neg_hess(const std::string& family, int y, int n_trials,
                                double eta, double phi) {
  if (family == "poisson") return neg_hess_log_lik_poisson(y, eta);
  if (family == "negbin") return neg_hess_log_lik_negbin(y, eta, phi);
  if (family == "binomial") return neg_hess_log_lik_binomial(y, n_trials, eta);
  return 0.0;
}

inline double compute_log_lik(const std::string& family, int y, int n_trials,
                               double eta, double phi) {
  if (family == "poisson") return log_lik_poisson(y, eta);
  if (family == "negbin") return log_lik_negbin(y, eta, phi);
  if (family == "binomial") return log_lik_binomial(y, n_trials, eta);
  return 0.0;
}

} // namespace scatr

#endif // SCATR_LIKELIHOOD_H
