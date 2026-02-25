// covariance.h
// Matern covariance functions for Gaussian process spatial models
// Adapted from numdenom v1.3.0 hmc_svc.h

#ifndef SPJAM_COVARIANCE_H
#define SPJAM_COVARIANCE_H

#include <cmath>

namespace spjam {

// Exponential covariance (Matern nu = 0.5)
// C(d) = sigma^2 * exp(-d / range)
inline double cov_exponential(double d, double range, double sigma) {
  if (d == 0.0) return sigma * sigma;
  return sigma * sigma * std::exp(-d / range);
}

// Matern 3/2 covariance (nu = 1.5)
// C(d) = sigma^2 * (1 + sqrt(3)*d/range) * exp(-sqrt(3)*d/range)
inline double cov_matern32(double d, double range, double sigma) {
  if (d == 0.0) return sigma * sigma;
  double r = std::sqrt(3.0) * d / range;
  return sigma * sigma * (1.0 + r) * std::exp(-r);
}

// Matern 5/2 covariance (nu = 2.5)
// C(d) = sigma^2 * (1 + sqrt(5)*d/range + 5*d^2/(3*range^2)) * exp(-sqrt(5)*d/range)
inline double cov_matern52(double d, double range, double sigma) {
  if (d == 0.0) return sigma * sigma;
  double r = std::sqrt(5.0) * d / range;
  return sigma * sigma * (1.0 + r + r * r / 3.0) * std::exp(-r);
}

// Dispatch by nu
inline double cov_matern(double d, double range, double sigma, double nu) {
  if (nu <= 0.75) return cov_exponential(d, range, sigma);
  if (nu <= 2.0) return cov_matern32(d, range, sigma);
  return cov_matern52(d, range, sigma);
}

// Build full dense covariance matrix for GP
// coords: [n x 2], row-major
// C_out: [n x n], row-major (must be pre-allocated)
inline void build_cov_matrix(const double* coords, int n,
                              double range, double sigma, double nu,
                              double* C_out) {
  for (int i = 0; i < n; i++) {
    C_out[i * n + i] = sigma * sigma;
    for (int j = i + 1; j < n; j++) {
      double dx = coords[i * 2] - coords[j * 2];
      double dy = coords[i * 2 + 1] - coords[j * 2 + 1];
      double d = std::sqrt(dx * dx + dy * dy);
      double c = cov_matern(d, range, sigma, nu);
      C_out[i * n + j] = c;
      C_out[j * n + i] = c;
    }
  }
}

} // namespace spjam

#endif // SPJAM_COVARIANCE_H
