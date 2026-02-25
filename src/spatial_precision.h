// spatial_precision.h
// Spatial prior precision structures (ICAR, BYM2, GP, SPDE)
// Adapted from numdenom v1.3.0 laplace_core.cpp + hmc_gp.h

#ifndef SPJAM_SPATIAL_PRECISION_H
#define SPJAM_SPATIAL_PRECISION_H

#include <vector>
#include <cmath>

namespace spjam {

// =====================================================================
// ICAR prior: x' Q x where Q = diag(n_neighbors) - Adjacency
// =====================================================================

// Compute ICAR quadratic form: 0.5 * tau * sum_{i~j} (x_i - x_j)^2
inline double icar_quadratic(const double* x, int n,
                              const int* adj_row_ptr,
                              const int* adj_col_idx,
                              double tau) {
  double qf = 0.0;
  for (int i = 0; i < n; i++) {
    for (int k = adj_row_ptr[i]; k < adj_row_ptr[i + 1]; k++) {
      int j = adj_col_idx[k];
      if (j > i) {  // Count each edge once
        double diff = x[i] - x[j];
        qf += diff * diff;
      }
    }
  }
  return 0.5 * tau * qf;
}

// ICAR gradient: d/dx_i [-0.5 * tau * x'Qx] = -tau * (n_i * x_i - sum_{j~i} x_j)
inline void icar_gradient(const double* x, double* grad, int n,
                           const int* adj_row_ptr,
                           const int* adj_col_idx,
                           const int* n_neighbors,
                           double tau) {
  for (int i = 0; i < n; i++) {
    double sum_neighbors = 0.0;
    for (int k = adj_row_ptr[i]; k < adj_row_ptr[i + 1]; k++) {
      sum_neighbors += x[adj_col_idx[k]];
    }
    grad[i] -= tau * (n_neighbors[i] * x[i] - sum_neighbors);
  }
}

// ICAR Hessian diagonal contribution: -d^2/dx_i^2 = tau * n_neighbors[i]
// Off-diagonal i,j: -tau (if i~j)
inline void icar_hessian_diag(double* H_diag, int n,
                               const int* n_neighbors,
                               double tau) {
  for (int i = 0; i < n; i++) {
    H_diag[i] += tau * n_neighbors[i];
  }
}

// =====================================================================
// BYM2 prior (Riebler parameterization)
// x = sigma * (sqrt(rho) * phi_structured + sqrt(1-rho) * theta_iid)
// where phi_structured is ICAR-scaled
// =====================================================================

// BYM2 quadratic form on the structured component
// scale_factor: geometric mean of marginal variances of ICAR
inline double bym2_quadratic(const double* phi_struct, const double* theta_iid,
                              int n, const int* adj_row_ptr,
                              const int* adj_col_idx,
                              double sigma, double rho,
                              double scale_factor) {
  // Structured part: ICAR precision on phi_struct
  // scaled by 1/(sigma^2 * rho * scale_factor)
  double tau_struct = 1.0 / (sigma * sigma * rho * scale_factor + 1e-10);
  double qf_struct = 0.0;
  for (int i = 0; i < n; i++) {
    for (int k = adj_row_ptr[i]; k < adj_row_ptr[i + 1]; k++) {
      int j = adj_col_idx[k];
      if (j > i) {
        double diff = phi_struct[i] - phi_struct[j];
        qf_struct += diff * diff;
      }
    }
  }

  // IID part: precision 1/(sigma^2 * (1-rho))
  double tau_iid = 1.0 / (sigma * sigma * (1.0 - rho) + 1e-10);
  double qf_iid = 0.0;
  for (int i = 0; i < n; i++) {
    qf_iid += theta_iid[i] * theta_iid[i];
  }

  return 0.5 * tau_struct * qf_struct + 0.5 * tau_iid * qf_iid;
}

// =====================================================================
// GP prior: x ~ N(0, C(range, sigma))
// Precision Q = C^{-1}
// For dense GP: -0.5 * x' * Q * x - 0.5 * log|C|
// =====================================================================

// GP log-prior using dense covariance (for moderate n_spatial)
// Requires pre-computed Cholesky L of covariance matrix C
// log p(x) = -0.5 * x' C^{-1} x - 0.5 * log|C| - n/2 * log(2pi)
// = -0.5 * ||L^{-1} x||^2 - sum(log(diag(L))) - n/2 * log(2pi)
inline double gp_log_prior_dense(const double* x, const double* L_flat,
                                  int n) {
  // Forward solve: z = L^{-1} x
  std::vector<double> z(n);
  for (int i = 0; i < n; i++) {
    double sum = x[i];
    for (int k = 0; k < i; k++) {
      sum -= L_flat[i * n + k] * z[k];
    }
    z[i] = sum / L_flat[i * n + i];
  }

  // ||z||^2
  double qf = 0.0;
  double log_det = 0.0;
  for (int i = 0; i < n; i++) {
    qf += z[i] * z[i];
    log_det += std::log(L_flat[i * n + i]);
  }

  return -0.5 * qf - log_det - 0.5 * n * std::log(2.0 * M_PI);
}

// GP gradient: d/dx [-0.5 x' Q x] = -Q x = -C^{-1} x
// Computed via Cholesky: grad = -L'^{-1} L^{-1} x
inline void gp_gradient_dense(const double* x, const double* L_flat,
                               double* grad, int n) {
  // Forward solve: z = L^{-1} x
  std::vector<double> z(n);
  for (int i = 0; i < n; i++) {
    double sum = x[i];
    for (int k = 0; k < i; k++) {
      sum -= L_flat[i * n + k] * z[k];
    }
    z[i] = sum / L_flat[i * n + i];
  }

  // Back solve: q = L'^{-1} z (= C^{-1} x)
  std::vector<double> q(n);
  for (int i = n - 1; i >= 0; i--) {
    double sum = z[i];
    for (int k = i + 1; k < n; k++) {
      sum -= L_flat[k * n + i] * q[k];
    }
    q[i] = sum / L_flat[i * n + i];
  }

  // grad -= Q * x
  for (int i = 0; i < n; i++) {
    grad[i] -= q[i];
  }
}

// GP Hessian contribution: the precision matrix Q = C^{-1}
// For dense case, add Q to the Hessian
inline void gp_hessian_dense(const double* L_flat, double* H_flat,
                              int n, int H_stride, int offset) {
  // Compute Q = L'^{-1} L^{-1} column by column
  // Q_flat stored separately, then added to H at offset
  std::vector<double> e(n, 0.0);
  std::vector<double> z(n);
  std::vector<double> q(n);

  for (int col = 0; col < n; col++) {
    // Unit vector
    std::fill(e.begin(), e.end(), 0.0);
    e[col] = 1.0;

    // Forward solve: z = L^{-1} e
    for (int i = 0; i < n; i++) {
      double sum = e[i];
      for (int k = 0; k < i; k++) {
        sum -= L_flat[i * n + k] * z[k];
      }
      z[i] = sum / L_flat[i * n + i];
    }

    // Back solve: q = L'^{-1} z
    for (int i = n - 1; i >= 0; i--) {
      double sum = z[i];
      for (int k = i + 1; k < n; k++) {
        sum -= L_flat[k * n + i] * q[k];
      }
      q[i] = sum / L_flat[i * n + i];
    }

    // Add column of Q^{-1} to H at [offset, offset]
    for (int row = 0; row < n; row++) {
      H_flat[(offset + row) * H_stride + (offset + col)] += q[row];
    }
  }
}

// =====================================================================
// SPDE precision: Q(kappa, tau) = tau^2 * (kappa^4 * C + 2*kappa^2 * G + G * C^{-1} * G)
// For alpha=2 (Matern nu=1): Q = tau^2 * (kappa^4 * C + 2*kappa^2 * G1 + G2)
// Passed as sparse matrices from R
// =====================================================================

// SPDE quadratic form: 0.5 * x' Q x
// Q stored as sparse CSR
inline double spde_quadratic(const double* x, int n,
                              const int* Q_row_ptr,
                              const int* Q_col_idx,
                              const double* Q_vals) {
  double qf = 0.0;
  for (int i = 0; i < n; i++) {
    for (int k = Q_row_ptr[i]; k < Q_row_ptr[i + 1]; k++) {
      int j = Q_col_idx[k];
      if (j == i) {
        qf += Q_vals[k] * x[i] * x[i];
      } else if (j > i) {
        qf += 2.0 * Q_vals[k] * x[i] * x[j];
      }
    }
  }
  return 0.5 * qf;
}

// SPDE gradient: -Q * x
inline void spde_gradient(const double* x, double* grad, int n,
                           const int* Q_row_ptr,
                           const int* Q_col_idx,
                           const double* Q_vals) {
  for (int i = 0; i < n; i++) {
    double Qx_i = 0.0;
    for (int k = Q_row_ptr[i]; k < Q_row_ptr[i + 1]; k++) {
      Qx_i += Q_vals[k] * x[Q_col_idx[k]];
    }
    grad[i] -= Qx_i;
  }
}

} // namespace spjam

#endif // SPJAM_SPATIAL_PRECISION_H
