// laplace_joint.cpp
// Joint Laplace approximation engine for preferential sampling models
// Adapted from numdenom v1.3.0 laplace_core.cpp (spatial Newton-Raphson)

#include "laplace_joint.h"
#include "likelihood.h"
#include "linalg_utils.h"
#include "covariance.h"
#include "spatial_precision.h"
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

namespace spjam {

JointLaplaceResult joint_laplace_mode(
    const std::vector<int>& y,
    const std::vector<int>& n_trials,
    const std::vector<double>& X_eco_flat,
    const std::string& family,
    double phi_overdispersion,
    const std::vector<double>& Z_samp_flat,
    const std::vector<double>& Z_quad_flat,
    const std::vector<double>& quad_weights,
    const std::vector<int>& quad_spatial_idx,
    int spatial_type,
    const std::vector<int>& obs_spatial_idx,
    const std::vector<double>& gp_L_shared_flat,
    const std::vector<double>& gp_L_eco_flat,
    const std::vector<double>& gp_L_samp_flat,
    const std::vector<int>& adj_row_ptr,
    const std::vector<int>& adj_col_idx,
    const std::vector<int>& n_neighbors,
    double tau_shared,
    double tau_eco,
    double tau_samp,
    const std::vector<int>& spde_Q_row_ptr,
    const std::vector<int>& spde_Q_col_idx,
    const std::vector<double>& spde_Q_shared_vals,
    const std::vector<double>& spde_Q_eco_vals,
    const std::vector<double>& spde_Q_samp_vals,
    int latent_type,
    double beta_share,
    const JointParamLayout& layout,
    double tau_beta,
    int max_iter,
    double tol
) {
  int N = y.size();
  int N_quad = quad_weights.size();
  int n_x = layout.total_params;
  int p_eco = layout.p_eco;
  int p_samp = layout.p_samp;
  int n_sp = layout.n_spatial;
  bool has_shared = (latent_type == LATENT_SHARED);

  JointLaplaceResult result;
  result.converged = false;
  result.n_iter = 0;

  // Initialize latent vector to zero
  std::vector<double> x(n_x, 0.0);
  std::vector<double> grad(n_x, 0.0);
  std::vector<double> H_flat(n_x * n_x, 0.0);

  for (int iter = 0; iter < max_iter; iter++) {
    // Zero gradient and Hessian
    std::fill(grad.begin(), grad.end(), 0.0);
    std::fill(H_flat.begin(), H_flat.end(), 0.0);

    // Pointers to parameter blocks
    const double* beta_eco = &x[layout.beta_eco_start];
    const double* beta_samp = &x[layout.beta_samp_start];
    const double* S_shared = has_shared ? &x[layout.s_shared_start] : nullptr;
    const double* S_eco = &x[layout.s_eco_start];
    const double* S_samp = &x[layout.s_samp_start];

    // ==================================================================
    // 1. ECOLOGICAL LIKELIHOOD contributions
    // ==================================================================
    for (int i = 0; i < N; i++) {
      int si = obs_spatial_idx[i];

      // Linear predictor: eta_eco = X_eco * beta_eco + S_shared[si] + S_eco[si]
      double eta_eco = spjam_linalg::dot_product(
        &X_eco_flat[i * p_eco], beta_eco, p_eco);
      if (has_shared) eta_eco += S_shared[si];
      eta_eco += S_eco[si];

      double g_i = compute_grad(family, y[i], n_trials[i], eta_eco,
                                 phi_overdispersion);
      double h_i = compute_neg_hess(family, y[i], n_trials[i], eta_eco,
                                     phi_overdispersion);

      // Gradient w.r.t. beta_eco
      for (int j = 0; j < p_eco; j++) {
        grad[layout.beta_eco_start + j] += g_i * X_eco_flat[i * p_eco + j];
      }

      // Gradient w.r.t. S_shared (if shared)
      if (has_shared) {
        grad[layout.s_shared_start + si] += g_i;
      }

      // Gradient w.r.t. S_eco
      grad[layout.s_eco_start + si] += g_i;

      // Hessian: beta_eco block
      for (int j = 0; j < p_eco; j++) {
        for (int k = 0; k < p_eco; k++) {
          H_flat[(layout.beta_eco_start + j) * n_x +
                 (layout.beta_eco_start + k)] +=
            h_i * X_eco_flat[i * p_eco + j] * X_eco_flat[i * p_eco + k];
        }
      }

      // Hessian: S_eco diagonal
      H_flat[(layout.s_eco_start + si) * n_x +
             (layout.s_eco_start + si)] += h_i;

      // Hessian: cross beta_eco - S_eco
      for (int j = 0; j < p_eco; j++) {
        int bj = layout.beta_eco_start + j;
        int sj = layout.s_eco_start + si;
        H_flat[bj * n_x + sj] += h_i * X_eco_flat[i * p_eco + j];
        H_flat[sj * n_x + bj] += h_i * X_eco_flat[i * p_eco + j];
      }

      if (has_shared) {
        int ss = layout.s_shared_start + si;

        // Hessian: S_shared diagonal
        H_flat[ss * n_x + ss] += h_i;

        // Cross: beta_eco - S_shared
        for (int j = 0; j < p_eco; j++) {
          int bj = layout.beta_eco_start + j;
          H_flat[bj * n_x + ss] += h_i * X_eco_flat[i * p_eco + j];
          H_flat[ss * n_x + bj] += h_i * X_eco_flat[i * p_eco + j];
        }

        // Cross: S_shared - S_eco
        int se = layout.s_eco_start + si;
        H_flat[ss * n_x + se] += h_i;
        H_flat[se * n_x + ss] += h_i;
      }
    }

    // ==================================================================
    // 2. SAMPLING LGCP LIKELIHOOD contributions
    // ==================================================================

    // 2a. Point log-intensity: sum_i log(lambda(s_i)) = sum_i eta_samp[i]
    for (int i = 0; i < N; i++) {
      int si = obs_spatial_idx[i];

      // eta_samp = Z_samp * beta_samp + beta_share * S_shared[si] + S_samp[si]
      double eta_samp = spjam_linalg::dot_product(
        &Z_samp_flat[i * p_samp], beta_samp, p_samp);
      if (has_shared) eta_samp += beta_share * S_shared[si];
      eta_samp += S_samp[si];

      // Gradient of log lambda(s_i) = eta_samp
      // d/d(beta_samp) = Z_samp[i,:]
      for (int j = 0; j < p_samp; j++) {
        grad[layout.beta_samp_start + j] += Z_samp_flat[i * p_samp + j];
      }

      // d/d(S_samp[si]) = 1
      grad[layout.s_samp_start + si] += 1.0;

      // d/d(S_shared[si]) = beta_share
      if (has_shared) {
        grad[layout.s_shared_start + si] += beta_share;
      }

      // No Hessian from the log-intensity term (it's linear in parameters)
    }

    // 2b. Negative integral: -sum_q w_q * exp(eta_quad[q])
    for (int q = 0; q < N_quad; q++) {
      int sq = quad_spatial_idx[q];

      double eta_quad = spjam_linalg::dot_product(
        &Z_quad_flat[q * p_samp], beta_samp, p_samp);
      if (has_shared) eta_quad += beta_share * S_shared[sq];
      eta_quad += S_samp[sq];

      double lam_q = spjam_linalg::safe_exp(eta_quad);
      double wlam = quad_weights[q] * lam_q;

      // Gradient: -w_q * exp(eta_q) * d(eta_q)/d(param)
      for (int j = 0; j < p_samp; j++) {
        grad[layout.beta_samp_start + j] -=
          wlam * Z_quad_flat[q * p_samp + j];
      }
      grad[layout.s_samp_start + sq] -= wlam;
      if (has_shared) {
        grad[layout.s_shared_start + sq] -= wlam * beta_share;
      }

      // Hessian from -w_q * exp(eta_q):
      // Second derivative = w_q * exp(eta_q) * d(eta)/d(a) * d(eta)/d(b)
      // (positive because negative of second deriv of the negative term)

      // beta_samp block
      for (int j = 0; j < p_samp; j++) {
        for (int k = 0; k < p_samp; k++) {
          H_flat[(layout.beta_samp_start + j) * n_x +
                 (layout.beta_samp_start + k)] +=
            wlam * Z_quad_flat[q * p_samp + j] * Z_quad_flat[q * p_samp + k];
        }
      }

      // S_samp diagonal
      H_flat[(layout.s_samp_start + sq) * n_x +
             (layout.s_samp_start + sq)] += wlam;

      // Cross: beta_samp - S_samp
      for (int j = 0; j < p_samp; j++) {
        int bj = layout.beta_samp_start + j;
        int sj = layout.s_samp_start + sq;
        H_flat[bj * n_x + sj] += wlam * Z_quad_flat[q * p_samp + j];
        H_flat[sj * n_x + bj] += wlam * Z_quad_flat[q * p_samp + j];
      }

      if (has_shared) {
        int ss = layout.s_shared_start + sq;

        // S_shared contribution (with beta_share factor)
        H_flat[ss * n_x + ss] += wlam * beta_share * beta_share;

        // Cross: beta_samp - S_shared
        for (int j = 0; j < p_samp; j++) {
          int bj = layout.beta_samp_start + j;
          H_flat[bj * n_x + ss] +=
            wlam * Z_quad_flat[q * p_samp + j] * beta_share;
          H_flat[ss * n_x + bj] +=
            wlam * Z_quad_flat[q * p_samp + j] * beta_share;
        }

        // Cross: S_shared - S_samp
        int smp = layout.s_samp_start + sq;
        H_flat[ss * n_x + smp] += wlam * beta_share;
        H_flat[smp * n_x + ss] += wlam * beta_share;
      }
    }

    // ==================================================================
    // 3. SPATIAL PRIOR contributions
    // ==================================================================

    if (spatial_type == SPATIAL_GP) {
      // GP prior: gradient and Hessian via dense Cholesky
      if (has_shared && gp_L_shared_flat.size() > 0) {
        gp_gradient_dense(S_shared, gp_L_shared_flat.data(),
                          &grad[layout.s_shared_start], n_sp);
        gp_hessian_dense(gp_L_shared_flat.data(), H_flat.data(),
                         n_sp, n_x, layout.s_shared_start);
      }
      if (gp_L_eco_flat.size() > 0) {
        gp_gradient_dense(S_eco, gp_L_eco_flat.data(),
                          &grad[layout.s_eco_start], n_sp);
        gp_hessian_dense(gp_L_eco_flat.data(), H_flat.data(),
                         n_sp, n_x, layout.s_eco_start);
      }
      if (gp_L_samp_flat.size() > 0) {
        gp_gradient_dense(S_samp, gp_L_samp_flat.data(),
                          &grad[layout.s_samp_start], n_sp);
        gp_hessian_dense(gp_L_samp_flat.data(), H_flat.data(),
                         n_sp, n_x, layout.s_samp_start);
      }

    } else if (spatial_type == SPATIAL_ICAR || spatial_type == SPATIAL_BYM2) {
      // ICAR prior on each spatial field
      if (has_shared) {
        icar_gradient(S_shared, &grad[layout.s_shared_start], n_sp,
                      adj_row_ptr.data(), adj_col_idx.data(),
                      n_neighbors.data(), tau_shared);
        icar_hessian_diag(&H_flat[0], n_sp, n_neighbors.data(), tau_shared);
        // Add off-diagonal ICAR terms
        for (int i = 0; i < n_sp; i++) {
          for (int k = adj_row_ptr[i]; k < adj_row_ptr[i + 1]; k++) {
            int j = adj_col_idx[k];
            int ri = layout.s_shared_start + i;
            int rj = layout.s_shared_start + j;
            H_flat[ri * n_x + rj] += -tau_shared;  // Off-diagonal = -tau
          }
        }
      }

      // Ecological spatial field
      icar_gradient(S_eco, &grad[layout.s_eco_start], n_sp,
                    adj_row_ptr.data(), adj_col_idx.data(),
                    n_neighbors.data(), tau_eco);
      for (int i = 0; i < n_sp; i++) {
        int ri = layout.s_eco_start + i;
        H_flat[ri * n_x + ri] += tau_eco * n_neighbors[i];
        for (int k = adj_row_ptr[i]; k < adj_row_ptr[i + 1]; k++) {
          int j = adj_col_idx[k];
          int rj = layout.s_eco_start + j;
          H_flat[ri * n_x + rj] += -tau_eco;
        }
      }

      // Sampling spatial field
      icar_gradient(S_samp, &grad[layout.s_samp_start], n_sp,
                    adj_row_ptr.data(), adj_col_idx.data(),
                    n_neighbors.data(), tau_samp);
      for (int i = 0; i < n_sp; i++) {
        int ri = layout.s_samp_start + i;
        H_flat[ri * n_x + ri] += tau_samp * n_neighbors[i];
        for (int k = adj_row_ptr[i]; k < adj_row_ptr[i + 1]; k++) {
          int j = adj_col_idx[k];
          int rj = layout.s_samp_start + j;
          H_flat[ri * n_x + rj] += -tau_samp;
        }
      }

    } else if (spatial_type == SPATIAL_SPDE) {
      // SPDE precision
      if (has_shared && spde_Q_shared_vals.size() > 0) {
        spde_gradient(S_shared, &grad[layout.s_shared_start], n_sp,
                      spde_Q_row_ptr.data(), spde_Q_col_idx.data(),
                      spde_Q_shared_vals.data());
        // Add SPDE precision to Hessian
        for (int i = 0; i < n_sp; i++) {
          for (int k = spde_Q_row_ptr[i]; k < spde_Q_row_ptr[i + 1]; k++) {
            int j = spde_Q_col_idx[k];
            int ri = layout.s_shared_start + i;
            int rj = layout.s_shared_start + j;
            H_flat[ri * n_x + rj] += spde_Q_shared_vals[k];
          }
        }
      }

      spde_gradient(S_eco, &grad[layout.s_eco_start], n_sp,
                    spde_Q_row_ptr.data(), spde_Q_col_idx.data(),
                    spde_Q_eco_vals.data());
      for (int i = 0; i < n_sp; i++) {
        for (int k = spde_Q_row_ptr[i]; k < spde_Q_row_ptr[i + 1]; k++) {
          int j = spde_Q_col_idx[k];
          H_flat[(layout.s_eco_start + i) * n_x +
                 (layout.s_eco_start + j)] += spde_Q_eco_vals[k];
        }
      }

      spde_gradient(S_samp, &grad[layout.s_samp_start], n_sp,
                    spde_Q_row_ptr.data(), spde_Q_col_idx.data(),
                    spde_Q_samp_vals.data());
      for (int i = 0; i < n_sp; i++) {
        for (int k = spde_Q_row_ptr[i]; k < spde_Q_row_ptr[i + 1]; k++) {
          int j = spde_Q_col_idx[k];
          H_flat[(layout.s_samp_start + i) * n_x +
                 (layout.s_samp_start + j)] += spde_Q_samp_vals[k];
        }
      }
    }

    // ==================================================================
    // 4. FIXED EFFECTS PRIOR
    // ==================================================================
    for (int j = 0; j < p_eco; j++) {
      int idx = layout.beta_eco_start + j;
      grad[idx] -= tau_beta * x[idx];
      H_flat[idx * n_x + idx] += tau_beta;
    }
    for (int j = 0; j < p_samp; j++) {
      int idx = layout.beta_samp_start + j;
      grad[idx] -= tau_beta * x[idx];
      H_flat[idx * n_x + idx] += tau_beta;
    }

    // ==================================================================
    // 5. NEWTON STEP: solve H * delta = grad
    // ==================================================================
    std::vector<double> L(n_x * n_x);
    double log_det;
    spjam_linalg::cholesky_dense(H_flat, L, log_det, n_x);

    std::vector<double> delta(n_x);
    spjam_linalg::solve_cholesky(L, grad, delta, n_x);

    // Check for NaN/Inf
    double max_delta = 0.0;
    bool has_nan = false;
    for (int j = 0; j < n_x; j++) {
      if (!std::isfinite(delta[j])) {
        has_nan = true;
        break;
      }
      max_delta = std::max(max_delta, std::abs(delta[j]));
    }

    if (has_nan) {
      // Try with damped step
      for (int j = 0; j < n_x; j++) {
        if (std::isfinite(delta[j])) {
          x[j] += 0.1 * delta[j];
        }
      }
      continue;
    }

    // Update
    for (int j = 0; j < n_x; j++) {
      x[j] += delta[j];
    }

    result.n_iter = iter + 1;

    // Check convergence
    if (max_delta < tol) {
      result.converged = true;
      break;
    }
  }

  // ==================================================================
  // 6. COMPUTE FINAL QUANTITIES AT MODE
  // ==================================================================

  // Recompute Hessian at mode (already have it from last iteration if converged)
  // Store mode and Hessian
  result.mode = x;

  // Recompute final Hessian cleanly
  std::fill(H_flat.begin(), H_flat.end(), 0.0);

  const double* beta_eco_final = &x[layout.beta_eco_start];
  const double* beta_samp_final = &x[layout.beta_samp_start];
  const double* S_shared_final = has_shared ? &x[layout.s_shared_start] : nullptr;
  const double* S_eco_final = &x[layout.s_eco_start];
  const double* S_samp_final = &x[layout.s_samp_start];

  // Ecological likelihood at mode
  double log_lik_eco = 0.0;
  for (int i = 0; i < N; i++) {
    int si = obs_spatial_idx[i];
    double eta_eco = spjam_linalg::dot_product(
      &X_eco_flat[i * p_eco], beta_eco_final, p_eco);
    if (has_shared) eta_eco += S_shared_final[si];
    eta_eco += S_eco_final[si];

    log_lik_eco += compute_log_lik(family, y[i], n_trials[i], eta_eco,
                                    phi_overdispersion);
    double h_i = compute_neg_hess(family, y[i], n_trials[i], eta_eco,
                                   phi_overdispersion);

    // Hessian contributions (same structure as above)
    for (int j = 0; j < p_eco; j++) {
      for (int k = 0; k < p_eco; k++) {
        H_flat[(layout.beta_eco_start + j) * n_x +
               (layout.beta_eco_start + k)] +=
          h_i * X_eco_flat[i * p_eco + j] * X_eco_flat[i * p_eco + k];
      }
    }
    H_flat[(layout.s_eco_start + si) * n_x + (layout.s_eco_start + si)] += h_i;
    for (int j = 0; j < p_eco; j++) {
      int bj = layout.beta_eco_start + j;
      int sj = layout.s_eco_start + si;
      H_flat[bj * n_x + sj] += h_i * X_eco_flat[i * p_eco + j];
      H_flat[sj * n_x + bj] += h_i * X_eco_flat[i * p_eco + j];
    }
    if (has_shared) {
      int ss = layout.s_shared_start + si;
      H_flat[ss * n_x + ss] += h_i;
      for (int j = 0; j < p_eco; j++) {
        int bj = layout.beta_eco_start + j;
        H_flat[bj * n_x + ss] += h_i * X_eco_flat[i * p_eco + j];
        H_flat[ss * n_x + bj] += h_i * X_eco_flat[i * p_eco + j];
      }
      int se = layout.s_eco_start + si;
      H_flat[ss * n_x + se] += h_i;
      H_flat[se * n_x + ss] += h_i;
    }
  }

  // Sampling likelihood at mode
  double log_lik_samp = 0.0;
  for (int i = 0; i < N; i++) {
    int si = obs_spatial_idx[i];
    double eta_samp = spjam_linalg::dot_product(
      &Z_samp_flat[i * p_samp], beta_samp_final, p_samp);
    if (has_shared) eta_samp += beta_share * S_shared_final[si];
    eta_samp += S_samp_final[si];
    log_lik_samp += eta_samp;  // log-intensity at observed points
  }

  for (int q = 0; q < N_quad; q++) {
    int sq = quad_spatial_idx[q];
    double eta_quad = spjam_linalg::dot_product(
      &Z_quad_flat[q * p_samp], beta_samp_final, p_samp);
    if (has_shared) eta_quad += beta_share * S_shared_final[sq];
    eta_quad += S_samp_final[sq];

    double lam_q = spjam_linalg::safe_exp(eta_quad);
    double wlam = quad_weights[q] * lam_q;
    log_lik_samp -= wlam;  // Negative integral

    // Hessian from quadrature
    for (int j = 0; j < p_samp; j++) {
      for (int k = 0; k < p_samp; k++) {
        H_flat[(layout.beta_samp_start + j) * n_x +
               (layout.beta_samp_start + k)] +=
          wlam * Z_quad_flat[q * p_samp + j] * Z_quad_flat[q * p_samp + k];
      }
    }
    H_flat[(layout.s_samp_start + sq) * n_x + (layout.s_samp_start + sq)] += wlam;
    for (int j = 0; j < p_samp; j++) {
      int bj = layout.beta_samp_start + j;
      int sj = layout.s_samp_start + sq;
      H_flat[bj * n_x + sj] += wlam * Z_quad_flat[q * p_samp + j];
      H_flat[sj * n_x + bj] += wlam * Z_quad_flat[q * p_samp + j];
    }
    if (has_shared) {
      int ss = layout.s_shared_start + sq;
      H_flat[ss * n_x + ss] += wlam * beta_share * beta_share;
      for (int j = 0; j < p_samp; j++) {
        int bj = layout.beta_samp_start + j;
        H_flat[bj * n_x + ss] += wlam * Z_quad_flat[q * p_samp + j] * beta_share;
        H_flat[ss * n_x + bj] += wlam * Z_quad_flat[q * p_samp + j] * beta_share;
      }
      int smp = layout.s_samp_start + sq;
      H_flat[ss * n_x + smp] += wlam * beta_share;
      H_flat[smp * n_x + ss] += wlam * beta_share;
    }
  }

  // Add spatial prior to final Hessian (same as in Newton loop)
  if (spatial_type == SPATIAL_GP) {
    if (has_shared && gp_L_shared_flat.size() > 0) {
      gp_hessian_dense(gp_L_shared_flat.data(), H_flat.data(), n_sp, n_x,
                       layout.s_shared_start);
    }
    if (gp_L_eco_flat.size() > 0) {
      gp_hessian_dense(gp_L_eco_flat.data(), H_flat.data(), n_sp, n_x,
                       layout.s_eco_start);
    }
    if (gp_L_samp_flat.size() > 0) {
      gp_hessian_dense(gp_L_samp_flat.data(), H_flat.data(), n_sp, n_x,
                       layout.s_samp_start);
    }
  } else if (spatial_type == SPATIAL_ICAR || spatial_type == SPATIAL_BYM2) {
    if (has_shared) {
      for (int i = 0; i < n_sp; i++) {
        int ri = layout.s_shared_start + i;
        H_flat[ri * n_x + ri] += tau_shared * n_neighbors[i];
        for (int k = adj_row_ptr[i]; k < adj_row_ptr[i + 1]; k++) {
          int rj = layout.s_shared_start + adj_col_idx[k];
          H_flat[ri * n_x + rj] += -tau_shared;
        }
      }
    }
    for (int i = 0; i < n_sp; i++) {
      int ri = layout.s_eco_start + i;
      H_flat[ri * n_x + ri] += tau_eco * n_neighbors[i];
      for (int k = adj_row_ptr[i]; k < adj_row_ptr[i + 1]; k++) {
        int rj = layout.s_eco_start + adj_col_idx[k];
        H_flat[ri * n_x + rj] += -tau_eco;
      }
    }
    for (int i = 0; i < n_sp; i++) {
      int ri = layout.s_samp_start + i;
      H_flat[ri * n_x + ri] += tau_samp * n_neighbors[i];
      for (int k = adj_row_ptr[i]; k < adj_row_ptr[i + 1]; k++) {
        int rj = layout.s_samp_start + adj_col_idx[k];
        H_flat[ri * n_x + rj] += -tau_samp;
      }
    }
  } else if (spatial_type == SPATIAL_SPDE) {
    if (has_shared && spde_Q_shared_vals.size() > 0) {
      for (int i = 0; i < n_sp; i++) {
        for (int k = spde_Q_row_ptr[i]; k < spde_Q_row_ptr[i + 1]; k++) {
          int j = spde_Q_col_idx[k];
          H_flat[(layout.s_shared_start + i) * n_x +
                 (layout.s_shared_start + j)] += spde_Q_shared_vals[k];
        }
      }
    }
    for (int i = 0; i < n_sp; i++) {
      for (int k = spde_Q_row_ptr[i]; k < spde_Q_row_ptr[i + 1]; k++) {
        int j = spde_Q_col_idx[k];
        H_flat[(layout.s_eco_start + i) * n_x +
               (layout.s_eco_start + j)] += spde_Q_eco_vals[k];
      }
    }
    for (int i = 0; i < n_sp; i++) {
      for (int k = spde_Q_row_ptr[i]; k < spde_Q_row_ptr[i + 1]; k++) {
        int j = spde_Q_col_idx[k];
        H_flat[(layout.s_samp_start + i) * n_x +
               (layout.s_samp_start + j)] += spde_Q_samp_vals[k];
      }
    }
  }

  // Fixed effects prior in final Hessian
  for (int j = 0; j < p_eco; j++) {
    H_flat[(layout.beta_eco_start + j) * n_x +
           (layout.beta_eco_start + j)] += tau_beta;
  }
  for (int j = 0; j < p_samp; j++) {
    H_flat[(layout.beta_samp_start + j) * n_x +
           (layout.beta_samp_start + j)] += tau_beta;
  }

  // Final Cholesky
  std::vector<double> L_final(n_x * n_x);
  double log_det_final;
  spjam_linalg::cholesky_dense(H_flat, L_final, log_det_final, n_x);

  // Log prior for spatial effects
  double log_prior_spatial = 0.0;
  if (spatial_type == SPATIAL_GP) {
    if (has_shared && gp_L_shared_flat.size() > 0) {
      log_prior_spatial += gp_log_prior_dense(S_shared_final,
                                               gp_L_shared_flat.data(), n_sp);
    }
    log_prior_spatial += gp_log_prior_dense(S_eco_final,
                                             gp_L_eco_flat.data(), n_sp);
    log_prior_spatial += gp_log_prior_dense(S_samp_final,
                                             gp_L_samp_flat.data(), n_sp);
  } else if (spatial_type == SPATIAL_ICAR || spatial_type == SPATIAL_BYM2) {
    if (has_shared) {
      log_prior_spatial -= icar_quadratic(S_shared_final, n_sp,
                                           adj_row_ptr.data(), adj_col_idx.data(),
                                           tau_shared);
    }
    log_prior_spatial -= icar_quadratic(S_eco_final, n_sp,
                                         adj_row_ptr.data(), adj_col_idx.data(),
                                         tau_eco);
    log_prior_spatial -= icar_quadratic(S_samp_final, n_sp,
                                         adj_row_ptr.data(), adj_col_idx.data(),
                                         tau_samp);
    // Add normalizing constants for ICAR (proportional to log(tau))
    // ICAR normalizing constants (proportional to log(tau))
    if (has_shared) {
      log_prior_spatial += 0.5 * (n_sp - 1) * std::log(tau_shared / (2.0 * M_PI));
    }
    log_prior_spatial += 0.5 * (n_sp - 1) * std::log(tau_eco / (2.0 * M_PI));
    log_prior_spatial += 0.5 * (n_sp - 1) * std::log(tau_samp / (2.0 * M_PI));
  } else if (spatial_type == SPATIAL_SPDE) {
    if (has_shared) {
      log_prior_spatial -= spde_quadratic(S_shared_final, n_sp,
                                           spde_Q_row_ptr.data(),
                                           spde_Q_col_idx.data(),
                                           spde_Q_shared_vals.data());
    }
    log_prior_spatial -= spde_quadratic(S_eco_final, n_sp,
                                         spde_Q_row_ptr.data(),
                                         spde_Q_col_idx.data(),
                                         spde_Q_eco_vals.data());
    log_prior_spatial -= spde_quadratic(S_samp_final, n_sp,
                                         spde_Q_row_ptr.data(),
                                         spde_Q_col_idx.data(),
                                         spde_Q_samp_vals.data());
  }

  // Log prior for fixed effects
  double log_prior_beta = 0.0;
  for (int j = 0; j < p_eco; j++) {
    log_prior_beta -= 0.5 * tau_beta * x[layout.beta_eco_start + j] *
                      x[layout.beta_eco_start + j];
  }
  for (int j = 0; j < p_samp; j++) {
    log_prior_beta -= 0.5 * tau_beta * x[layout.beta_samp_start + j] *
                      x[layout.beta_samp_start + j];
  }

  // Laplace marginal: log p(y|theta) ≈ log p(y,x*|theta) - 0.5*log|H| + n_x/2*log(2pi)
  result.log_lik_eco = log_lik_eco;
  result.log_lik_samp = log_lik_samp;
  result.log_det_H = log_det_final;
  result.log_marginal = log_lik_eco + log_lik_samp + log_prior_spatial +
                        log_prior_beta - 0.5 * log_det_final +
                        0.5 * n_x * std::log(2.0 * M_PI);

  result.H_flat = H_flat;
  result.L_flat = L_final;

  return result;
}

} // namespace spjam
