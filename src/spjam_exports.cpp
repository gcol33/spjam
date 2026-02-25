// spjam_exports.cpp
// Rcpp export wrappers for joint Laplace functions

#include "laplace_joint.h"
#include "linalg_utils.h"
#include <Rcpp.h>
#include <random>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpp_joint_laplace_fit(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector n_trials,
    Rcpp::NumericMatrix X_eco,
    std::string family,
    double phi_overdispersion,
    Rcpp::NumericMatrix Z_samp,
    Rcpp::NumericMatrix Z_quad,
    Rcpp::NumericVector quad_weights,
    Rcpp::IntegerVector quad_spatial_idx,
    int spatial_type,
    Rcpp::IntegerVector obs_spatial_idx,
    Rcpp::NumericVector gp_L_shared_flat,
    Rcpp::NumericVector gp_L_eco_flat,
    Rcpp::NumericVector gp_L_samp_flat,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors_vec,
    double tau_shared,
    double tau_eco,
    double tau_samp,
    Rcpp::IntegerVector spde_Q_row_ptr,
    Rcpp::IntegerVector spde_Q_col_idx,
    Rcpp::NumericVector spde_Q_shared_vals,
    Rcpp::NumericVector spde_Q_eco_vals,
    Rcpp::NumericVector spde_Q_samp_vals,
    int latent_type,
    double beta_share,
    int p_eco,
    int p_samp,
    int n_spatial,
    double tau_beta,
    int max_iter,
    double tol
) {
  // Build layout
  spjam::JointParamLayout layout;
  layout.p_eco = p_eco;
  layout.p_samp = p_samp;
  layout.n_spatial = n_spatial;
  layout.latent_type = static_cast<spjam::LatentType>(latent_type);

  layout.beta_eco_start = 0;
  layout.beta_samp_start = p_eco;

  if (latent_type == spjam::LATENT_SHARED) {
    layout.s_shared_start = p_eco + p_samp;
    layout.s_eco_start = p_eco + p_samp + n_spatial;
    layout.s_samp_start = p_eco + p_samp + 2 * n_spatial;
    layout.total_params = p_eco + p_samp + 3 * n_spatial;
  } else {
    layout.s_shared_start = -1;
    layout.s_eco_start = p_eco + p_samp;
    layout.s_samp_start = p_eco + p_samp + n_spatial;
    layout.total_params = p_eco + p_samp + 2 * n_spatial;
  }

  // Convert Rcpp types to std::vector
  int N = y.size();
  int N_quad = quad_weights.size();

  std::vector<int> y_vec(y.begin(), y.end());
  std::vector<int> nt_vec(n_trials.begin(), n_trials.end());
  std::vector<double> qw_vec(quad_weights.begin(), quad_weights.end());
  std::vector<int> qs_vec(quad_spatial_idx.begin(), quad_spatial_idx.end());
  std::vector<int> os_vec(obs_spatial_idx.begin(), obs_spatial_idx.end());

  // Design matrices (row-major flat)
  std::vector<double> X_eco_flat(N * p_eco);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < p_eco; j++) {
      X_eco_flat[i * p_eco + j] = X_eco(i, j);
    }
  }

  std::vector<double> Z_samp_flat(N * p_samp);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < p_samp; j++) {
      Z_samp_flat[i * p_samp + j] = Z_samp(i, j);
    }
  }

  std::vector<double> Z_quad_flat(N_quad * p_samp);
  for (int i = 0; i < N_quad; i++) {
    for (int j = 0; j < p_samp; j++) {
      Z_quad_flat[i * p_samp + j] = Z_quad(i, j);
    }
  }

  // GP Cholesky (may be empty)
  std::vector<double> gp_Ls(gp_L_shared_flat.begin(), gp_L_shared_flat.end());
  std::vector<double> gp_Le(gp_L_eco_flat.begin(), gp_L_eco_flat.end());
  std::vector<double> gp_Lp(gp_L_samp_flat.begin(), gp_L_samp_flat.end());

  // ICAR adjacency
  std::vector<int> arp(adj_row_ptr.begin(), adj_row_ptr.end());
  std::vector<int> aci(adj_col_idx.begin(), adj_col_idx.end());
  std::vector<int> nn(n_neighbors_vec.begin(), n_neighbors_vec.end());

  // SPDE
  std::vector<int> sqrp(spde_Q_row_ptr.begin(), spde_Q_row_ptr.end());
  std::vector<int> sqci(spde_Q_col_idx.begin(), spde_Q_col_idx.end());
  std::vector<double> sqsv(spde_Q_shared_vals.begin(), spde_Q_shared_vals.end());
  std::vector<double> sqev(spde_Q_eco_vals.begin(), spde_Q_eco_vals.end());
  std::vector<double> sqpv(spde_Q_samp_vals.begin(), spde_Q_samp_vals.end());

  // Run joint Laplace
  spjam::JointLaplaceResult res = spjam::joint_laplace_mode(
    y_vec, nt_vec, X_eco_flat, family, phi_overdispersion,
    Z_samp_flat, Z_quad_flat, qw_vec, qs_vec,
    spatial_type, os_vec,
    gp_Ls, gp_Le, gp_Lp,
    arp, aci, nn, tau_shared, tau_eco, tau_samp,
    sqrp, sqci, sqsv, sqev, sqpv,
    latent_type, beta_share,
    layout, tau_beta, max_iter, tol
  );

  // Return as R list
  return Rcpp::List::create(
    Named("mode") = NumericVector(res.mode.begin(), res.mode.end()),
    Named("H_flat") = NumericVector(res.H_flat.begin(), res.H_flat.end()),
    Named("L_flat") = NumericVector(res.L_flat.begin(), res.L_flat.end()),
    Named("log_det_H") = res.log_det_H,
    Named("log_marginal") = res.log_marginal,
    Named("log_lik_eco") = res.log_lik_eco,
    Named("log_lik_samp") = res.log_lik_samp,
    Named("n_iter") = res.n_iter,
    Named("converged") = res.converged,
    Named("n_params") = layout.total_params
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_joint_laplace_sample(
    Rcpp::NumericVector mode,
    Rcpp::NumericVector L_flat,
    int n_params,
    int n_samples,
    int seed
) {
  int n = n_params;
  NumericMatrix samples(n_samples, n);

  std::mt19937 rng(seed);
  std::normal_distribution<double> std_normal(0.0, 1.0);

  for (int s = 0; s < n_samples; s++) {
    // Generate standard normal vector z
    std::vector<double> z(n);
    for (int j = 0; j < n; j++) {
      z[j] = std_normal(rng);
    }

    // Solve L' * w = z to get w ~ N(0, H^{-1})
    // (L is Cholesky of H, so L'^{-1} z ~ N(0, H^{-1}))
    std::vector<double> w(n);
    for (int j = n - 1; j >= 0; j--) {
      double sum = z[j];
      for (int k = j + 1; k < n; k++) {
        sum -= L_flat[k * n + j] * w[k];
      }
      w[j] = sum / L_flat[j * n + j];
    }

    // Sample = mode + w
    for (int j = 0; j < n; j++) {
      samples(s, j) = mode[j] + w[j];
    }
  }

  return samples;
}
