// laplace_joint.h
// Joint Laplace approximation for preferential sampling models
// Core data structures and function declarations

#ifndef SPJAM_LAPLACE_JOINT_H
#define SPJAM_LAPLACE_JOINT_H

#include <Rcpp.h>
#include <vector>
#include <string>

namespace spjam {

// Spatial type enumeration
enum SpatialType {
  SPATIAL_GP = 0,
  SPATIAL_ICAR = 1,
  SPATIAL_BYM2 = 2,
  SPATIAL_SPDE = 3
};

// Latent structure type
enum LatentType {
  LATENT_SHARED = 0,
  LATENT_INDEPENDENT = 1
};

// Layout of the latent parameter vector x
// For SHARED: x = [beta_eco | beta_samp | S_shared | S_eco_only | S_samp_only]
// For INDEPENDENT: x = [beta_eco | beta_samp | S_eco | S_samp]
struct JointParamLayout {
  int p_eco;           // Number of ecological fixed effects (incl. intercept)
  int p_samp;          // Number of sampling fixed effects (incl. intercept)
  int n_spatial;       // Number of spatial units

  int beta_eco_start;
  int beta_samp_start;
  int s_shared_start;   // -1 if independent
  int s_eco_start;
  int s_samp_start;
  int total_params;

  LatentType latent_type;
};

// Result of joint Laplace mode finding
struct JointLaplaceResult {
  std::vector<double> mode;     // Mode of latent field
  std::vector<double> H_flat;   // Dense Hessian at mode (row-major)
  std::vector<double> L_flat;   // Cholesky of Hessian
  double log_det_H;             // Log determinant of Hessian
  double log_marginal;          // Laplace marginal likelihood
  double log_lik_eco;           // Ecological log-likelihood at mode
  double log_lik_samp;          // Sampling log-likelihood at mode
  int n_iter;
  bool converged;
};

// Joint Laplace mode finding for the coupled preferential sampling model
//
// Given hyperparameters (range, sigma per field, beta_share, phi_overdispersion),
// find the mode of the latent field x = [betas | spatial effects] by Newton-Raphson.
//
// Returns the mode, Hessian, and Laplace marginal likelihood.
JointLaplaceResult joint_laplace_mode(
    // Ecological data
    const std::vector<int>& y,
    const std::vector<int>& n_trials,
    const std::vector<double>& X_eco_flat,  // [N x p_eco] row-major
    const std::string& family,
    double phi_overdispersion,

    // Sampling data
    const std::vector<double>& Z_samp_flat, // [N x p_samp] row-major
    const std::vector<double>& Z_quad_flat, // [N_quad x p_samp] row-major
    const std::vector<double>& quad_weights,
    const std::vector<int>& quad_spatial_idx,

    // Spatial structure
    int spatial_type,
    const std::vector<int>& obs_spatial_idx,

    // GP data (if spatial_type == GP)
    const std::vector<double>& gp_L_shared_flat,
    const std::vector<double>& gp_L_eco_flat,
    const std::vector<double>& gp_L_samp_flat,

    // ICAR/BYM2 data (if spatial_type == ICAR or BYM2)
    const std::vector<int>& adj_row_ptr,
    const std::vector<int>& adj_col_idx,
    const std::vector<int>& n_neighbors,
    double tau_shared,
    double tau_eco,
    double tau_samp,

    // SPDE data (if spatial_type == SPDE)
    const std::vector<int>& spde_Q_row_ptr,
    const std::vector<int>& spde_Q_col_idx,
    const std::vector<double>& spde_Q_shared_vals,
    const std::vector<double>& spde_Q_eco_vals,
    const std::vector<double>& spde_Q_samp_vals,

    // Shared structure
    int latent_type,
    double beta_share,

    // Layout
    const JointParamLayout& layout,

    // Prior precision for fixed effects
    double tau_beta,

    // Control
    int max_iter,
    double tol
);

} // namespace spjam

#endif // SPJAM_LAPLACE_JOINT_H
