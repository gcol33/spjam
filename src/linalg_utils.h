// linalg_utils.h
// Linear algebra utilities for scatR Laplace engine
// Adapted from numdenom v1.3.0 linalg_fast.h

#ifndef SCATR_LINALG_UTILS_H
#define SCATR_LINALG_UTILS_H

#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>

namespace scatr_linalg {

// Numerically safe exp (clamp to prevent overflow)
inline double safe_exp(double x) {
  if (x > 500.0) return std::exp(500.0);
  if (x < -500.0) return 0.0;
  return std::exp(x);
}

// Dot product with loop unrolling
inline double dot_product(const double* x, const double* y, int n) {
  double sum = 0.0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    sum += x[i] * y[i] + x[i+1] * y[i+1] +
           x[i+2] * y[i+2] + x[i+3] * y[i+3];
  }
  for (; i < n; i++) {
    sum += x[i] * y[i];
  }
  return sum;
}

// Vector sum
inline double vector_sum(const double* x, int n) {
  double sum = 0.0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    sum += x[i] + x[i+1] + x[i+2] + x[i+3];
  }
  for (; i < n; i++) {
    sum += x[i];
  }
  return sum;
}

// axpy: y = a*x + y
inline void axpy(double a, const double* x, double* y, int n) {
  int i = 0;
  for (; i + 3 < n; i += 4) {
    y[i] += a * x[i];
    y[i+1] += a * x[i+1];
    y[i+2] += a * x[i+2];
    y[i+3] += a * x[i+3];
  }
  for (; i < n; i++) {
    y[i] += a * x[i];
  }
}

// Matrix-vector multiply: y = X * beta
// X is N x p stored row-major
inline void matvec(const double* X_flat, const double* beta,
                   double* y, int N, int p) {
  for (int i = 0; i < N; i++) {
    y[i] = dot_product(&X_flat[i * p], beta, p);
  }
}

// Transpose matrix-vector: y = X' * z (accumulate)
// X is N x p row-major, z is N-vector, y is p-vector
inline void matvec_transpose_add(const double* X_flat, const double* z,
                                  double* y, int N, int p) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < p; j++) {
      y[j] += z[i] * X_flat[i * p + j];
    }
  }
}

// Dense Cholesky factorization (lower triangular)
// H is n x n stored as flat row-major, L is output (same layout)
// Returns log determinant. Returns false if not positive definite.
inline bool cholesky_dense(const std::vector<double>& H,
                            std::vector<double>& L,
                            double& log_det, int n) {
  L.assign(n * n, 0.0);
  log_det = 0.0;

  for (int j = 0; j < n; j++) {
    double sum = H[j * n + j];
    for (int k = 0; k < j; k++) {
      sum -= L[j * n + k] * L[j * n + k];
    }
    if (sum <= 0.0) {
      // Add small jitter and retry
      sum = 1e-6;
    }
    L[j * n + j] = std::sqrt(sum);
    log_det += std::log(L[j * n + j]);

    for (int i = j + 1; i < n; i++) {
      sum = H[i * n + j];
      for (int k = 0; k < j; k++) {
        sum -= L[i * n + k] * L[j * n + k];
      }
      L[i * n + j] = sum / L[j * n + j];
    }
  }

  log_det *= 2.0;  // log|H| = 2 * sum(log(diag(L)))
  return true;
}

// Solve L * z = b (lower triangular, dense, row-major)
inline void solve_lower_dense(const std::vector<double>& L,
                               const std::vector<double>& b,
                               std::vector<double>& z, int n) {
  z = b;
  for (int j = 0; j < n; j++) {
    z[j] /= L[j * n + j];
    for (int i = j + 1; i < n; i++) {
      z[i] -= L[i * n + j] * z[j];
    }
  }
}

// Solve L' * x = z (upper triangular from Cholesky, dense, row-major)
inline void solve_upper_dense(const std::vector<double>& L,
                               const std::vector<double>& z,
                               std::vector<double>& x, int n) {
  x = z;
  for (int j = n - 1; j >= 0; j--) {
    for (int i = j + 1; i < n; i++) {
      x[j] -= L[i * n + j] * x[i];
    }
    x[j] /= L[j * n + j];
  }
}

// Solve H * x = b via Cholesky (H = L L')
inline void solve_cholesky(const std::vector<double>& L,
                            const std::vector<double>& b,
                            std::vector<double>& x, int n) {
  std::vector<double> z(n);
  solve_lower_dense(L, b, z, n);
  solve_upper_dense(L, z, x, n);
}

} // namespace scatr_linalg

#endif // SCATR_LINALG_UTILS_H
