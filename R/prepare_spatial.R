#' Prepare spatial structure for C++ backend
#'
#' S3 generic that dispatches on spatial specification class to build
#' C++-ready data structures.
#'
#' @keywords internal
prepare_spatial_structure <- function(spatial, model_data) {
  UseMethod("prepare_spatial_structure")
}

#' @keywords internal
prepare_spatial_structure.scatr_spatial_gp <- function(spatial, model_data) {
  locs <- model_data$unique_locations
  n <- nrow(locs)
  nu <- spatial$nu

  # For moderate n, use dense GP covariance
  # Build distance matrix
  dist_mat <- as.matrix(dist(locs))

  list(
    type = "gp",
    spatial_type = 0L,  # SPATIAL_GP
    n_spatial = n,
    coords = locs,
    dist_mat = dist_mat,
    nu = nu,
    # Adjacency fields (empty for GP)
    adj_row_ptr = integer(0),
    adj_col_idx = integer(0),
    n_neighbors = integer(0),
    # SPDE fields (empty for GP)
    spde_Q_row_ptr = integer(0),
    spde_Q_col_idx = integer(0)
  )
}

#' @keywords internal
prepare_spatial_structure.scatr_spatial_car <- function(spatial, model_data) {
  adj <- spatial$adjacency

  # Convert adjacency to CSR format
  csr <- adjacency_to_csr(adj)
  n <- nrow(adj)

  list(
    type = "car",
    spatial_type = 1L,  # SPATIAL_ICAR
    n_spatial = n,
    adj_row_ptr = csr$row_ptr,
    adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors,
    car_type = spatial$car_type,
    # GP fields (empty)
    coords = matrix(0, 0, 2),
    dist_mat = matrix(0, 0, 0),
    nu = 1.5,
    # SPDE fields (empty)
    spde_Q_row_ptr = integer(0),
    spde_Q_col_idx = integer(0)
  )
}

#' @keywords internal
prepare_spatial_structure.scatr_spatial_bym2 <- function(spatial, model_data) {
  adj <- spatial$adjacency

  csr <- adjacency_to_csr(adj)
  n <- nrow(adj)

  # Compute BYM2 scaling factor (geometric mean of marginal variances)
  scale_factor <- compute_bym2_scale(adj)

  list(
    type = "bym2",
    spatial_type = 2L,  # SPATIAL_BYM2
    n_spatial = n,
    adj_row_ptr = csr$row_ptr,
    adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors,
    scale_factor = scale_factor,
    # GP fields (empty)
    coords = matrix(0, 0, 2),
    dist_mat = matrix(0, 0, 0),
    nu = 1.5,
    # SPDE fields (empty)
    spde_Q_row_ptr = integer(0),
    spde_Q_col_idx = integer(0)
  )
}

#' @keywords internal
prepare_spatial_structure.scatr_spatial_spde <- function(spatial, model_data) {
  locs <- model_data$unique_locations
  n <- nrow(locs)

  # Build SPDE mesh and matrices
  if (!is.null(spatial$mesh)) {
    mesh <- spatial$mesh
  } else {
    if (!requireNamespace("fmesher", quietly = TRUE)) {
      stop("fmesher package required for spatial_spde(). Install with: ",
           "install.packages('fmesher')", call. = FALSE)
    }
    mesh <- fmesher::fm_mesh_2d(loc = locs, max.edge = c(0.1, 0.3) * max(dist(locs)))
  }

  # Build FEM matrices from mesh
  fem <- fmesher::fm_fem(mesh)
  C_mat <- fem$c0  # Mass matrix
  G_mat <- fem$g1  # Stiffness matrix

  # Projection matrix (maps mesh vertices to observation locations)
  A_mat <- fmesher::fm_basis(mesh, loc = locs)

  # SPDE precision structure: Q = tau^2 * (kappa^4 * C + 2*kappa^2 * G1 + G2)
  # We pass the sparse matrices and let the hyperparameter loop build Q
  # For now, store the components

  # Convert to CSR for the C++ side
  # We'll build the actual Q in the R backend given hyperparameters
  n_mesh <- mesh$n

  list(
    type = "spde",
    spatial_type = 3L,  # SPATIAL_SPDE
    n_spatial = n_mesh,
    mesh = mesh,
    C_mat = C_mat,
    G_mat = G_mat,
    A_mat = A_mat,
    # Adjacency fields (empty)
    adj_row_ptr = integer(0),
    adj_col_idx = integer(0),
    n_neighbors = integer(0),
    # GP fields (empty)
    coords = locs,
    dist_mat = matrix(0, 0, 0),
    nu = 1.5,
    # SPDE CSR will be built per hyperparameter evaluation
    spde_Q_row_ptr = integer(0),
    spde_Q_col_idx = integer(0)
  )
}

#' Convert adjacency matrix to CSR format
#' @keywords internal
adjacency_to_csr <- function(adj) {
  if (inherits(adj, "Matrix")) {
    adj <- as.matrix(adj)
  }
  n <- nrow(adj)
  row_ptr <- integer(n + 1)
  col_idx <- integer(0)
  n_neighbors <- integer(n)

  idx <- 0L
  for (i in seq_len(n)) {
    row_ptr[i] <- idx
    neighbors <- which(adj[i, ] != 0 & seq_len(n) != i)
    n_neighbors[i] <- length(neighbors)
    col_idx <- c(col_idx, neighbors - 1L)  # 0-based
    idx <- idx + length(neighbors)
  }
  row_ptr[n + 1] <- idx

  list(
    row_ptr = row_ptr,
    col_idx = col_idx,
    n_neighbors = n_neighbors
  )
}

#' Compute BYM2 scaling factor
#' @keywords internal
compute_bym2_scale <- function(adj) {
  n <- nrow(adj)
  # Build ICAR precision matrix Q = D - A
  D <- diag(rowSums(adj != 0))
  A <- adj
  diag(A) <- 0
  Q <- D - A

  # Add small jitter for invertibility (ICAR is rank-deficient)
  Q_reg <- Q + diag(1e-6, n)

  # Geometric mean of marginal variances = exp(mean(log(diag(solve(Q_reg)))))
  marginal_var <- diag(solve(Q_reg))
  scale_factor <- exp(mean(log(pmax(marginal_var, 1e-10))))

  scale_factor
}

#' Build GP covariance Cholesky for given hyperparameters
#' @keywords internal
build_gp_cholesky <- function(dist_mat, range, sigma, nu) {
  n <- nrow(dist_mat)
  C <- matrix(0, n, n)

  for (i in seq_len(n)) {
    C[i, i] <- sigma^2 + 1e-6  # Small jitter for numerical stability
    for (j in seq_len(i - 1)) {
      d <- dist_mat[i, j]
      if (nu <= 0.75) {
        # Exponential
        C[i, j] <- sigma^2 * exp(-d / range)
      } else if (nu <= 2.0) {
        # Matern 3/2
        r <- sqrt(3) * d / range
        C[i, j] <- sigma^2 * (1 + r) * exp(-r)
      } else {
        # Matern 5/2
        r <- sqrt(5) * d / range
        C[i, j] <- sigma^2 * (1 + r + r^2 / 3) * exp(-r)
      }
      C[j, i] <- C[i, j]
    }
  }

  L <- tryCatch(
    t(chol(C)),
    error = function(e) {
      # Add more jitter
      diag(C) <- diag(C) + 1e-4
      t(chol(C))
    }
  )

  as.numeric(t(L))  # Row-major flat for C++
}

#' Build SPDE precision matrix for given hyperparameters
#' @keywords internal
build_spde_precision <- function(spde_info, range, sigma) {
  # kappa = sqrt(8 * nu) / range, for nu = 1 (alpha = 2)
  kappa <- sqrt(8) / range
  tau <- 1 / (sigma * sqrt(4 * pi * kappa^2))

  # Q = tau^2 * (kappa^4 * C + 2 * kappa^2 * G)
  Q <- tau^2 * (kappa^4 * spde_info$C_mat + 2 * kappa^2 * spde_info$G_mat)

  # Convert to CSR
  Q_sparse <- Matrix::drop0(Q, tol = 1e-15)
  Q_csc <- as(Q_sparse, "dgCMatrix")

  # dgCMatrix stores in CSC (compressed sparse column)
  # Convert to CSR by transposing
  Q_csr <- Matrix::t(Q_csc)

  list(
    row_ptr = Q_csr@p,
    col_idx = Q_csr@i,
    vals = Q_csr@x
  )
}
