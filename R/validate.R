#' Validate inputs to spjam()
#' @keywords internal
validate_spjam_inputs <- function(formula, data, locations, family,
                                   spatial, shared, sampling_spec) {
  # Formula

  if (!inherits(formula, "spjam_formula")) {
    stop("'formula' must be a spjam_formula object", call. = FALSE)
  }

  # Data
  if (!is.data.frame(data)) {
    stop(sprintf("'data' must be a data.frame, got %s", class(data)[1]),
         call. = FALSE)
  }
  if (nrow(data) == 0) {
    stop("'data' must not be empty", call. = FALSE)
  }

  # Check formula variables exist in data
  eco_vars <- all.vars(formula$ecological)
  samp_vars <- all.vars(formula$sampling)
  missing_eco <- setdiff(eco_vars, c(names(data), ""))
  missing_samp <- setdiff(samp_vars, c(names(data), ""))

  # locations may be referenced via formula
  if (inherits(locations, "formula")) {
    loc_vars <- all.vars(locations)
    missing_eco <- setdiff(missing_eco, loc_vars)
    missing_samp <- setdiff(missing_samp, loc_vars)
    missing_loc <- setdiff(loc_vars, names(data))
    if (length(missing_loc) > 0) {
      stop(sprintf("location variables not found in data: %s",
                   paste(missing_loc, collapse = ", ")), call. = FALSE)
    }
  }

  if (length(missing_eco) > 0) {
    stop(sprintf("ecological formula variables not found in data: %s",
                 paste(missing_eco, collapse = ", ")), call. = FALSE)
  }
  if (length(missing_samp) > 0) {
    stop(sprintf("sampling formula variables not found in data: %s",
                 paste(missing_samp, collapse = ", ")), call. = FALSE)
  }

  # Family
  if (!inherits(family, "spjam_family")) {
    stop("'family' must be a spjam_family object", call. = FALSE)
  }

  # Spatial
  if (!inherits(spatial, "spjam_spatial")) {
    stop("'spatial' must be a spjam_spatial object (spatial_gp, spatial_car, etc.)",
         call. = FALSE)
  }

  # Shared/latent

  if (!inherits(shared, "spjam_latent")) {
    stop("'shared' must be a spjam_latent object (latent_shared, latent_independent, etc.)",
         call. = FALSE)
  }

  # Locations
  if (inherits(locations, "formula")) {
    # OK
  } else if (is.matrix(locations)) {
    if (nrow(locations) != nrow(data)) {
      stop(sprintf("locations has %d rows but data has %d rows",
                   nrow(locations), nrow(data)), call. = FALSE)
    }
    if (ncol(locations) != 2) {
      stop("locations matrix must have exactly 2 columns (x, y)", call. = FALSE)
    }
  } else {
    stop("'locations' must be a formula (~ x + y) or a 2-column matrix",
         call. = FALSE)
  }

  # CAR/BYM2 require adjacency
  if (inherits(spatial, "spjam_spatial_car") ||
      inherits(spatial, "spjam_spatial_bym2")) {
    adj <- spatial$adjacency
    if (is.null(adj)) {
      stop("spatial_car() and spatial_bym2() require an adjacency matrix",
           call. = FALSE)
    }
  }

  invisible(TRUE)
}
