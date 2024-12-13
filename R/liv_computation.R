# R/liv_computation.R

#' Compute LIV Estimates over a Grid
#'
#' Computes Local Instrumental Variable (LIV) estimates at specified instrument grid points using non-parametric regression.
#'
#' @param data A data frame containing residualized data with outcome, treatment, and instrument residuals.
#' @param instrument_grid A numeric vector specifying grid points for the residualized instrument.
#' @param outcome_resid_var A string specifying the name of the residualized outcome variable (e.g., "Y_resid").
#' @param treatment_resid_var A string specifying the name of the residualized treatment variable (e.g., "D_resid").
#' @param instrument_resid_var A string specifying the name of the residualized instrument variable (e.g., "Z_resid").
#' @param bw Bandwidth for non-parametric regression. Default is \code{0.05}.
#' @param kernel Kernel type for non-parametric regression. Options include \code{"epanechnikov"}, \code{"gaussian"}, etc. Default is \code{"epanechnikov"}.
#' @return A numeric vector of LIV estimates corresponding to each point in \code{instrument_grid}.
#' @examples
#' \dontrun{
#' LIV_estimates <- compute_LIV_curve(
#'   data = mroz_filtered,
#'   instrument_grid = instrument_grid,
#'   outcome_resid_var = "Y_resid",
#'   treatment_resid_var = "D_resid",
#'   instrument_resid_var = "Z_resid",
#'   bw = 0.05,
#'   kernel = "epanechnikov"
#' )
#' }
#' @import purrr
#' @import np
#' @export
compute_LIV_curve <- function(data, instrument_grid, outcome_resid_var, treatment_resid_var, instrument_resid_var, bw = 0.05, kernel = "epanechnikov") {
  # Check if required variables exist in the data
  required_vars <- c(outcome_resid_var, treatment_resid_var, instrument_resid_var)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("The following required variables are missing in the data:", paste(missing_vars, collapse = ", ")))
  }
  
  LIV_estimates <- map_dbl(instrument_grid, function(z_val) {
    # Compute derivative of outcome with respect to instrument
    deriv_Y <- npreg(
      tydat = data[[outcome_resid_var]],
      txdat = data[[instrument_resid_var]],
      exdat = data.frame(Z_resid = z_val),
      gradients = TRUE,
      bws = bw,
      ckertype = kernel
    )
    
    # Compute derivative of treatment with respect to instrument
    deriv_D <- npreg(
      tydat = data[[treatment_resid_var]],
      txdat = data[[instrument_resid_var]],
      exdat = data.frame(Z_resid = z_val),
      gradients = TRUE,
      bws = bw,
      ckertype = kernel
    )
    
    dE_Y_dZ <- deriv_Y$grad
    dE_D_dZ <- deriv_D$grad
    epsilon <- 1e-8
    
    if (abs(dE_D_dZ) < epsilon || is.na(dE_D_dZ)) {
      return(NA_real_)
    } else {
      return(dE_Y_dZ / dE_D_dZ)
    }
  })
  
  return(LIV_estimates)
}
