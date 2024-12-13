# R/liv_computation.R

#' Bootstrap LIV Estimates
#'
#' Performs bootstrapping to obtain confidence intervals for LIV estimates.
#'
#' @param data A data frame containing residualized data with outcome, treatment, and instrument residuals.
#' @param instrument_grid A numeric vector specifying grid points for the residualized instrument.
#' @param outcome_resid_var A string specifying the name of the residualized outcome variable (e.g., "Y_resid").
#' @param treatment_resid_var A string specifying the name of the residualized treatment variable (e.g., "D_resid").
#' @param instrument_resid_var A string specifying the name of the residualized instrument variable (e.g., "Z_resid").
#' @param bw Bandwidth for non-parametric regression. Default is \code{0.05}.
#' @param kernel Kernel type for non-parametric regression. Options include \code{"epanechnikov"}, \code{"gaussian"}, etc. Default is \code{"epanechnikov"}.
#' @return A numeric vector of bootstrapped LIV estimates corresponding to each point in \code{instrument_grid}.
#' @examples
#' \dontrun{
#' boot_LIV_estimates <- compute_LIV_boot(
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
compute_LIV_boot <- function(data, instrument_grid, outcome_resid_var, treatment_resid_var, instrument_resid_var, bw = 0.05, kernel = "epanechnikov") {
  # Resample the data based on indices
  # Note: The actual resampling is handled by the boot package, so this function assumes data is already resampled.
  
  # Check if required variables exist in the data
  required_vars <- c(outcome_resid_var, treatment_resid_var, instrument_resid_var)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("The following required variables are missing in the data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Compute LIV estimates
  LIV_boot_estimates <- compute_LIV_curve(
    data = data,
    instrument_grid = instrument_grid,
    outcome_resid_var = outcome_resid_var,
    treatment_resid_var = treatment_resid_var,
    instrument_resid_var = instrument_resid_var,
    bw = bw,
    kernel = kernel
  )
  
  return(LIV_boot_estimates)
}
