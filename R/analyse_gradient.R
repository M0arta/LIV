# R/analyse_gradients.R

#' Analyze Gradients for Selected Bandwidths
#'
#' Fits non-parametric regression models using selected bandwidths and extracts the gradients.
#'
#' @param bw_list A named list of bandwidth objects as returned by \code{select_bandwidths}.
#' @param data A data frame containing the residualized treatment and instrument variables.
#' @param treatment_resid_var A string specifying the name of the residualized treatment variable (e.g., "D_resid").
#' @param instrument_resid_var A string specifying the name of the residualized instrument variable (e.g., "Z_resid").
#' @return A named list of gradient vectors corresponding to each bandwidth method.
#' @examples
#' \dontrun{
#' gradients <- analyse_gradients(
#'   bw_list = bw_list,
#'   data = all_tot,
#'   treatment_resid_var = "D_resid",
#'   instrument_resid_var = "Z_resid"
#' )
#' }
#' @import np
#' @export
analyse_gradients <- function(bw_list, data, treatment_resid_var, instrument_resid_var) {
  # Validate that bandwidth list is not empty
  if (length(bw_list) == 0) {
    stop("Bandwidth list is empty. Please provide valid bandwidth objects.")
  }
  
  # Validate that required variables exist in the data
  required_vars <- c(treatment_resid_var, instrument_resid_var)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("The following required variables are missing in the data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Initialize list to store gradient vectors
  grad_list <- list()
  
  # Iterate over each bandwidth and fit the model
  for (name in names(bw_list)) {
    bw <- bw_list[[name]]
    model <- npreg(
      bws = bw,
      tydat = data[[treatment_resid_var]],
      txdat = data[[instrument_resid_var]],
      gradients = TRUE
    )
    grad <- model$grad[, 1]  # Extract the gradient
    grad_list[[name]] <- grad
  }
  
  return(grad_list)
}
