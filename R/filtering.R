# R/filtering.R

#' Filter Data for Stable Region with Diagnostics
#'
#' Filters the dataset to include only observations where the absolute residualized treatment exceeds a specified threshold.
#' Provides diagnostic outputs about the filtering process.
#'
#' @param data A data frame containing residualized data with a residualized treatment variable.
#' @param treatment_resid_var A string specifying the name of the residualized treatment variable (e.g., "D_resid").
#' @param epsilon A numeric threshold for filtering residualized treatment. Default is \code{1e-3}.
#' @return A filtered data frame containing only observations where \code{|treatment_resid_var| > epsilon}.
#' @examples
#' \dontrun{
#' filtered_data <- filter_stable_region(data = mroz_resid, treatment_resid_var = "D_resid", epsilon = 1e-3)
#' }
#' @import dplyr
#' @export
filter_stable_region <- function(data, treatment_resid_var, epsilon = 1e-3) {
  # Check if the treatment_resid_var exists in the data
  if (!treatment_resid_var %in% names(data)) {
    stop(paste("The specified treatment residualized variable", treatment_resid_var, "does not exist in the data."))
  }
  
  initial_n <- nrow(data)
  data_filtered <- data %>%
    filter(abs(.data[[treatment_resid_var]]) > epsilon)
  final_n <- nrow(data_filtered)
  
  cat("Number of observations before stable region filter:", initial_n, "\n")
  cat("Number of observations after stable region filter:", final_n, "\n")
  cat("Fraction of data excluded:", round((initial_n - final_n) / initial_n, 4), "\n\n")
  
  return(data_filtered)
}
