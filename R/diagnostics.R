# R/diagnostics.R

#' Plot Correlation Between Residualized D and Z
#'
#' Generates a scatter plot with a fitted regression line to assess the correlation between residualized treatment and instrument variables.
#'
#' @param data A data frame containing residualized treatment and instrument variables.
#' @param treatment_resid_var A string specifying the name of the residualized treatment variable.
#' @param instrument_resid_var A string specifying the name of the residualized instrument variable.
#' @return A ggplot object visualizing the correlation.
#' @examples
#' \dontrun{
#' diagnostic_plot <- plot_correlation_D_Z(
#'   data = mroz_filtered,
#'   treatment_resid_var = "D_resid",
#'   instrument_resid_var = "Z_resid"
#' )
#' print(diagnostic_plot)
#' }
#' @import ggplot2
#' @export
plot_correlation_D_Z <- function(data, treatment_resid_var, instrument_resid_var) {
  # Check if variables exist
  required_vars <- c(treatment_resid_var, instrument_resid_var)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("The following required variables are missing in the data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Compute correlation
  correlation <- cor(data[[treatment_resid_var]], data[[instrument_resid_var]], use = "complete.obs")
  cat("Correlation between", treatment_resid_var, "and", instrument_resid_var, ":", correlation, "\n")
  
  # Generate plot
  p <- ggplot(data, aes_string(x = instrument_resid_var, y = treatment_resid_var)) +
    geom_point(alpha = 0.5, color = "darkgreen") +
    geom_smooth(method = "lm", color = "red") +
    labs(
      title = paste("Correlation Between", treatment_resid_var, "and", instrument_resid_var),
      x = paste("Residualized Instrument (", instrument_resid_var, ")", sep = ""),
      y = paste("Residualized Treatment (", treatment_resid_var, ")", sep = "")
    ) +
    theme_minimal()
  
  return(p)
}

