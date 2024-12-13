# R/plot_gradient_distributions.R

#' Plot Gradient Distributions for Different Bandwidths
#'
#' Generates histograms and scatter plots to compare gradient distributions across different bandwidth selection methods.
#'
#' @param grad_list A named list of gradient vectors as returned by \code{analyze_gradients}.
#' @param instrument_grid A numeric vector specifying the grid points for the instrument.
#' @return A list containing ggplot objects for histograms and scatter plots.
#' @examples
#' \dontrun{
#' plots <- plot_gradient_distributions(
#'   grad_list = gradients,
#'   instrument_grid = instrument_grid
#' )
#' print(plots$hist_lscv)
#' print(plots$hist_aic)
#' print(plots$hist_default)
#' 
#' print(plots$scatter_lscv)
#' print(plots$scatter_aic)
#' print(plots$scatter_default)
#' }
#' @import ggplot2
#' @import dplyr
#' @export
plot_gradient_distributions <- function(grad_list, instrument_grid) {
  # Validate that grad_list is not empty
  if (length(grad_list) == 0) {
    stop("Gradient list is empty. Please provide valid gradient vectors.")
  }
  
  # Initialize lists to store plots
  hist_plots <- list()
  scatter_plots <- list()
  
  # Iterate over each gradient vector and create plots
  for (name in names(grad_list)) {
    grad <- grad_list[[name]]
    
    # Create a data frame for plotting
    plot_df <- data.frame(
      Z_resid = instrument_grid,
      Gradient = grad
    )
    
    # Histogram of Gradients
    hist_p <- ggplot(plot_df, aes(x = Gradient)) +
      geom_histogram(binwidth = 0.001, fill = "lightblue", color = "black") +
      labs(
        title = paste("Gradient Distribution:", name, "Bandwidth"),
        x = "Gradient dE[D|Z]/dZ",
        y = "Frequency"
      ) +
      theme_minimal()
    
    hist_plots[[paste0("hist_", name)]] <- hist_p
    
    # Scatter Plot of Gradients vs Z_resid with Lowess Smoother
    scatter_p <- ggplot(plot_df, aes(x = Z_resid, y = Gradient)) +
      geom_point(alpha = 0.5, color = "darkgreen") +
      geom_smooth(method = "loess", color = "red") +
      labs(
        title = paste("Gradients vs Z_resid:", name, "Bandwidth"),
        x = "Residualized Instrument (Z_resid)",
        y = "Gradient dE[D|Z]/dZ"
      ) +
      theme_minimal()
    
    scatter_plots[[paste0("scatter_", name)]] <- scatter_p
  }
  
  return(list(
    histograms = hist_plots,
    scatter_plots = scatter_plots
  ))
}
