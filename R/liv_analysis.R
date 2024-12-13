# R/liv_analysis.R

#' Run LIV Analysis for a Given Outcome with Diagnostics
#'
#' Performs the complete Local Instrumental Variable (LIV) analysis pipeline for a specified outcome variable, including optional bandwidth selection and gradient diagnostics.
#'
#' @param outcome_var A string specifying the name of the outcome variable (e.g., "Y").
#' @param data A data frame containing the dataset with necessary variables.
#' @param treatment_var A string specifying the name of the treatment variable (e.g., "D").
#' @param instrument_var A string specifying the name of the instrument variable (e.g., "Z").
#' @param covariates A string vector specifying the names of covariate variables.
#' @param fixed_effects A string vector specifying the names of fixed effect variables.
#' @param R Number of bootstrap replications. Default is \code{1000}.
#' @param seed Random seed for bootstrapping. Default is \code{123}.
#' @param bw Bandwidth for non-parametric regression. Default is \code{0.05}.
#' @param kernel Kernel type for non-parametric regression. Options include \code{"epanechnikov"}, \code{"gaussian"}, etc. Default is \code{"epanechnikov"}.
#' @param epsilon Threshold for filtering residualized treatment. Default is \code{1e-3}.
#' @param perform_diagnostics Logical indicating whether to perform bandwidth selection and gradient diagnostics. Default is \code{FALSE}.
#' @param diagnostic_methods A character vector specifying bandwidth selection methods for diagnostics. Default is \code{c("default", "cv.aic", "cv.ls")}.
#' @param workers Number of parallel workers for bootstrapping. Default is the number of available cores.
#' @return A list containing:
#' \item{outcome}{Name of the outcome variable.}
#' \item{model_Y}{The fitted fixed effects model for the outcome variable.}
#' \item{LIV_estimates}{Vector of LIV estimates over the instrument grid.}
#' \item{instrument_grid}{Numeric vector of instrument grid points.}
#' \item{ci_lower}{Lower bounds of the 95\% confidence intervals for LIV estimates.}
#' \item{ci_upper}{Upper bounds of the 95\% confidence intervals for LIV estimates.}
#' \item{plot}{A ggplot object visualizing the LIV curve with confidence intervals.}
#' \item{diagnostics}{A list containing diagnostic plots if \code{perform_diagnostics = TRUE}.}
#' @examples
#' \dontrun{
#' results <- run_LIV_for_outcome(
#'   outcome_var = "hours",
#'   data = mroz_data,
#'   treatment_var = "wage",
#'   instrument_var = "kids",
#'   covariates = c("age", "educ", "married"),
#'   fixed_effects = c("region"),
#'   R = 1000,
#'   seed = 123,
#'   bw = 0.05,
#'   kernel = "epanechnikov",
#'   epsilon = 1e-3,
#'   perform_diagnostics = TRUE,
#'   diagnostic_methods = c("default", "cv.aic", "cv.ls"),
#'   workers = 4
#' )
#' 
#' # Plot the LIV curve
#' print(results$plot)
#' 
#' # View diagnostic histograms
#' print(results$diagnostics$histograms$hist_default)
#' print(results$diagnostics$histograms$hist_cv.aic)
#' print(results$diagnostics$histograms$hist_cv.ls)
#' 
#' # View diagnostic scatter plots
#' print(results$diagnostics$scatter_plots$scatter_default)
#' print(results$diagnostics$scatter_plots$scatter_cv.aic)
#' print(results$diagnostics$scatter_plots$scatter_cv.ls)
#' }
#' @import dplyr
#' @import purrr
#' @import fixest
#' @import np
#' @import boot
#' @import ggplot2
#' @import future
#' @import furrr
#' @export
run_LIV_for_outcome <- function(outcome_var,
                                data,
                                treatment_var,
                                instrument_var,
                                covariates,
                                fixed_effects,
                                R = 1000,
                                seed = 123,
                                bw = 0.05,
                                kernel = "epanechnikov",
                                epsilon = 1e-3,
                                perform_diagnostics = FALSE,
                                diagnostic_methods = c("default", "cv.aic", "cv.ls"),
                                workers = NULL) {
  message("Processing outcome: ", outcome_var)
  
  # Check if specified variables exist in the data
  required_vars <- c(outcome_var, treatment_var, instrument_var, covariates, fixed_effects)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("The following required variables are missing in the data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Step 1: Residualize Treatment and Instrument
  data_resid <- residualize_D_Z(
    data = data,
    treatment_var = treatment_var,
    instrument_var = instrument_var,
    covariates = covariates,
    fixed_effects = fixed_effects
  )
  
  # Step 2: Residualize Outcome
  res_out <- residualize_outcome(
    data = data_resid,
    outcome_var = outcome_var,
    covariates = covariates,
    fixed_effects = fixed_effects
  )
  data_resid <- res_out$data
  model_Y <- res_out$model_Y
  
  # Step 3: Filter for Stable Region
  data_filtered <- filter_stable_region(
    data = data_resid,
    treatment_resid_var = "D_resid",
    epsilon = epsilon
  )
  
  # Step 4: Define Instrument Grid
  Z_min <- quantile(data_filtered$Z_resid, 0.05, na.rm = TRUE)
  Z_max <- quantile(data_filtered$Z_resid, 0.95, na.rm = TRUE)
  instrument_grid <- seq(Z_min, Z_max, length.out = 50)
  
  # Step 5: Compute LIV Estimates
  LIV_estimates <- compute_LIV_curve(
    data = data_filtered,
    instrument_grid = instrument_grid,
    outcome_resid_var = "Y_resid",
    treatment_resid_var = "D_resid",
    instrument_resid_var = "Z_resid",
    bw = bw,
    kernel = kernel
  )
  
  # Initialize diagnostics list
  diagnostics <- NULL
  
  # Step 6: Perform Diagnostics if Requested
  if (perform_diagnostics) {
    message("Performing diagnostics: Bandwidth selection and gradient analysis.")
    
    # a. Select Bandwidths using different methods
    bw_list <- select_bandwidths(
      formula = D_resid ~ Z_resid,
      data = data_filtered,
      regtype = "ll",
      bwmethods = diagnostic_methods
    )
    
    # b. Analyze Gradients for each bandwidth
    grad_list <- analyze_gradients(
      bw_list = bw_list,
      data = data_filtered,
      treatment_resid_var = "D_resid",
      instrument_resid_var = "Z_resid"
    )
    
    # c. Plot Gradient Distributions
    diagnostics_plots <- plot_gradient_distributions(
      grad_list = grad_list,
      instrument_grid = instrument_grid
    )
    
    diagnostics <- diagnostics_plots
  }
  
  # Step 7: Define Bootstrap Function
  compute_LIV_boot <- function(data, indices) {
    # Resample the data
    d_boot <- data[indices, ]
    d_boot <- d_boot %>% filter(abs(D_resid) > epsilon)
    
    # Compute LIV estimates
    LIV_boot_estimates <- compute_LIV_curve(
      data = d_boot,
      instrument_grid = instrument_grid,
      outcome_resid_var = "Y_resid",
      treatment_resid_var = "D_resid",
      instrument_resid_var = "Z_resid",
      bw = bw,
      kernel = kernel
    )
    
    return(LIV_boot_estimates)
  }
  
  # Step 8: Set Up Parallel Processing if Needed
  if (!is.null(workers)) {
    future::plan(future::multisession, workers = workers)
    message("Running bootstrapping in parallel using ", workers, " workers.")
  } else {
    future::plan(future::sequential)
  }
  
  # Step 9: Perform Bootstrapping
  set.seed(seed)
  boot_result <- boot::boot(
    data = data_filtered,
    statistic = compute_LIV_boot,
    R = R,
    parallel = ifelse(future::plan() == "multisession", "multicore", "no"),
    ncpus = ifelse(future::plan() == "multisession", workers, 1)
  )
  
  # Reset to sequential plan
  future::plan(future::sequential)
  
  # Step 10: Check for NAs in Bootstrap Results
  na_counts <- apply(boot_result$t, 2, function(x) sum(is.na(x)))
  cat("Number of NA bootstrap estimates per grid point:\n")
  print(na_counts)
  
  # Step 11: Compute Confidence Intervals
  ci_list <- purrr::map_dfr(seq_along(instrument_grid), function(i) {
    boot_estimates <- boot_result$t[, i]
    boot_estimates <- boot_estimates[!is.na(boot_estimates)]
    if (length(boot_estimates) > 0) {
      ci <- quantile(boot_estimates, probs = c(0.025, 0.975), na.rm = TRUE)
      tibble(
        Z_resid = instrument_grid[i],
        ci_lower = ci[1],
        ci_upper = ci[2]
      )
    } else {
      tibble(Z_resid = instrument_grid[i], ci_lower = NA_real_, ci_upper = NA_real_)
    }
  })
  
  # Step 12: Merge LIV Estimates with Confidence Intervals
  plot_data <- tibble(
    Z_resid = instrument_grid,
    LIV = LIV_estimates
  ) %>%
    left_join(ci_list, by = "Z_resid")
  
  # Step 13: Plot LIV Curve with Confidence Intervals
  p <- ggplot(plot_data, aes(x = Z_resid, y = LIV)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "lightblue") +
    labs(
      title = paste("Local Instrumental Variables (LIV) Curve for", outcome_var),
      x = "Residualized Instrument (Z_resid)",
      y = "LIV Estimate"
    ) +
    theme_minimal()
  
  # Step 14: Correlation Between D_resid and Z_resid
  correlation <- cor(data_filtered$D_resid, data_filtered$Z_resid, use = "complete.obs")
  cat("\nWithin-stable-region correlation between D_resid and Z_resid:", correlation, "\n\n")
  
  # Step 15: Generate Diagnostic Correlation Plot (if not already done)
  if (perform_diagnostics) {
    diagnostic_plot <- plot_correlation_D_Z(
      data = data_filtered,
      treatment_resid_var = "D_resid",
      instrument_resid_var = "Z_resid"
    )
  } else {
    # Optional: Generate a basic diagnostic plot
    diagnostic_plot <- ggplot(data_filtered, aes(x = Z_resid, y = D_resid)) +
      geom_point(alpha = 0.5, color = "darkgreen") +
      geom_smooth(method = "lm", color = "red") +
      labs(
        title = "Correlation Between Residualized D and Z",
        x = "Residualized Instrument (Z_resid)",
        y = "Residualized Treatment (D_resid)"
      ) +
      theme_minimal()
  }
  
  # Print the diagnostic plot
  print(diagnostic_plot)
  
  # Compile Results
  result_list <- list(
    outcome = outcome_var,
    model_Y = model_Y,
    LIV_estimates = LIV_estimates,
    instrument_grid = instrument_grid,
    ci_lower = ci_list$ci_lower,
    ci_upper = ci_list$ci_upper,
    plot = p
  )
  
  if (perform_diagnostics) {
    result_list$diagnostics <- diagnostics
  }
  
  return(result_list)
}
