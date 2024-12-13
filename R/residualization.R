# R/residualization.R

#' Residualize Outcome Variable Y
#'
#' Residualizes the outcome variable by regressing it on specified covariates and fixed effects.
#'
#' @param data A data frame containing the dataset.
#' @param outcome_var A string specifying the outcome variable name.
#' @param covariates A string vector specifying covariate names.
#' @param fixed_effects A string vector specifying fixed effect variables.
#' @return A list containing the modified data frame with residuals and the fitted model.
#' @import dplyr
#' @import fixest
#' @export
residualize_outcome <- function(data, outcome_var, covariates, fixed_effects) {
  # Create the formula dynamically
  covariates_formula <- paste(covariates, collapse = " + ")
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")
  
  formula_Y <- as.formula(paste0(
    outcome_var, 
    " ~ ", 
    covariates_formula, 
    " | ", 
    fixed_effects_formula
  ))
  
  model_Y <- feols(formula_Y, data = data)
  data <- data %>% mutate(Y_resid = resid(model_Y))
  return(list(data = data, model_Y = model_Y))
}



# R/residualization.R

#' Residualize Treatment and Instrument Variables
#'
#' Residualizes the treatment (D) and instrument (Z) variables by regressing them on specified covariates and fixed effects.
#'
#' @param data A data frame containing the dataset.
#' @param treatment_var A string specifying the treatment variable name.
#' @param instrument_var A string specifying the instrument variable name.
#' @param covariates A string vector specifying covariate names.
#' @param fixed_effects A string vector specifying fixed effect variables.
#' @return A data frame with residualized treatment and instrument variables.
#' @import dplyr
#' @import fixest
#' @export
residualize_D_Z <- function(data, treatment_var, instrument_var, covariates, fixed_effects) {
  # Residualize Treatment (D)
  covariates_formula <- paste(covariates, collapse = " + ")
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")
  
  formula_D <- as.formula(paste0(
    treatment_var, 
    " ~ ", 
    covariates_formula, 
    " | ", 
    fixed_effects_formula
  ))
  
  model_D <- feols(formula_D, data = data)
  data <- data %>% mutate(D_resid = resid(model_D))
  
  # Residualize Instrument (Z)
  formula_Z <- as.formula(paste0(
    instrument_var, 
    " ~ ", 
    covariates_formula, 
    " | ", 
    fixed_effects_formula
  ))
  
  model_Z <- feols(formula_Z, data = data)
  data <- data %>% mutate(Z_resid = resid(model_Z))
  
  return(data)
}
