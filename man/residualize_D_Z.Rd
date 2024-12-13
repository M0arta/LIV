% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residualization.R
\name{residualize_D_Z}
\alias{residualize_D_Z}
\title{Residualize Treatment and Instrument Variables}
\usage{
residualize_D_Z(data, treatment_var, instrument_var, covariates, fixed_effects)
}
\arguments{
\item{data}{A data frame containing the dataset.}

\item{treatment_var}{A string specifying the treatment variable name.}

\item{instrument_var}{A string specifying the instrument variable name.}

\item{covariates}{A string vector specifying covariate names.}

\item{fixed_effects}{A string vector specifying fixed effect variables.}
}
\value{
A data frame with residualized treatment and instrument variables.
}
\description{
Residualizes the treatment (D) and instrument (Z) variables by regressing them on specified covariates and fixed effects.
}
