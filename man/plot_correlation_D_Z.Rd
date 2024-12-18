% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{plot_correlation_D_Z}
\alias{plot_correlation_D_Z}
\title{Plot Correlation Between Residualized D and Z}
\usage{
plot_correlation_D_Z(data, treatment_resid_var, instrument_resid_var)
}
\arguments{
\item{data}{A data frame containing residualized treatment and instrument variables.}

\item{treatment_resid_var}{A string specifying the name of the residualized treatment variable.}

\item{instrument_resid_var}{A string specifying the name of the residualized instrument variable.}
}
\value{
A ggplot object visualizing the correlation.
}
\description{
Generates a scatter plot with a fitted regression line to assess the correlation between residualized treatment and instrument variables.
}
\examples{
\dontrun{
diagnostic_plot <- plot_correlation_D_Z(
  data = mroz_filtered,
  treatment_resid_var = "D_resid",
  instrument_resid_var = "Z_resid"
)
print(diagnostic_plot)
}
}
