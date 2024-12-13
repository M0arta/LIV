# R/select_bandwidth.R

#' Select Bandwidths Using Different Methods
#'
#' Selects bandwidth for non-parametric regression using various bandwidth selection methods.
#'
#' @param formula A formula specifying the relationship between the residualized treatment and instrument (e.g., \code{D_resid ~ Z_resid}).
#' @param data A data frame containing the necessary variables.
#' @param regtype Type of regression. Default is \code{"ll"} for local linear.
#' @param bwmethods A character vector specifying bandwidth selection methods. Options include \code{"default"}, \code{"cv.aic"}, \code{"cv.ls"}.
#' @return A named list of bandwidth objects corresponding to each selected method.
#' @examples
#' \dontrun{
#' bw_list <- select_bandwidth(
#'   formula = D_resid ~ Z_resid,
#'   data = all_tot,
#'   regtype = "ll",
#'   bwmethods = c("default", "cv.aic", "cv.ls")
#' )
#' }
#' @import np
#' @export
select_bandwidth <- function(formula, data, regtype = "ll", bwmethods = c("default", "cv.aic", "cv.ls")) {
  # Validate bwmethods
  allowed_methods <- c("default", "cv.aic", "cv.ls")
  invalid_methods <- setdiff(bwmethods, allowed_methods)
  if (length(invalid_methods) > 0) {
    stop(paste("Invalid bandwidth methods:", paste(invalid_methods, collapse = ", "),
               ". Allowed methods are:", paste(allowed_methods, collapse = ", "), "."))
  }
  
  # Initialize list to store bandwidth objects
  bw_list <- list()
  
  # Iterate over each method and select bandwidth
  for (method in bwmethods) {
    if (method == "default") {
      bw <- npregbw(formula = formula, data = data, regtype = regtype)
      bw_list[["default"]] <- bw
    } else if (method == "cv.aic") {
      bw <- npregbw(formula = formula, data = data, regtype = regtype, bwmethod = "cv.aic")
      bw_list[["cv.aic"]] <- bw
    } else if (method == "cv.ls") {
      bw <- npregbw(formula = formula, data = data, regtype = regtype, bwmethod = "cv.ls")
      bw_list[["cv.ls"]] <- bw
    }
  }
  
  return(bw_list)
}
