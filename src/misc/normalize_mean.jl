#' @title normalize_mean
#' 
#' @description Scales around the mean
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return Normalized vector

normalize_mean <- function(x)
{
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}