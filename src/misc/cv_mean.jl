#' @title cv_mean
#' 
#' @description Calculates coefficient of variation for a mean
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

cv_mean <- function(x)
{
    sd(x) / mean(x)
}