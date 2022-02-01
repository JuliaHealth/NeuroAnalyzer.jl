#' @title rms
#' 
#' @description Calculates root mean square
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

rms <- function(x)
{
    sqrt(mean(x^2))
}