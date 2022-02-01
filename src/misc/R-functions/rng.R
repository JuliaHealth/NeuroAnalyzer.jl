#' @title rng
#' 
#' @description Calculates range
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

rng <- function(x)
{
    max(x) - min(x)
}