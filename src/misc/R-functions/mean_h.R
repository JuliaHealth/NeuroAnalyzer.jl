#' @title mean_h
#' 
#' @description Calculates harmonic mean
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

mean_h <- function(x)
{
    length(x) / sum(1 / x)
}