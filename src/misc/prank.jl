#' @title prank
#' 
#' @description Calculates percentile rank
#' 
#' @param x The vector to analyze
#' @param y The value to rank as a percentile
#' @export
#' 
#' @return

prank <- function(x, y)
{
    length(x[x <= y]) / length(x) * 100
}