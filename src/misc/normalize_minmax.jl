#' @title normalize_minmax
#' 
#' @description Normalizes to 0…1
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return Normalized vector

normalize_minmax <- function(x)
{
    (x - min(x)) / (max(x) - min(x))
}