#' @title euclidean_distance
#' 
#' @description Calculates Euclidean distance between two vectors
#' 
#' @param x The vector to analyze
#' @param y The vector to analyze
#' @export
#' 
#' @return

euclidean_distance <- function(x, y)
{
    sqrt(sum((x - y)^2))
}