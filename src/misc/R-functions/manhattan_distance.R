#' @title manhattan_distance
#' 
#' @description Calculates Manhattan distance between two vectors
#' 
#' @param x The vector to analyze
#' @param y The vector to analyze
#' @export
#' 
#' @return

manhattan_distance <- function(x, y)
{
    sum(x - y)
}