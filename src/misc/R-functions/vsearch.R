#' @title vsearch
#' 
#' @description Gets the index of particular value
#' 
#' @param x The vector to analyze
#' @param y The value to search
#' @export
#' 
#' @return

vsearch <- function(x, y)
{
    which(abs(x - y) == min(abs(x - y)))
}