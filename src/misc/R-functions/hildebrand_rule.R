#' @title hildebrand_rule
#' 
#' @description Calculates Hildebrand rule for symmetry
#' @description H < 0.2 means symmetry
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

hildebrand_rule <- function(x)
{
    (mean(x) - median(x))/sd(x)
}