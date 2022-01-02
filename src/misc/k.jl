#' @title k
#' 
#' @description Calculates number of categories for a given sample size
#' 
#' @param n Sample size
#' @export
#' 
#' @return

k <- function(n)
{
    c(sqrt(n), (1 + 3.222 * log10(n)))
}