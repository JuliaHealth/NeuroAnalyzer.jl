#' @title moe
#' 
#' @description Calculates margin of error for given sample size
#' 
#' @param n Sample size
#' @export
#' 
#' @return

moe <- function(n)
{
    1 / sqrt(n)
}