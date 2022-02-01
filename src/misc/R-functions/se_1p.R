#' @title se_1p
#' 
#' @description Calculates standard error of one proportions
#' 
#' @param p Proportion
#' @param n Sample size 
#' @export
#' 
#' @return

se_1p <- function(p, n)
{
    sqrt((p * (1 - p)) / n)
}