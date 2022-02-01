#' @title se_2p
#' 
#' @description Calculates standard error of two proportions
#' 
#' @param p1 Proportion 1
#' @param p2 Proportion 2
#' @param n Sample size 
#' @export
#' 
#' @return

se_2p <- function(p1, p2, n)
{
    sqrt((p1 * p2) / n)
}