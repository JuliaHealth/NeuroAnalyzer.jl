#' @title se_2p_diff
#' 
#' @description Calculates standard error for the difference in two proportions
#' 
#' @param p1 Proportion 1
#' @param p2 Proportion 2
#' @param n1 Sample size 1
#' @param n2 Sample size 1
#' @export
#' 
#' @return

se_2p_diff <- function(p1, p2, n1, n2)
{
    sqrt(((p1 * (1 - p1)) / n1) + ((p2 * (1 - p2)) / n2))
}