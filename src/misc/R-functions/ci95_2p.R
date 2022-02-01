#' @title ci95_2p
#' 
#' @description Calculates 95%CI for a difference in proportions
#' 
#' @param p1 Proportion 1
#' @param p2 Proportion 2
#' @param n1 Sample size 1
#' @param n2 Sample size 2
#' @export
#' 
#' @return Vector containing 2.5% and 97.5% boundaries

ci95_2p <- function(p1, n1, p2, n2)
{
    ci.l <- (p1 - p2) - qnorm(0.975) * sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    ci.h <- (p1 - p2) + qnorm(0.975) * sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    return(c(ci.l, ci.h))
}