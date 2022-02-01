#' @title z_prob
#' 
#' @description Calculates z-test for two probabilities
#' 
#' @param p1 Probability 1
#' @param p2 Probability 2
#' @param n1 Sample size 1
#' @param n2 Sample size 2
#' @export
#' 
#' @return

z_prob <- function(p1, p2, n1, n2)
{
    numerator <- (p1 / n1) - (p2 / n2)
    p.common <- (p1 + p2) / (n1 + n2)
    denominator <- sqrt(p.common * (1 - p.common) * (1 / n1 + 1 / n2))
    numerator / denominator
}