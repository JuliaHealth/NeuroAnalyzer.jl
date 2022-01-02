#' @title ci95_p
#' 
#' @description Calculates 95%CI of a proportion
#' 
#' @param p Proportion
#' @param n Sample size
#' @export
#' 
#' @return Vector containing 2.5% and 97.5% boundaries

ci95_p <- function(p, n)
{
    ci.l <- p - qnorm(0.975) * sqrt(p * (1 - p) / n)
    ci.h <- p + qnorm(0.975) * sqrt(p * (1 - p) / n)
    return(c(ci.l, ci.h))
}