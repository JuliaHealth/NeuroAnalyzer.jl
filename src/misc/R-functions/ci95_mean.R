#' @title ci95_mean
#' 
#' @description Calculates 95%CI of mean
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return Vector containing 2.5% and 97.5% boundaries

ci95_mean <- function(x)
{
    ci.l <- mean(x) - qt(0.975, df = length(x)-1) * sd(x) / sqrt(length(x))
    ci.h <- mean(x) + qt(0.975, df = length(x) - 1) * sd(x) / sqrt(length(x))
    return(c(ci.l, ci.h))
}