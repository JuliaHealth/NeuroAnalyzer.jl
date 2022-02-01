#' @title ci95_2mean
#' 
#' @description Calculates 95%CI for a difference in means
#' 
#' @param x The vector to analyze
#' @param y The vector to analyze
#' @export
#' 
#' @return Vector containing 2.5% and 97.5% boundaries

ci95_2mean <- function(x, y)
{
    #calculate pooled variance
    pvar = ((length(x) - 1) * sd(x)^2 + (length(y) - 1) * sd(y)^2) / (length(x) + length(y) - 2)
    
    #calculate margin of error
    margin <- qt(0.975, df = length(x) + length(y) - 1) * sqrt(pvar / length(x) + pvar / length(y))

    ci.l <- (mean(x) - mean(y)) - margin
    ci.h <- (mean(x) - mean(y)) + margin
    return(c(ci.l, ci.h))
}