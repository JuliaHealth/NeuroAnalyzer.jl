#' @title ci95_median
#' 
#' @description Calculates 95%CI of median
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return Vector containing 2.5% and 97.5% boundaries

ci95_median <- function(x)
{
    ci.l <- median(x) - 1.58 * (IQR(x) / sqrt(length(x)))
    ci.h <- median(x) + 1.58 * (IQR(x) / sqrt(length(x)))
    return(c(ci.l, ci.h))
}