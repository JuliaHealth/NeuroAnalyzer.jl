#' @title binom_stat
#' 
#' @description Calculates mean and standard deviation for a given probability
#' 
#' @param p Probability
#' @param n Number of events
#' @export
#' 
#' @return A vector containing mean and standard deviation

binom_stat <- function(p, n)
{
    binom_mean <- n * p
    binom_sd <- sqrt(n * p * (1 - p))
    return(c(binom_mean, binom_sd))
}