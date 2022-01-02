#' @title binom.eq
#' 
#' @description Calculates probability of a proportion of success for a given number of successes in a given number of trials
#' 
#' @param p     Proportion of successes
#' @param r     Number of successes
#' @param n     Number of trials
#' @export
#' 
#' @return

binom.eq <- function(p, r, n)
{
    (factorial(n)/(factorial(n) * factorial(n-r))) * (p^r) * ((1 - p)^(n - r))
}