#' @title binom_prob
#' 
#' @description Calculates the probability of a proportion for a given number of times in a number of trials
#' 
#' @param p Probability
#' @param y Number of times
#' @param n Number of trials
#' @export
#' 
#' @return

binom_prob <- function(p, y, n)
{
    (factorial(n) / factorial(y) * factorial(n - y)) * p^y * (1 - p)^(n - y)
}