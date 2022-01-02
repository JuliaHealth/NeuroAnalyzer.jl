#' @title cv
#' 
#' @description Calculates coefficient of variation for a statistic
#' 
#' @param se Standard error
#' @param s  Statistics (e.g. mean)
#' @export
#' 
#' @return

cv <- function(se, s)
{
    100 * (se / s)
}