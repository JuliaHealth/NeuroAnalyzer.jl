#' @title mean_g
#' 
#' @description Calculates geometric mean
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

mean_g <- function(x)
{
    exp(mean(log(x[x > 0])))
}