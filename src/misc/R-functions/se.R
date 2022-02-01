#' @title se
#' 
#' @description Calculates standard error
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

sem <- function(x)
{
    sd(x) / sqrt(length(x))
}