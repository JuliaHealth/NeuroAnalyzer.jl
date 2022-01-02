#' @title cv_median
#' 
#' @description Calculates coefficient of variation for a median
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

cv_median <- function(x)
{
    cv <- ((quantile(x, 0.75) - quantile(x, 0.25)) / 2) / median(x)
    return(cv[[1]])
}