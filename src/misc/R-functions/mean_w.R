#' @title mean_w
#' 
#' @description Calculates weighted mean
#' 
#' @param x The vector to analyze
#' @param w The weights vector
#' @export
#' 
#' @return
 
mean.w <- function(x, w)
{
    sum(x * w) / sum(w)
}