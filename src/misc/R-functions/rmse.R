#' @title rmse
#' 
#' @description Calculates root mean square error
#' 
#' @param x The vector to analyze
#' @param y The vector to analyze
#' @export
#' 
#' @return

rmse <- function(x, y)
{
    sqrt(mean((x - y)^2))
}