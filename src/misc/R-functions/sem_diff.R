#' @title sem_diff
#' 
#' @description Calculates SEM (standard error of the mean) for the difference of two means
#' 
#' @param x The vector to analyze
#' @param y The vector to analyze
#' @export
#' 
#' @return

sem_diff <- function(x, y)
{
    sqrt((sd(x)^2 / sqrt(length(x))) + (sd(y)^2 / sqrt(length(y))))
}