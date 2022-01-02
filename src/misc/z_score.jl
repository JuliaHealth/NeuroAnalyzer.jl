#' @title z_score
#' 
#' @description Calculates Z-Score
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

z_score <- function(x)
{
    (x - mean(x)) / sd(x)
}