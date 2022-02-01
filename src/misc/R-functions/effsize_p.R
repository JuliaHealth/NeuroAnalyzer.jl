#' @title effsize_p
#' 
#' @description Calculates effect size for two proportions
#' 
#' @param p1 Proportion 1
#' @param p2 Proportion 2
#' @export
#' 
#' @return

effsize_h <- function(p1, p2)
{
    2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))
}