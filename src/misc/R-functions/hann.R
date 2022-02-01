#' @title hann
#' 
#' @description Returns the N-point symmetric Hann window in a column vector
#' 
#' @param n The length of the window to generate 
#' @export
#' 
#' @return

hann <- function(n)
{
    0.5 * (1 - cos(2 * pi * seq(0, 1, length = n)))
}