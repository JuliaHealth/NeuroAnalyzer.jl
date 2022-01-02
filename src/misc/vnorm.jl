#' @title vnorm
#' 
#' @description Returns a normalized vector: x / |x|
#' 
#' @param x The vector to normalize 
#' @param L Type of norm 2: Euclidean (default), 1: Manhattan, Inf: Max, 3: RMS, 4: x / |x|
#' @export
#' 
#' @return Normalized vector

vnorm <- function(x, L = 2)
{
    if (L == 2)
    {
        n <- magnitude(x)
    }
    if (L == 1)
    {
        n <- sum(abs(x))
    }
    if (L == Inf)
    {
        n <- max(abs(x))
    }
    if (L == 3)
    {
        n <- as.vector(magnitude(x) / sqrt(length(x)))
    }
    if (L == 4)
    {
        n <- x / as.vector(magnitude(x))
    }
    return(n)
}