#' @title magnitude
#' 
#' @description Calculates Euclidean norm (magnitude) of the vector or matrix
#' 
#' @param x The vector or matrix to analyze
#' @export
#' 
#' @return

magnitude <- function(x)
{
    if (is.vector(x) | is.ts(x))
    {
        m <- as.vector(sqrt(crossprod(x)))
    }
    else if (is.matrix(x))
    {
        m <- norm(x)
    }
    return(m)
}