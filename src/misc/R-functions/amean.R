#' @title amean
#' 
#' @description Calculates mean for a multi-dimensional x×y×z×… array
#' 
#' @param a The array to analyze 
#' @export
#' 
#' @return m Matrix m×n

amean <- function(a)
{
    m <- matrix(0, nrow = nrow(a), ncol = ncol(a))
    z = dim(a)[3]
    for (x in 1:nrow(a))
    {
        for (y in 1:ncol(a))
        {
            m[x, y] = mean(a[x,y,1:z], na.rm = TRUE)
        }
    }
    return(as.matrix(m))
}