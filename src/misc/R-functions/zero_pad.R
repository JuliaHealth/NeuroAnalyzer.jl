#' @title z_score
#' 
#' @description Pads the matrix M with zeros to make M square
#' 
#' @param M The matrix to pad
#' @export
#' 
#' @return Padded matrix

zero_pad <- function(M)
{
    nr <- nrow(M)
    nc <- ncol(M)
    if (nr > nc)
    {
        Mp <- matrix(0, nrow = nr, ncol = nr - nc)
        Mp <- cbind(M, Mp)
    }
    else if (nr < nc)
    {
        Mp <- matrix(0, nrow = nc - nr, ncol = nc)
        Mp <- rbind(M, Mp)
    }
    else if (nr == nc)
    {
        Mp <- M
    }
    return(Mp)
}