#' @title angle
#' 
#' @description Returns the phase angles, in radians, of a matrix with complex elements
#' 
#' @param h The vector to analyze
#' @export
#' 
#' @return 

angle <- function(h)
{
    atan2(Im(h), Re(h))
}