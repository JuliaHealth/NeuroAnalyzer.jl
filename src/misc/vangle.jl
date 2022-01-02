#' @title vangle
#' 
#' @description Calculates the angle between two vectors
#' 
#' @param x The vector to analyze
#' @param y The vector to analyze
#' @export
#' 
#' @return List containing angle in radians and degrees

vangle <- function(x, y)
{
    angle_r <- acos(x %*% y / (vnorm(x) * vnorm(y)))
    angle_d <- angle_r * (180 / pi)
    return(list(angle_rad = as.vector(angle_r), angle_deg = as.vector(angle_d)))
}