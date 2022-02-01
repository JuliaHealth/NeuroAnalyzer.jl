#' @title ci95_boot
#' 
#' @description Calculates 95%CI of mean using bootstrapping
#' 
#' @param x The vector to analyze
#' @param r Sampling number (default r = 1000)
#' @export
#' 
#' @return Vector containing 2.5% and 97.5% boundaries

ci95_boot <- function(x, r = 1000)
{
    mn <- numeric()
  
    for (i in 1:r)
    {
        s <- sample(x, replace = TRUE)
        mn[i] <- mean(s)
    }
    ci.l <- quantile(mn, 0.025)
    ci.h <- quantile(mn, 0.975)
    return(c(ci.l, ci.h))    
}