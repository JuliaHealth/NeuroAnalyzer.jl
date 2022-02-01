#' @title prank
#' 
#' @description Calculates the prediction interval (95% CI adjusted for sample size)
#' 
#' @param n Sample size
#' @export
#' 
#' @return

pred_int <- function(n)
{
    if (n == 0)
    {
        stop("n cannot be 0")
        return()
    }
    pred_int <- 0
    if (n == 25) pred_int <- 2.10
    if (n == 30) pred_int <- 2.08
    if (n == 35) pred_int <- 2.06
    if (n == 40) pred_int <- 2.05
    if (n == 50) pred_int <- 2.03
    if (n == 60) pred_int <- 2.02
    if (n == 70) pred_int <- 2.01
    if (n == 80) pred_int <- 2.00
    if (n == 90) pred_int <- 2.00
    if (n == 100) pred_int <- 1.99
    if (n == 200) pred_int <- 1.98
    if (n > 0 & n < 21)
    {
        pred_int.v <- c(NaN, 15.56, 4.97, 3.56, 3.04, 2.78, 2.62, 2.51, 2.43, 2.37, 2.33, 2.29, 2.26, 2.24, 2.22, 2.18, 2.17, 2.16, 2.10)
        pred_int <- pred_int.v[n]
    }
    if (pred_int == 0)
    {
        warning("Result may not be accurate")
        pred_int <- 1.96
    }
    return(pred_int)
}