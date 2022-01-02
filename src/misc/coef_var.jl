#' @title coef_var
#' 
#' @description Calculates coefficient of variation
#' 
#' @param x The vector to analyze
#' @export
#' 
#' @return

coef_var <- function(x)
{
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}