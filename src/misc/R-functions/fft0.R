#' @title fft0
#' 
#' @description Zero-padded FFT
#' 
#' @param x The vector for FFT
#' @param padlength Number of zeros to pad
#' @export
#' 
#' @return

fft0 <- function(x, padlength)
{
    fft(c(x, rep(0, padlength)))
}