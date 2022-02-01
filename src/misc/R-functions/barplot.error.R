#' @title barplot.error
#' 
#' @description Plot barplot with error bars
#' 
#' @param x		Factor object
#' @param y		Number of successes
#' @param ...	Additional arguments to be passed to barplot function
#' @export
#' 
#' @return

barplot.error <- function(x, y, ...)
{
    mod <- lm(y ~ x)
    reps <- sqrt(length(y) / length(levels(x)))
    sem <- sigma(mod) / reps
    means <- tapply(y, x, mean)
    upper <- max(means) + sem
    lev <- levels(x)
    barpl <- barplot(means, ...)
    invisible(sapply(1:length(barpl),
                     function(i) arrows(barpl[i],
                                        means[i] + sem,
                                        barpl[i],
                                        means[i] - sem,
                                        angle = 90,
                                        code = 3,
                                        length = 0.08)))
}