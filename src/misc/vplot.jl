#' @title vplot
#' 
#' @description Plots 2d vector(s) (e.g. c(2,2)) with (0,0) as the origin
#' 
#' @param v The vector to draw
#' @param l The X- and Y-axis min and max value, default 10
#' @export
#' 
#' @return 

library(ggplot2)

vplot <- function(v, l = 10)
{
    vx <- numeric(length(v) / 2)
    vy <- numeric(length(v) / 2)

    for (i in 1:(length(v) / 2))
    {
        vx[i] <- v[2 * i - 1]
        vy[i] <- v[2 * i]
    }

    ggplot() + 
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_segment(aes(x = 0, y = 0, xend = vx, yend = vy)) +
    scale_x_continuous("x", limits = c(-l, l), expand = c(0, 0)) +
    scale_y_continuous("y", limits = c(-l, l), expand = c(0, 0))
}