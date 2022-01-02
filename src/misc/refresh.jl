#' @title refresh
#' 
#' @description Sources all local R-functions
#' 
#' @param fdir Path to R-functions directory, default ~/Documents/Code/R-functions
#' @export
#' 
#' @return

refresh <- function(fdir = "~/Documents/Code/R-functions")
{
    current_dir <- getwd()
    setwd(fdir)
    rfunctions <- list.files(pattern = "*.R")
    for (f in rfunctions)
    {
        source(f)
    }
    setwd(current_dir)
    cat("Reloaded", length(rfunctions), "functions..\n")
}