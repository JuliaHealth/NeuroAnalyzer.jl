#' @title current_density
#' 
#' @description Calculates current density for tDCS stimulation [A/m2^2]
#' 
#' @param current Current in mA
#' @param area Electrode size in cm^2
#' @export
#' 
#' @return

current_density <- function(current, area)
{
    # convert mA to A
    current = current / 1000

    # convert cm^2 to m^2
    area = area / 10000

    current / area
}