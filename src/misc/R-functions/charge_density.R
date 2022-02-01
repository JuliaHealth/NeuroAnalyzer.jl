#' @title charge_density
#' 
#' @description Calculates charge density for tDCS stimulation [A×s/m2^2]
#' 
#' @param current Current in mA
#' @param area Electrode size in cm^2
#' @param time Time of stimulation in s
#' @export
#' 
#' @return

charge_density <- function(current, area, time)
{
    # convert mA to A
    current = current / 1000

    # convert cm^2 to m^2
    area = area / 10000

    (current * time) / area
}