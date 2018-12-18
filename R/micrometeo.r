#' get roughness length. From page 123 in Moene and van Dam
#' @param surface is one of the following strings: water, smooth ice, rough ice, snow, soils, short grass, long grass, low mature crops, high mature crops, continuous bushland, mature pine forest, tropical forest, deciduous forest, dense low buildings, regularly built town
#' @export
get_z0 <- function(surface)
{
    return z0$z0[surface==z0$surface]
}

#' get displancement height. From page 123 in Moene and van Dam
#' @param surface is one of the following strings: water, smooth ice, rough ice, snow, soils, short grass, long grass, low mature crops, high mature crops, continuous bushland, mature pine forest, tropical forest, deciduous forest, dense low buildings, regularly built town
#' @export
get_d <- function(surface)
{
    d <- z0$d[surface==z0$surface]
    d[is.na(d)] <- 0
    return(d)
}

#' accelaration of gravity
#' @export
g_acceleration <- function()
{
    return(9.8)
}

#' Specific flux of virtual sensible heat
#' from Katul and Parlange 1992 WWR
#' @export
virtual_H <- function(H,LE,Ta,q)
{
    return(H+0.61*Ta*spec_heat_cp_air(q)*LE)
}

#' Obukhov length L
#' from Katul and Parlange 1992 WWR
#' @export
length_L <- function(ustar, Hv, Ta,q)
{
    num=ustar^3
    den=vonKarman()*g_acceleration()*(Hv/(mass_density_air()*spec_heat_cp_air(q)*Ta))
    return(-num/den)
}

#' calculate aerodynamic resistance for natural surfaces
#' from http://www.fao.org/docrep/x0490e/x0490e00.htm, Allen, Pereira, Raes, Smith, Crop evapotranspiration - Guidelines for computing crop water requirements - FAO Irrigation and drainage paper 56, FAO, Rome, 1998
#' you can check z0 and d in Moene and Van Dam, Transport in the Atmosphere-Vegetation-Soil Continuum, Cambridge University Press, page 123
#' @param uz is the wind speed at height zm [m.s^-1]
#' @param zm is the height of the wind speed measurements
#' @param zh is the height of the humidity measurements
#' @param d is the displacement height
#' @param z0 is the aerodynamic roughness length
#' @param h is the height of the reference crop, ie clipped grass
#' @importFrom micrometeo vonKarman
#' @export
res_aero <- function(uz,zh=2,zm=2,d=NULL,z0=NULL,h=0.12,urban=FALSE)
{
    if(urban)
    {
        return(1/(65*0.001)) # s.m-1, from Grimmond, C. S. B., and T. R. Oke. “Aerodynamic Properties of Urban Areas Derived from Analysis of Surface Form.” Journal of Applied Meteorology 38, no. 9 (September 1, 1999): 1262–92. https://doi.org/10.1175/1520-0450(1999)038<1262:APOUAD>2.0.CO;2.
    } else
    {
        if(is.null(d) == FALSE & is.null(z0) == FALSE)
        {
            zom=z0
        } else
        {
            d <- h*2/3
            zom <- 0.123*h
        }

        zoh <- 0.1*zom
        num <- log((zm-d)/zom)*log((zh-d)/zoh)
        den <- (vonKarman()^2)*uz
        ra <- num/den
        return(ra)
    }
}


#' calculate stomatal resistance for reference crop, ie clipped grass
#' @param rl bulk stomatal resistance of well illuminated leaf [s.m^-1]
#' @param LAIact active (sunlit) leaf area index [m^2 (leaf area) . m^-2 (soil surface)]is the wind speed at height h [m.s^-1]
#' @param h is the height of the canopy
#' @export
res_stom <- function(h=0.12,rl=100,LAI=NULL)
{
    if(is.null(LAI))
    {
        LAI <- 24*h # for clipped grass
    }
        LAIact <- 0.5*LAI #active (sunlit) leaf area index [m^2 (leaf area) . m^-2 (soil surface)]is the wind speed at height h [m.s^-1]
        rs <- rl/LAIact
    return(rs)
}
