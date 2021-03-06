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

#' calculates the potential et from mowene and van dam
#' @param T mean daily temperature [K]
#' @param q mean daily specific humidity [kg.kg^-1]
#' @param p mean daily air pressure at reference height [kPa]
#' @param vpd vapour pressure deficit at reference height [kPa]
#' @param Rn net radiation [W.m^-2]
#' @param G ground heat flux [W.m^-2]
#' @importFrom micrometeo latent_heat_vap slope_sat_vpressure psychr spec_heat_cp_air mass_density_air
#' @export
penmon_day <- function(T,q,p,vpd,Rn,G,ra,rs)
{
    num1 = slope_sat_vpressure(T)*(Rn - G)
    num2 = vpd *mass_density_air(T,q,p)*spec_heat_cp_air(q)/ra

    den = slope_sat_vpressure(T) + psychr(p)*(1+rs/ra)

    E = 24*60*60*(num1+num2)/den
    return(E/latent_heat_vap())
}


#' calculates the Penman evapotranspiration from a water surface. from Moene and Van Dam, page 258
#' @param T mean daily temperature [K]
#' @param q mean daily specific humidity [kg.kg^-1]
#' @param p mean daily air pressure at reference height [kPa]
#' @param vpd vapour pressure deficit at reference height [kPa]
#' @param Rn net radiation [W.m^-2]
#' @param G ground heat flux [W.m^-2]
#' @importFrom micrometeo latent_heat_vap slope_sat_vpressure psychr spec_heat_cp_air mass_density_air
#' @export
penman <- function(T,q,p,vpd,Rn,G,ra)
{
    den <- slope_sat_vpressure(T) + psychr(p)
    num1 <- slope_sat_vpressure(T) * (Rn-G)
#    num2 <- vpd*mass_density_air(T,q,p)*spec_heat()/ra

    num2 <- vpd*mass_density_air(T,q,p)*spec_heat_cp_air(q)/ra

    et <- 24*60*60*(num1/den+num2/den)/latent_heat_vap() # convert latent heat flux into [mm]
    return(et)
}

#' Solves the Penman
