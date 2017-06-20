
#' thermal diffusivity in m^2/s. Table in page 50 of Moene and van Dam, "Transport in the Atmosphere-Vegetation-Soil Continuum "
#' @export
thermal_diff <- function()
{
    return(1.2*10^-6)
}


#' pressure used for calculating psychrometric constant in kPa
#' @export
p_coruche <- function()
{
    return(101.5)
}

#' latent heat of vaporization (lambda) at normal pressure MJ/kg
#' @export
latent_heat_vap <- function()
{
    return(2.45) 
}

#' ratio molecular weight of water vapour/dry air = 0.622
#' @export
ratio_mol_w <- function()
{
    return(0.622)
}

#' specific heat cp at constant pressure MJ*kg^-1*degreeC^-1
#' @export
spec_heat <- function()
{

    return(1.013*10^-3)
}

#' mass density of water (rho) in kg/m3 
#' @export
mass_density_h2o <- function()
{
    return(999) 
}
 

#' psychrometric constant in kPa/K for local pressure conditions
#' check http://www.fao.org/docrep/X0490E/x0490e07.htm#psychrometric%20constant%20(g)
#' @export
psychr <- function(p)
    {
        psy <- spec_heat()*p/(ratio_mol_w()*latent_heat_vap())
        return(psy)
    }

#' vapour pressure at saturation in kPa after Moene and Van Dam Cambridge University Press 2014 "Transport in the Atmosphere-Vegetation-Soil Continuum " page 353 eq B.20
#' @param T in [K]
#' @param esat in [kPa]
#' @export
sat_vpressure <- function(T)
{
    T <- T+273.15
    esat <- 0.6112*exp(17.62*(T-273.15)/(-0.53+T))
    return(esat)
}

#' vapour pressure in Pa after Moene and Van Dam Cambridge University Press 2014 "Transport in the Atmosphere-Vegetation-Soil Continuum " page 353 eq B.19
#' @param q in [kg.kg^-1]
#' @param p in kPa
#' @param e in kPa
#' @export
specific_hum2vpressure <- function(q,p)
{
    e <- q*p*8/5
    return(e)
}

#' vapour pressure in kPa from wet bulb after Moene and Van Dam Cambridge University Press 2014 "Transport in the Atmosphere-Vegetation-Soil Continuum " page 354
#' @param Twet in [K]
#' @param Tdry in [K]
#' @param p in kPa
#' @export
wetbrh2ulb2vpressure <- function(Twet,Tdry,p)
    {
        e <- sat_vpressure(Twet)-psychr(p)*(Tdry-Twet)     #psychr is in kPa/K
        return(e)
    }



#' relative humidity after Moene and Van Dam Cambridge University Press 2014 "Transport in the Atmosphere-Vegetation-Soil Continuum " page 352
#' @param q in [kg.kg^-1]
#' @param p in Pa
#' @param T in [K]
#' @export
specific_hum2rh <- function(q,T,p)
{    
    return(specific_hum2vpressure(q,p)/sat_vpressure(T))
}

#' relative humidity to vapour pressure
#' @param rh in [-]
#' @param T in [K]
#' @export
rh2vpressure <- function(rh,T)
{    
    return(rh*sat_vpressure(T))
}


#' bowen ratio
#' @export
bowen <- function(Ta2,Ta1,ev2,ev1,p)
{
    B <- psychr(p)*(Ta1-Ta2)/(ev1-ev2)
    return(B)
}

#' ET from latent heat flux (le) in mm
#' @export
ET <- function(le,rho)
{
    et <- le/(rho*latent_heat_vap())
    return(et)
}

#' latent heat flux from the energy balance
#' @export
LE <- function(R,dQdt,G,B)
{
    le <- (R-G-dQdt)/(1+B)
    return(le)
}

#' sensible heat flux from bowen ratio 
#' @export
H <- function(le,B)
{
    h <- B*le
    return(h)
}



#' virtual temperature after Arya, S. P. (2001). Introduction to Micrometeorology. International Geophysics Series (Vol. 79, p. 420). San Diego: Academic Press. doi:10.4043/13298-MS)
#' @export
Tvirtual <- function(T,Q) #virtual temperature
    {
        Tv <- T*(1+0.61*Q)
        return(Tv)
    }

#' specific humidity Mw/(Md+Mw) where Md is mass of dry air and Mw is mass of water vapour. to a good approximation Q=mw*ev/md*p where ma and mw are the mean molecular masses. After Arya, S. P. (2001) Introduction to Micrometeorology. International Geophysics Series (Vol. 79, p. 420). San Diego: Academic Press. doi:10.4043/13298-MS
#' @param ev is water vapour pressure
#' @param p is total pressure.
#' @export
Q <- function(ev,p) 
    {
        Q <- 0.622*ev/p # approximation 
        return(Q)
    }

#' potential temperature,  k=R/cp, after Arya, S. P. (2001). Introduction to Micrometeorology. International Geophysics Series (Vol. 79, p. 420). San Diego: Academic Press. doi:10.4043/13298-MS))
#' @export
O <- function(T,p) 
    {
        p0 <- 1000 #in mbar
        k <- 0.286 
        O <- T*(p0/p)^k
        return(O)
    }

#' volumetric heat capacity of a soil. From the encyclopedia of soil science, page 306. Unit is J.m^-3.K^-1
#' @param xw is the volumetric water content
#' @param xorg is the organic carbon content of the soil
#' @param xsolid is the mineral content of the soil
#' @export
vol_heat_capacity_soil <- function(xw,xorg,xsolid) #this is the volumetric heat capacity!!! organic content in soil xorg and water content xw and mineral content xsolid
    {
        Cs <- 2.45*xsolid+2.45*xorg+4.186*xw 
        return(Cs*1e6)
    }



richardson <- function(Ri_in)
    {
        #gradient richardson number
        Ov1 <- Ri_in$Ov1
        Ov2 <- Ri_in$Ov2
        v1 <- Ri_in$v1
        v2 <- Ri_in$v2
        z1 <- Ri_in$z1
        z2 <- Ri_in$z2
        Tv <- Ri_in$Tv
        
        g <- 9.8
        num <- (Ov1-Ov2) / abs(z1-z2)
        den <- ((v1-v2)/(z1-z2))^2
        return(g*num/(Tv*den))
    }

daytimeAggr <- function(xtsObj,R,aggrFunc="mean")
    {
        #xtsObj is the time series to aggregate, R is the intradaily radiation time series, 10 W/sqm is the threshold for defining day
        thresh <- 10
        xtsObj <- xtsObj[!is.na(xtsObj)]
        func <- paste('xtsObj[time(R>thresh)]')
        daytime <- eval(parse(text = func))
        return(makeDaily(daytime,aggrFunc))
    }

nighttimeAggr <- function(xtsObj,R,aggrFunc="mean")
    {
        #xtsObj is the time series to aggregate, R is the intradaily radiation time series, 10 W/sqm is the threshold for defining day
        thresh <- 10
        xtsObj <- xtsObj[!is.na(xtsObj)]
        func <- paste('xtsObj[time(R<thresh)]')
        nighttime <- eval(parse(text = func))
        index(nighttime) <- index(nighttime)+12*60*60
        return(makeDaily(nighttime,aggrFunc))
    }


sync <- function(syncTS,TS)
    {
        #syncTS is a zoo object that serves as reference for syncronyzing the time series (TS will be returned with the same time steps as syncTS)
        tmp <- xts(x=rep(NA,length(syncTS)),order.by=time(syncTS))
        tmp2 <- merge(tmp,TS)
        tmp3 <- na.approx(tmp2[,2])
        tmp4 <- merge(tmp3,tmp,join='inner')
        TS <- tmp4[,1]
        return(TS)
    }





getGolm <- function(param,t0,tf,station)
    {
        query <- paste0("SELECT SampleTimeStart,SampleValue FROM TBL_Sample WHERE LocatID='",station,"' AND ParamID='",param,"' AND SampleTimeStart>'",t0,"' AND SampleTimeEnd<'",tf,"';")
        rs <- dbSendQuery(con, query)
        d1 <- fetch(rs,n=-1)
        paramts <- xts(d1[,2],order.by=as.POSIXct(d1[,1]))
        dbClearResult(dbListResults(con)[[1]])
        return(paramts)
    }


getEVAP_PT <- function(param,t0,tf,station)
    {
        query <- paste0("SELECT SampleTimeStart,SampleValue FROM TBL_Sample WHERE LocatID='",station,"' AND ParamID='",param,"' AND SampleTimeStart>'",t0,"' AND SampleTimeEnd<'",tf,"';")
        rs <- dbSendQuery(con, query)
        d1 <- fetch(rs,n=-1)
        paramts <- xts(d1[,2],order.by=as.POSIXct(d1[,1]))
        dbClearResult(dbListResults(con)[[1]])
        return(paramts)
    }



slope_sat_vpressure <- function(T)
    {
        D <- 4098*(0.6108*exp(17.27*T/(T+237.3)))/(T+237.3)^2
        return(D) #in kPa
    }


xts2df <- function(xtsObj)
    {
        xtsObj <- xtsObj[!is.na(xtsObj)]
        year <- format(time(xtsObj),format = "%Y")
        month <- format(time(xtsObj), format = "%m")
        day <- format(time(xtsObj), format = "%d")
        hour <- format(time(xtsObj), format = "%H")
        minute <- format(time(xtsObj), format = "%M")
        second <- format(time(xtsObj), format = "%S")
        df <- data.frame(Value=coredata(xtsObj),Year=year,Month=month,Day=day,Hour=hour,Minute=minute,Second=second)
        colnames(df) <- c("Value","year","month","day","hour","minute","second")
        return(df)
    }

makeHourly <- function(xtsObj,aggrFunc)
    {
        xtsObj <- xtsObj[!is.na(xtsObj)]
        func <- paste('period.apply(xtsObj, endpoints(xtsObj, "hours"),',aggrFunc,')')
        hourly <- eval(parse(text = func))
        time(hourly) <- as.POSIXct(trunc(time(hourly),"hours"))
        return(hourly)
    }

makeDaily <- function(xtsObj,aggrFunc)
    {
        xtsObj <- xtsObj[!is.na(xtsObj)]
        func <- paste('period.apply(xtsObj, endpoints(xtsObj, "days"),',aggrFunc,')')
        daily <- eval(parse(text = func))
        time(daily) <- as.POSIXct(trunc(time(daily),"days"))
        return(daily)
    }


#rho_air <- function(p,T)
#    {
#        R <- 287.058 #J/kg.K
#        rho <- p*100/(R*(273.16+T))
#        return(rho)
#    }

putConstantToSelectedDay <- function(indexObj,xtsObj,value)
    {
        subset <- eval(parse(text=paste0(substitute(xtsObj),"['",as.Date(indexObj),"']")))
        subset <- value
        return(subset)
    }

replaceDayForNA <- function(xtsObj,xtsWorm)
    {
        #Worm will dig into Obj and replace all the values in one day by NA
        worm <- as.POSIXct(trunc(index(xtsWorm),"days"))
        for(i in seq(1,length(worm)))
            {
        #        cat(i)
                subset <- eval(parse(text=paste0("xtsObj['",trunc(worm[i],"day"),"']")))
                xtsObj[index(subset)] <- NA
            }
        return(xtsObj)
    }

replaceDailyIntoXts <- function(xtsObj,xtsWorm)
    {
        #Worm will dig into Obj and replace all the values in one day by the daily  mean
        index(xtsWorm) <- as.POSIXct(trunc(index(xtsWorm),"days"))
        subset <- as.data.frame(index(xtsWorm))
        t <- apply(subset,1,strftime)
        
        for(i in t)
            {
                                  #as.Date(eval(parse(text=paste0("index(",substitute(xtsWorm),"['",strftime(subset[i]),"'])"))))
                #print(t)
                xtsObj[i] <- xtsWorm[i]
            }
        return(xtsObj)
    }

detectLargeVaporPressureDeficit <- function(ev)
    {
#        Tdiff2 <- apply.daily(abs(T$Tdry2["T12:00/T16:00"]-T$Twet2["T12:00/T16:00"]),"mean")
 #       Tdiff1 <- apply.daily(abs(T$Tdry1["T12:00/T16:00"]-T$Twet1["T12:00/T16:00"]),"mean")
  #      error1 <- index(Tdiff1[Tdiff1 < 0.1]) #if difference between wet and dry temperature is less than 0.1 degree replace value by NA
   #     error2 <- index(Tdiff2[Tdiff2 < 0.1])
        dev <- abs(ev$ev1-ev$ev2)
        error <- index(dev[dev>4])        
        daily  <- trunc(error,"days")
        return(daily)
    }

findLargeTheta <- function(theta)
    {
        theta90 <- quantile(theta,c(0.99))
    }

if(FALSE)
    {
ShuttlWallUnder <- function(param)
    {
        #this is the Shuttleworth and Wallace equation for evapotranspiration from the understorey. check page 38 of Güntner, A. (2002). Large-Scale Hydrological Modelling in the Semi-Arid North-East of Brazil. University of Potsdam.Güntner, A. (2002). Large-Scale Hydrological Modelling in the Semi-Arid North-East of Brazil. University of Potsdam.
        
        E <- (t/lambda)*(delta*As+rho*cp*Dm/rsa)/(delta+psych*(1+rss/rsa))
        return(E)
    }
}



full_data <- function(met_data)
    {
        
        prec <- makeHourly(met_data$prec,"sum")
        R <- makeHourly(met_data$R,"mean")
        b <- makeHourly(met_data$B,"mean")

        G <- makeHourly(met_data$G,"mean")
        dQdt <- makeHourly(met_data$dQdt,"mean")
        theta <- makeHourly(met_data$theta,"mean")
        T <- makeHourly(met_data$T,"mean")
        vpd <- makeHourly(met_data$vpd,"mean")
        vpdday <- daytimeAggr(met_data$vpd,met_data$R,aggrFunc="mean")
        b <- replaceDayForNA(b,vpdday[vpdday<0.1])

        vpd <- replaceDayForNA(vpd,vpdday[vpdday<0.1])
        
        ev <- makeHourly(met_data$ev,"mean")
        ev <- replaceDayForNA(ev,vpdday[vpdday<0.1])

        u <- makeHourly(met_data$u,"mean")
        u <- replaceDayForNA(u,vpdday[vpdday<0.1])

        Td <- apply.daily(met_data$T,"mean")
        
        
        R[R < -500] <- NA 
        R[R > 1500] <- NA
        
        
        le <- LE(R,dQdt,G,b)
        h <- H(le,b)

        le[le < -1500] <- NA 
        le[le > 1500] <- NA
        h[h < -1500] <- NA 
        h[h > 1500] <- NA
              
        evap <- ET(le,rho,lambda)*3600*1000 #in mm
        evap[evap<0] <- NA

        met_full <- merge(R,T,G,theta,h,le,b,evap,vpd,ev,u)
        colnames(met_full) <- c("R","T","G","theta","h","le","b","evap","vpd","ev","u")
        return(met_full)
    }

hydro_data <- function(met_data)
    {
        
        prec <- makeHourly(met_data$prec,"sum")
        R <- makeHourly(met_data$R,"mean")
        b <- makeHourly(met_data$B,"mean")
        G <- makeHourly(met_data$G,"mean")
        dGdt <- makeHourly(met_data$dGdt,"mean")
        theta <- makeHourly(met_data$theta,"mean")
        T <- makeHourly(met_data$T,"mean")
        Td <- apply.daily(met_data$T,"mean")
                
        R[R < -500] <- NA 
        R[R > 1500] <- NA
        
        
        le <- LE(R,dGdt,G,b)
        h <- H(le,b)

        le[le < -1500] <- NA 
        le[le > 1500] <- NA
        h[h < -1500] <- NA 
        h[h > 1500] <- NA
        
        evap <- ET(le,rho,lambda)*3600*1000 #in mm
        evap[evap<0] <- NA

        hydrodata <- merge(prec,theta,evap)
        colnames(hydrodata) <- c("prec","theta","ET")
        return(hydrodata)
    }

simple_data <- function(met_data)
    {
        
        prec <- makeHourly(met_data$prec,"sum")
        R <- makeHourly(met_data$R,"mean")
        b <- makeHourly(met_data$B,"mean")
        G <- makeHourly(met_data$G,"mean")
        dGdt <- makeHourly(met_data$dGdt,"mean")
        theta <- makeHourly(met_data$theta,"mean")
        T <- makeHourly(met_data$T,"mean")
        Td <- apply.daily(met_data$T,"mean")
        
#        index(Td) <- as.POSIXct(trunc(index(Td),"days"))
#        Tdaily <- replaceDailyIntoXts(T,Td)
        
        R[R < -500] <- NA 
        R[R > 1500] <- NA
        
        
        le <- LE(R,dGdt,G,b)
        h <- H(le,b)

        le[le < -1500] <- NA 
        le[le > 1500] <- NA
        h[h < -1500] <- NA 
        h[h > 1500] <- NA
        
        evap <- ET(le,rho,lambda)*3600*1000 #in mm
        evap[evap<0] <- NA

        met_simple <- merge(R,T,evap)
        
        return(met_simple)
    }

soil_data <- function(met_data)
    {
        
        prec <- makeHourly(met_data$prec,"sum")
        G <- makeHourly(met_data$G,"mean")
        dGdt <- makeHourly(met_data$dGdt,"mean")
        theta <- makeHourly(met_data$theta,"mean")
        T <- makeHourly(met_data$T,"mean")
        Td <- apply.daily(met_data$T,"mean")
        
        index(Td) <- as.POSIXct(trunc(index(Td),"days"))
        Tdaily <- replaceDailyIntoXts(T,Td)
                        
 
        met_soil <- merge(prec,G,dGdt,theta,Tdaily)
        colnames(met_soil) <- c("Prec","G","dGdt","theta","Td")
        
        return(met_soil)
    }

