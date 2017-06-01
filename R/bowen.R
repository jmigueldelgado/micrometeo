require('RMySQL')
require('xts')
require('lattice')
require('numDeriv')
require('stats')
require('ggplot2')
require('reshape2')
require('scales')
pressure <- 101.5 #kPa
lambda <- 2260000 #latent heat of vaporization at normal pressure J/kg
rho <- 999 #kg/m3 mass density of water
 

t0 <- '2014-04-01'
tf <- '2014-09-30'


bowen <- function(Ta2,Ta1,ev2,ev1,pressure)
    {
        B <- psychr(pressure)*(Ta1-Ta2)/(ev1-ev2)
        return(B)
    }

Tvirtual <- function(T,Q) #virtual temperature
    {
        Tv <- T*(1+0.61*Q)
         #      after Arya, S. P. (2001). Introduction to Micrometeorology. International Geophysics Series (Vol. 79, p. 420). San Diego: Academic Press. doi:10.4043/13298-MS)
        return(Tv)
    }


Q <- function(ev,p) #specific humidity Mw/(Md+Mw) where Md is mass of dry air and Mw is mass of water vapour
    {
# ev is water vapour pressure and p is total pressure. to a good approximation Q=mw*ev/md*p where ma and mw are the mean molecular masses
        Q <- 0.622*ev/p # approximation 
        return(Q)

  #      after Arya, S. P. (2001). Introduction to Micrometeorology. International Geophysics Series (Vol. 79, p. 420). San Diego: Academic Press. doi:10.4043/13298-MS)) ))
    }

O <- function(T,p) #potential temperature
    {
        p0 <- 1000 #in mbar
        k <- 0.286 # k=R/cp, after Arya, S. P. (2001). Introduction to Micrometeorology. International Geophysics Series (Vol. 79, p. 420). San Diego: Academic Press. doi:10.4043/13298-MS))
        O <- T*(p0/p)^k
        return(O)
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

psychr <- function(pressure)
    {
        psy <- pressure*0.665*10^-3 #kPa/K pressure is given in kPa above
        return(psy)
    }



ET <- function(le,rho,lambda)
    {
        et <- le/(rho*lambda)
        return(et)
    }

LE <- function(R,dQdt,G,B)
    {
        le <- (R-G-dQdt)/(1+B)
        return(le)
    }

H <- function(le,B)
    {
        h <- B*le
        return(h)
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

heatCapacity <- function(xw,xorg,xsolid) #this is the volumetric heat capacity!!! organic content in soil xorg and water content xw and mineral content xsolid
    {
        Cs <- 2.45*xsolid+2.45*xorg+4.186*xw #from the encyclopedia of soil science, page 306 unit is MJ.m^-3.K^-1
        return(Cs*1e6) #unit in J/m3/K
    }

sat_vpressure <- function(T)
    {
        es <- 0.6108*exp(17.27*T/(T+237.3)) #http://www.fao.org/docrep/x0490e/x0490e07.htm in kPa
        #es <- 6.112*exp(17.62*T/(243.12+T)) #from Magnus-Formel (hier nach Sonntag 1990, siehe Literatur) http://de.wikipedia.org/wiki/S%C3%A4ttigungsdampfdruck#Literatur
        #es <- 6.11*10^(7.5*T/(237.3+T)) #http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf
        return(es)
       #in hPa
    }

slope_sat_vpressure <- function(T)
    {
        D <- 4098*(0.6108*exp(17.27*T/(T+237.3)))/(T+237.3)^2
        return(D) #in kPa
    }

vpressure <- function(Twet,Tdry)
    {
        e <- sat_vpressure(Twet)-psychr(pressure)*(Tdry-Twet)
        #in kPa
#        e <- 6.108*exp(17.27*Twet/(Twet+237.3)) #after WMO http://www.fao.org/docrep/x0490e/x0490e07.htm in hPa (originally in kPa)
        
        return(e) #hPa
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

