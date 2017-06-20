source("./bowen.R")
Dhapp <- 6.27e-7 #m^3/s

con <- dbConnect(MySQL(), user="guest",dbname="evap", host="localhost",password="")

station <- 'NSA'
R <- getEVAP_PT('NRLight2',t0,tf,station)
#G <- getEVAP_PT('HeatFlux',t0,tf,station)
Twet2 <- getEVAP_PT('Temp_Wet_TOP',t0,tf,station)
Twet1 <- getEVAP_PT('Temp_Wet_BOTTOM',t0,tf,station)
Tdry2 <- getEVAP_PT('Temp_Dry_TOP',t0,tf,station)
Tdry1 <- getEVAP_PT('Temp_Dry_BOTTOM',t0,tf,station) 
#theta <- getEVAP_PT('Soil_Moisture',t0,tf,station)
SoilTemp <- getEVAP_PT('Soil_Temp',t0,tf,station)

V1 <- getEVAP_PT('wind_BOTTOM',t0,tf,station)
V2 <- getEVAP_PT('wind_TOP',t0,tf,station)

if(FALSE)
    {
        station <- 'NSB_tower'

        SoilTemp3 <- getEVAP_PT('Soil_Temp1',t0,tf,station)
        SoilTemp2 <- getEVAP_PT('Soil_Temp2',t0,tf,station)
        SoilTemp1 <- getEVAP_PT('Soil_Temp3',t0,tf,station)
        SoilTemp <- merge.xts(SoilTemp1,SoilTemp2,SoilTemp3)
    }

station <- 'Tower'
prec <- getEVAP_PT('Precipitation',t0,tf,station)

#mysqlCloseConnection(con)

xsolid <- 0.6 #estimate based on probes!!!
xorg <- xsolid*mean(c(4.95,4.52,4.17,5.66,4.97,11.09,9.71,14.23))*0.001 #from the lab
xw <- theta # estimate, should be of class zoo or xts
Cs <- heatCapacity(xw,xorg,xsolid) #volumetric heat capacity: it should be around 2e6 J/m3/K
k <- mean(Dhapp*Cs)
G <- k*(SoilTemp1-SoilTemp3)/0.10 #0.10 is the distance between termometers. 10 *60 is seconds between measurements. Tmean is soil temperature obtaine


Cair <- 1010 # Specific Heat at constant pressure (J/kg.degreeC) "Tables of Thermal Properties of Gases", NBS Circular 564,1955


k <- mean(Dhapp*Cs)

#### this is not necessary for the energy balance in the surface layer!
#dQdt <- Cair*rho_air(pressure,Tdry1)*diff(Tdry1,lag=1)/(10*60) # the thickness of the active layer we are measuring is 1 meter (we assume the maximum height of the understorey vegetation to be 1 meter). In this case the specific heat used is not the volumetric heat capacity, but the specific heat capacity (per kg) under constant pressure.

#G <- k*(SoilTemp3-SoilTemp1)/0.10 #0.10 is the distance between termometers. 10 *60 is seconds between measurements. Tmean is soil temperature obtained from soil termometers in station NSBtower Unit is J/s, which is the same as W

dQdt <- xts(rep(0,times=length(Tdry1)),order.by=index(Tdry1))

dQdt_temp <- Cs*0.2*diff(SoilTemp,lag=1)/(10*60) #0.2 is the soil depth where stored heat is being measured. 10*60 is the time step of 10 minute turned into seconds

dQdt[index(dQdt_temp)] <- dQdt_temp


###check this again with WMO web page!!!
ev1 <- vpressure(Twet1,Tdry1) #hPa
ev2 <- vpressure(Twet2,Tdry2) #hPa
es1 <- sat_vpressure(Tdry1)
es2 <- sat_vpressure(Tdry2)
vpd1 <- es1-ev1
vpd2 <- es2-ev2


b <- bowen(Tdry2,Tdry1,ev2,ev1,pressure) #ev and pressure must have the same unit
b[b>10] <- NA #from savage 2009, J Hydrology
b[b< -10] <- NA
b[-1.25 < b & b < -0.75] <- NA #from savage 2009, J Hydrology


#detect days when the psychrometer's water container was empty ad erase them from the data. days where the average difference between dry and wet temp is lower than 0.1 degree

Tdiff <- apply.daily(abs(Tdry2["T12:00/T18:00"]-Twet2["T12:00/T18:00"]),"mean")
error <- index(Tdiff[Tdiff < 0.1]) #if difference between wet and dry temperature is less than 0.1 degree replace value by NA

T_error <- xts(x=rep(NA,times=length(error)),order.by=error)

b <- replaceDailyIntoXts(b,T_error)

met_data <- list(prec=prec,R=R,B=b,G=G,dQdt=dQdt,theta=theta,T=Tdry2,vpd=vpd2,ev=ev2,u=V2)

##see bowen.R
met_full <- full_data(met_data)
NSA_full <- met_full

                                        #met_simple <- simple_data(met_data)

Ov1 <- O(Tvirtual(273.15+Tdry1,Q(ev1,pressure)),pressure)
Ov2 <- O(Tvirtual(273.15+Tdry2,Q(ev2,pressure)),pressure)

Tv <- Tvirtual(273.15+(Tdry1+Tdry2)/2,Q((ev1+ev2)/2,pressure))

#local static stability

s <- (9.8/Tv)*(Ov1-Ov2)

#richardson number (dynamic stability)

Ri_in <- list(z1=1,z2=2,Ov1=Ov1,Ov2=Ov2,v1=V1,v2=V2,Tv=Tv)

Ri <- richardson(Ri_in)

#RiDF <- xts2df(Ri["2014-06"])

#RiDF <- RiDF[-2:-4]
#RiDF <- RiDF[-3:-4]

#upper.fence <- quantile(RiDF$Value)[4] + 4*IQR(RiDF$Value)
#lower.fence <- quantile(RiDF$Value)[2] - 4*IQR(RiDF$Value)

#upper.fence <- 150
#lower.fence <- -250

#ggplot(RiDF,aes(hour,Value)) + geom_boxplot() + coord_cartesian(ylim=c(lower.fence, upper.fence))
#ggsave(file="Ri_jun_NSA.pdf")

#colnames(met_full) <- c("R","T","G","H","LE","B","ET")
#colnames(met_simple) <- c("R","T","ET")

