require(micrometeo)
require(plyr)
require(data.table)
require("xts")
require("dplyr")

f1 <- list("COMBILOG29","COMBILOG28","COMBILOG26","COMBILOG25","COMBILOG23","COMBILOG22","COMBILOG21","COMBILOG20","COMBILOG18")

readCOMBILOG_HS <- function(f)
{
    
    x <- read.table(paste0("~/GEOECOLOGY/user/delgado/EVAP_PT/HS/",f,".LOG"),header=FALSE,skip=8,sep="\t")
    h <- readLines(paste0("~/GEOECOLOGY/user/delgado/EVAP_PT/HS/",f,".LOG"),n=6)[6]
    h <- unlist(strsplit(h,split="\t"))
    h <- make.names(h)
    colnames(x) <- h
    x <- x[!is.na(names(x))]
    tt <- as.POSIXct(x$TimeDate,format="%m/%d/%Y %I:%M:%S %p",tz="Portugal")
    xx <- cbind(data.frame(time=tt),x[,c(-1,-2)])
    
    x <- read.table(paste0("~/GEOECOLOGY/user/delgado/EVAP_PT/HS/",f,"_S.LOG"),header=FALSE,skip=8,sep="\t")
    h <- readLines(paste0("~/GEOECOLOGY/user/delgado/EVAP_PT/HS/",f,"_S.LOG"),n=6)[6]
    h <- unlist(strsplit(h,split="\t"))
    h <- make.names(h)
    colnames(x) <- h
    x <- x[!is.na(names(x))]
    tt <- as.POSIXct(x$TimeDate,format="%m/%d/%Y %I:%M:%S %p",tz="Portugal")
    xx_S <- cbind(data.frame(time=tt),x[,c(-1,-2)])
    
    dat <- inner_join(xx,xx_S,by="time")
    
    
    return(dat)
}

ll <- lapply(f1,readCOMBILOG_HS)

df <- rbindlist(ll)
df <- df[with(df,order(time)),]
setkey(df,"time")
df2 <- unique(df,by="time")

df2 <- rename(df2,Soil.Temp_tmp=Soil.Temp3) %>% rename(Soil.Temp3=Soil.Temp2,Soil.Temp2=Soil.Temp_tmp)


### for soil heat flux
xsolid <- 0.6 #estimate based on probes!!!
xorg <- xsolid*mean(c(4.95,4.52,4.17,5.66,4.97,11.09,9.71,14.23))*0.001 #from the lab
xw <-  0.1  #theta # estimate, should be of class zoo or xts
Cs <- vol_heat_capacity_soil(xw,xorg,xsolid) #volumetric heat capacity: it should be around 2e6 J/m3/K
k <- mean(Dhapp*Cs)
G <- k*(df2$Soil.Temp1-df2$Soil.Temp3)/0.10 #0.10 is the distance between termometers. 10 *60 is seconds between measurements. Tmean is soil temperature obtaine

k <- mean(thermal_diff()*Cs) # thermal conductivity

dQdt <- c(0,Cs*0.2*diff(df2$Soil.Temp3,lag=1)/(10*60)) #0.2 is the soil depth where stored heat is being measured. 10*60 is the time step of 10 minute turned into seconds


ev1 <- wetbulb2vpressure(df2$Temp.Wet.BOTTOM,df2$Temp.Dry.BOTTOM,p_coruche()) #hPa
ev2 <- wetbulb2vpressure(df2$Temp.Wet.TOP,df2$Temp.Dry.TOP,p_coruche()) #hPa

q1 <- Q(ev1,p_coruche())
q2 <- Q(ev2,p_coruche())

rh1 <- specific_hum2rh(q1,df2$Temp.Dry.BOTTOM,p_coruche())
rh2 <- specific_hum2rh(q2,df2$Temp.Dry.TOP,p_coruche())



b <- bowen(df2$Temp.Dry.TOP,df2$Temp.Dry.BOTTOM,ev2,ev1,p_coruche()) #ev and pressure must have the same unit
b[b>10] <- NA #from savage 2009, J Hydrology
b[b< -10] <- NA
b[-1.25 < b & b < -0.75] <- NA #from savage 2009, J Hydrology
