rm(list=ls())
library(deSolve)
library(zoo)                                                  # Kettl hole volume: 5747 m3  (700 m2)
library(xlsx)
library(Hmisc)
library(scatterplot3d)


# Running ODE models implemented in C
# ------------------------------------

setwd("~/Desktop/Docs/Academic/MODEL_FINAL_VERSION_ZALF")

#dyn.unload("C:/Users/omari/Desktop/MODEL/R_C_Models/ModelFinal.dll")

system("R CMD SHLIB ~/Desktop/Docs/Academic/MODEL_FINAL_VERSION_ZALF/ModelFinal.c")

# Creating DLL file 
#--------------------
dyn.load("~/Desktop/Docs/Academic/MODEL_FINAL_VERSION_ZALF/ModelFinal.so")

#Model initialisation
#--------------------
source('~/Desktop/Docs/Academic/MODEL_FINAL_VERSION_ZALF/ModelFinal_para.R')

#Compiling the model by calling it from C
#-----------------------------------------

# Time interval
#---------------
numpoint = 14000 #736
steps = 1
times <- seq(1,numpoint, by = steps )

Carbon <- function ( parms, times) {
  
  times <- seq(1,numpoint, by = steps )
  
  
  
  Y <- c(Algae=0.8, DWeed = 5, DO=0.05,
         POC=0.1 , DOC=19.7, DIC=74., TC=95,
         
         ON=0.7, NH4_N=0.01, NO3_N=0.35, TN= 3.,
         
         OP=0.15, PO4=0.5, TP= 0.7, TPOC=0.1, TON=2., TOP=0.23 ) 
  
  
  
  # Forcing function:
  # ----------------
  
  data <- read.table('~/Desktop/Docs/Academic/MODEL_FINAL_VERSION_ZALF/data_forcing3.csv', header=TRUE, dec = ",", sep=";")
  
  
  dat <- as.data.frame(list(Times = data[1:numpoint,1], temp = data[1:numpoint,4], H = data[1:numpoint,5],
                            I0 = data[1:numpoint,3], Ph = data[1:numpoint,6] , Wind = data[1:numpoint,7]))#, Ca = data[1:numpoint,13],Fe = data[1:numpoint,14]))
  
  Depth = c( 1.78,rep(NA,1545), 1.84,rep(NA,1545), 1.90,rep(NA,1545), 1.55,rep(NA,1545), 1.00, rep(NA,1545),
             0.80,rep(NA,1545), 0.65,rep(NA,1545), 0.60,rep(NA,1545), 0.70,rep(NA,1545), 0.82,rep(NA,84), 0.90)
  
  volume = c( 1024,rep(NA,1545),1098,rep(NA,1545), 1176,rep(NA,1545), 764,rep(NA,1545),  302, rep(NA,1545),
              178,rep(NA,1545), 100,rep(NA,1545), 81,rep(NA,1545), 122,rep(NA,1545), 183,rep(NA,84), 236)
  
  
  ph = c(6.72,rep(NA,1545),7.62,rep(NA,1545),6.56,rep(NA,1545),6.56,rep(NA,1545),6.36,rep(NA,1545),
         6.38,rep(NA,1545),6.57,rep(NA,1545),6.79,rep(NA,1545),7.08,rep(NA,1545),7.54,rep(NA,84),7.07)
  
  Mn = c( 1.96, rep(NA,1545), 0.05, rep(NA,1545), 0.80, rep(NA,1545), 0.35, rep(NA,1545),1.08, rep(NA,1545),
          1.07, rep(NA,1545), 0.96, rep(NA,1545),0.64, rep(NA,1545), 0.33, rep(NA,1545), 1.00, rep(NA,84), 1.8)
  
  Cal = c( 107.3, rep(NA,1545), 99.3, rep(NA,1545), 91, rep(NA,1545), 95, rep(NA,1545), 84, rep(NA,1545),
           81, rep(NA,1545), 71.33, rep(NA,1545),73, rep(NA,1545), 72, rep(NA,1545), 74.3, rep(NA,84), 105)
  
  Fer = c( 0.08, rep(NA,1545), 0.03, rep(NA,1545), 0.17, rep(NA,1545), 0.08, rep(NA,1545),0.09, rep(NA,1545),
           0.11, rep(NA,1545), 0.09, rep(NA,1545),0.08, rep(NA,1545), 0.07, rep(NA,1545), 0.14, rep(NA,84), 0.23 )
  
  
  Fer = na.approx(Fer,1:numpoint)
  Mn = na.approx(Mn,1:numpoint)
  Cal = na.approx(Cal,1:numpoint)
  dat$H = na.approx(Depth,1:numpoint) 
  dat$Ph = 8.7#na.approx(ph,1:numpoint)
  Volume = na.approx(volume,1:numpoint)*1e3  # multip by 1e3 to convert from m3 to L,  1m3 = 1e3 L
  
  Sigimp1   <- approxfun(dat$Times, dat$temp, rule = 2)
  Sigimp2   <- approxfun(dat$Times, dat$H, rule = 2)
  Sigimp3   <- approxfun(dat$Times, dat$I0, rule = 2)
  Sigimp4   <- approxfun(dat$Times, dat$Ph, rule = 2)
  Sigimp5   <- approxfun(dat$Times, Cal, rule = 2)
  Sigimp6   <- approxfun(dat$Times, Fer, rule = 2)   
  Sigimp7   <- approxfun(dat$Times, Mn, rule = 2)
  Sigimp8   <- approxfun(dat$Times, dat$Wind, rule = 2) 
  Sigimp9   <- approxfun(dat$Times, Volume, rule = 2)
  
  forcings  <- list(cbind(times, Sigimp1(times)),cbind(times, Sigimp2(times)),cbind(times, Sigimp3(times)),
                    cbind(times, Sigimp4(times)),cbind(times, Sigimp5(times)),cbind(times, Sigimp6(times))
                    ,cbind(times, Sigimp7(times)),cbind(times, Sigimp8(times))
                    ,cbind(times, Sigimp9(times)))  
  
  
  #Solving the ODEs:
  #-----------------
  #lsodes , hmax=1e-1
  
  out <- lsodes(Y, times, func = "derivs", parms = parms, rtol =1e-4, atol = 1e-4, maxsteps = 2500000,
                dllname = "ModelFinal", initforc = "forcc", forcings=forcings, initfunc = "initmod")
}

#Results and Plots:
#------------------

data <- read.table('~/Desktop/Docs/Academic/MODEL_FINAL_VERSION_ZALF/data_forcing3.csv', header=TRUE, dec=",", sep=";")

print(system.time(output <- Carbon( parms, times)))


#++++++++++++++++++++++++++++++++++++++++++++++plotting the results+++++++++++++++++++++++++++++++

TimeSampl <- strptime(data[1:10,9], "%m.%d.%Y %H:%M")

lab <- c(rep("mg/l",1),rep("gDW/m2",1),rep("mg/l",  12))

Col <- c(rep("red",1),rep("red",1),rep("red",12))

nameoutput <- names(output)

par(mfrow = c(4, 4),ask=F,mar=c(2,5,5,1))

for(i in 2:ncol(output)) 
{ 
  
  xlabel <- strptime(data[seq(1,numpoint, by=steps),2], "%m.%d.%Y %H:%M")
  
  plot(xlabel,output[,i], xaxt = "n" , ylim=range(output[,i],na.rm=T),
       main=nameoutput[i],ylab=lab[i-1],t="l",col="blue",cex.axis=1., cex.main=1.3, cex.lab=1.5)
  
  par(new=TRUE)

  #grid(nx=NA, ny=NULL)
  
  abline(v=axis.POSIXct(1, at=pretty(xlabel,n=15), format = "%b"),col = "lightgray", lty = "dotted", lwd = par("lwd"))
}



  