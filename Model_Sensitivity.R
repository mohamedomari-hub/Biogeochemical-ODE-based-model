rm(list=ls());
library(deSolve);
library(FME);
library(reshape2)
#library(zoo)
library(xlsx)

source('~/Desktop/Docs/Academic/MODEL_FINAL_VERSION_ZALF/ModelFinal_Control.R')

Y <- c(Algae=0.8, DWeed = 5, DO=0.05,
       POC=0.1 , DOC=19.7, DIC=74., TC=95,
       ON=0.7, NH4_N=0.01, NO3_N=0.35, TN= 3.,
       OP=0.15, PO4=0.5, TP= 0.7, TPOC=0.1, 
       TON=2., TOP=0.23 ) 



sensvar= names(Y)

senspar=names(parms)

#senspar= c("ArrhBM","ArrhGr")


#names(parms)#c("ArrhGr", "PCPhy", "NCPhy")
selec_par = 10#length(senspar) # should be less than length of senspar
num_par <- length(senspar)
num_var <- length(sensvar)

#n <- parms[c("ArrGrowth", "Arrhres","Arrhmor","arrhdoc", 
#             "ArrhresVgs","arrhdenit","Arrhnitr","ArrhdecC","arrhdon","ArrhdecP",
#             "Arrhpon", "Arrhpopset", "ArrhresVgr","ArrhdecN",
#             "TemeffDOC","Arrhpoc", "Arrhponset","TemeffPO","arrhdop","Arrhpop", "TemeffDIC",
#             "poro","decN","PH","ks3","nitSed", "mindon","PerPO4","decP","assi","DOsatsod","KBiomass","Gmax")]


print(system.time(SnsCarb <- sensFun(func = Carbon, parms=parms, sensvar= sensvar, senspar= senspar,
                                     varscale = NULL, parscale=NULL , tiny = 1e-1)))  

#SnsCarb <- sensFun(func = Carbon, parms=params, senspar= c("Arrhmor", "ArrGrowth"), 
#                                                            varscale = 1, tiny = 1e-1)  

for (i in 1:length(Y)){

plot(SnsCarb, which=names(Y)[i] ,xlab="time", lwd = 2,legpos="right")

}

#Plotting the sensitivity function:
if(FALSE){
resize.win <- function(Width=6, Height=6)
{
  # works for windows
  # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}
dev.off();
resize.win(5,5)
par(mar=c(6,5,2,2))
layout(matrix(c(1,1,1,1)), widths=c(3,3), heights=c(.75,2.25))
plot(SnsCarb,which=names(Y)[13:17],col = colors()[c(552,254,26,61,90,455,151,451,70,574,142,153,260,7,87,258,19,382)]
       #random(rainbow(18, s = 1, v = 1, start = 0, end = max(1, 18 - 1)/18, alpha = 1))
       #sample(colors(distinct=T), 18)
       ,legpos=NA)
plot.new()
par(xpd=TRUE)
legend("center",pch=15, bg="white",col=colors()[c(552,254,26,61,90,455,151,451,70,574,142,153,260,7,87,258,19,382)],
       legend=senspar)
#c(1,2,6,7,8,10,11,17,19,35,83)
}


Summary <- summary(SnsCarb, var = TRUE)

#write.xlsx2(Summary, file = "C:/Users/omari/Desktop/MODEL/Monte_Output/Result5/MonteCarloSensitivity3.xls", sheetName= "Sheet1", 
            col.names=TRUE, row.names=TRUE, append=F, showNA=TRUE)


#print(summary(SnsCarb, var = TRUE))

################################################### Ranking ##############################################################

#Summary <- read.xlsx("C:/Users/omari/Desktop/MODEL/Monte_Output/Result5/MonteCarloSensitivity3.xls",sheetIndex=1)


# Sorting the varibale from the summary function of senstivity analysis.
#-----------------------------------------------------------------------


d <- matrix(nrow = num_par, ncol = num_var)

for (k in 1:num_var) 
  {
  for (i in 0:(num_par-1))
  {
    d[i+1,k] <- Summary[(i*num_var) +k,5]
  }
  }

# ordring the variables
#--------------------------

m <- matrix(nrow = num_par, ncol = num_var)

for (i in 1:num_var)
  
  {
  m[,i] <- sort(d[,i], decreasing = T)
  }



#write.xlsx2(m, file = "C:/Users/omari/Dropbox/R/Monte_Output/Result5/MonteCarloParameterRank3_value.xls", sheetName= "Sheet1", 
#            col.names=TRUE, row.names=TRUE, append=F, showNA=TRUE)

#write.xlsx2(m[1:selec_par,], file = "C:/Users/omari/Dropbox/R/Monte_Output/Result5/MonteCarloParameterSelected_value.xls", sheetName= "Sheet1", 
#            col.names=TRUE, row.names=TRUE, append=F, showNA=TRUE)


#order.value <- as.data.frame(read.xlsx("C:/Users/omari/Dropbox/R/Monte_Output/Result5/MonteCarloParameterRank3_value.xls",sheetIndex=1))
#order.name <- as.data.frame(read.xlsx("C:/Users/omari/Dropbox/R/Monte_Output/Result5/MonteCarloParameterRank3.xls",sheetIndex=1))


parametrs <- data.frame()


# Ranking the parametrs
#----------------------


for (h in 1:num_var)
{
  
  #m is sens_order and Summary is summary
  for(i in 1:(nrow(m))) { 
    a = round(m[i,h],digits=12)
    for (j in 1:nrow(Summary)){
      b =  round(Summary[j,5], digits=12)
      if(a == b)  parametrs[i,h] <- Summary[j,1]
    }
  }
}

# writing the data.frame in order!
#--------------------------------


#write.xlsx2(parametrs, file = "C:/Users/omari/Dropbox/R/Monte_Output/Result5/MonteCarloParameterRank3.xls", sheetName= "Sheet1", 
#            col.names=TRUE, row.names=TRUE, append=F, showNA=TRUE)

#10 most influencing para Selected
#---------------------------------

empty <- c(rep(NA,num_var))

name_col <- Summary$var[1:num_var]

x <- as.matrix(rbind(levels(name_col),empty, round(m, digits=2)))

y <- as.matrix(rbind( x[1:(selec_par+2),], empty, parametrs[1:selec_par,]))

names(y) <- NULL

#suppressWarnings(write.xlsx2(y, file = "C:/Users/omari/Desktop/MODEL/Monte_Output/Result5/MonteCarloParameterSelected.xls", sheetName= "Sheet1", 
#            col.names=TRUE, row.names=TRUE, append=F, showNA=TRUE))

# parameters to estimate
#=======================

table<- read.xlsx("C:/Users/omari/Desktop/MODEL/Monte_Output/Result5/MonteCarloParameterSelected.xls",sheetIndex=1)
#table <- table[,-1]
sele.par <- data.frame()

for(j in 1:4)
  
{

for(i in 1:selec_par)
  
  {
  if (j<10) sele.par[i,j] <- substr(as.character(table[i + selec_par,j+1]), 1, nchar(as.character(table[i+selec_par,j+1]))- 1)
  
  else sele.par[i,j] <- substr(as.character(table[i + selec_par,j+1]), 1, nchar(as.character(table[i+selec_par,j+1]))- 2)
  
  }

}

para <- unique(unlist(sele.par))

#write(para, file = "C:/Users/omari/Desktop/MODEL/Monte_Output/Result5/para.txt")

#==============================================================================================================
#BindedData
if(FALSE){

fileName <- c(parametrs,m)  #parametrs dataframe has to have a column's name of parameters' names to be defined

df.list <- vector("list", length(fileName))

df.list2 <- vector("list", num_var)

for(i in 1:2){df.list[[i]] <- fileName[i] }


for(j in 2:(num_var + 1)){
  
  data.bind <- matrix( nrow=nrow(df.list[[1]]))
  
  for (k in 1:13)    
  {   
    
    data.bind <- cbind(data.bind, df.list[[k]][,j]) 
  }
  df.list2[[j-1]] <- data.bind        
}



#write.xlsx2(data.bind1, file = "C:/Users/omari/Desktop/MODEL/Monte_Output/Result5/MonteCarlo_ISOMAP.xls", sheetName= "Sheet1", 
            col.names=TRUE, row.names=TRUE, append=F, showNA=TRUE)
}


############################################################################################################
