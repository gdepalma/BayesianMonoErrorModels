
library(mixtools)
library(ggplot2)
library(reshape2)
library(DPpackage)
library(ggplot2)

# setwd("C:/Users/glend_000/Dropbox/Research Projects/LogisticDPPackage")
setwd("C:/Users/glend_000/Desktop/040617/LogisticDPPackage")

source(paste(getwd(),'/GenData.R',sep=''))
source(paste(getwd(),'/logisticFunctions3Parm.R',sep=''))
source(paste(getwd(),'/logisticApproachFunction3Parm.R',sep=''))
# dyn.load('splineFunctionsC2.dll')
# dyn.load('splineFunctionsC.so')


nsim=1
saveMICDensEst=matrix(nrow=nsim,ncol=1000)
saveFitEst=matrix(nrow=nsim,ncol=1000)


for(i in 1:nsim){
  
  if(i %% 2==0){
    cat('Iter: ',i,'\n')
  }
  
  parms=runLogistic3(genDataType=5)
  saveMICDensEst[i,]=parms$MICDens
  saveFitEst[i,]=parms$fit
}


tmp=list(saveMICDensEst=saveMICDensEst,saveFitEst=saveFitEst)
save(tmp,file=paste(getwd(),'/Runs/LogisticGen5_three.Rdata',sep=''))
