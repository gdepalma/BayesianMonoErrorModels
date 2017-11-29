
library(mixtools)
library(ggplot2)
library(reshape2)
library(DPpackage)
library(ggplot2)


source(paste(getwd(),'/GenData.R',sep=''))
source(paste(getwd(),'/splineFunctions.R',sep=''))
source(paste(getwd(),'/splineApproachFunction.R',sep=''))
dyn.load('splineFunctionsC.so')


nsim=200
saveMICDensEst=matrix(nrow=nsim,ncol=1000)
saveFitEst=matrix(nrow=nsim,ncol=1000)


for(i in 1:nsim){
  
  if(i %% 2==0){
    cat('Iter: ',i,'\n')
  }
  
  parms=runSpline(genDataType=1)
  saveMICDensEst[i,]=parms$MICDens
  saveFitEst[i,]=parms$fit
}


tmp=list(saveMICDensEst=saveMICDensEst,saveFitEst=saveFitEst)
save(tmp,file=paste(getwd(),'/Runs/LogisticGen2_four.Rdata',sep=''))
