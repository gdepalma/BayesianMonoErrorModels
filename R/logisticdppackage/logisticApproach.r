
library(mixtools)
library(reshape2)
library(DPpackage)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggthemes)


source('GenData.R')
source('logisticFunctions.R')
source('initialize_parms_log.R')
dyn.load('splineFunctionsC.so')

### Generate Data
xsig=.707; ysig=2.121
temp=genData(type=6,nobs=500)
xobs=temp$xobs; yobs=temp$yobs; xcens=temp$xcens; ycens=temp$ycens; 
trueFit=temp$true; xtrue=temp$xtrue; densTrue=temp$densTrue; xgrid=temp$xgrid
nobs=length(xobs)

parms=initialize_parms_log(xobs,yobs,xcens,ycens)
coefs=parms$coefs
xtrue=parms$xtrue
ytrue=getYtrue(coefs,xtrue)

### Update Dens
parms=updateDens(xtrue,state=0,first=TRUE,niter=1000,xgrid)
# mu=parms$mu; sigma=parms$sigma; p=parms$p; state=parms$state;
# dens=parms$densEst; k=length(p); groups=parms$groups
k=parms$k; densEst=parms$densEst; state=parms$state


### Set up variables
burnin=20000
numIter=30000
coefMat=matrix(nrow=numIter,ncol=length(coefs))
fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
MICDens=matrix(nrow=numIter-burnin,ncol=length(xgrid))

#### smooth coefficient parameters
a1=.5
a2=.03
a3=.02
a4=.03
  
# sigmaM=rep(NA,numIter)
# sigmaD=rep(NA,numIter)
acceptCoef=rep(NA,numIter)
llk=rep(NA,numIter-burnin)
popllk=rep(NA,numIter-burnin)
prior=rep(NA,numIter-burnin)
post=rep(NA,numIter-burnin)
coefMatSave=matrix(nrow=numIter-burnin,ncol=length(coefs))
kSave=rep(NA,numIter-burnin)

coefMonitor=matrix(nrow=numIter,ncol=4)


tmpTime=Sys.time()
begin=Sys.time()

for(iter in 1:numIter){
  
  if(iter %% 500==0){
    cat('Iter: ',iter,'\n')
    print(Sys.time()-tmpTime)
    tmpTime=Sys.time()
  }
  
  ## Update Density
  if(iter %% 500==0){
    parms=updateDens(xtrue,state,first=FALSE,niter=500,xgrid)
    # mu=parms$mu; sigma=parms$sigma; p=parms$p; state=parms$state; 
    # dens=parms$densEst; k=length(p); groups=parms$groups
    k=parms$k; densEst=parms$densEst; state=parms$state
  }
  
  ### Update xtrue (and corresponding ytrue)
  # if(iter %% 25==0 | iter==1){
    parms=updateXtrue(xtrue,ytrue,coefs,xcens,ycens,xsig,ysig,xobs,yobs,
         nobs,mu,sigma,k,p,groups,densEst)
    xtrue=parms$xtrue; ytrue=parms$ytrue
  # }
    
  
  ####update Coefs
  parms=updateLogCoef_2(iter,coefs,xobs,yobs,xtrue,xcens,ycens,ytrue,xsig,ysig,a1,a2,a3,a4)
  acceptCoef[iter]=parms$accept; coefs=parms$coefs; ytrue=parms$ytrue; coefMonitor[iter,]=parms$accepts
  
  ### Update Logistic Jump Proposals
  if(iter %% 100==0){
      print(c(a1,a2,a3,a4))
    
      if(mean(coefMonitor[(iter-99):iter,1]) > .4) a1=a1+.01
      if(mean(coefMonitor[(iter-99):iter,1]) < .2) a1=a1-.01
      
      if(mean(coefMonitor[(iter-99):iter,2]) > .4) a2=a2+.01
      if(mean(coefMonitor[(iter-99):iter,2]) < .2) a2=a2-.01
      if(a2<=0) a2=.01
      
      if(mean(coefMonitor[(iter-99):iter,3]) > .4) a3=a3+.005
      if(mean(coefMonitor[(iter-99):iter,3]) < .2) a3=a3-.005
      if(a3<=0) a3=.005
      
      if(mean(coefMonitor[(iter-99):iter,4]) > .4) a4=a4+.005
      if(mean(coefMonitor[(iter-99):iter,4]) < .2) a4=a4-.005
      if(a4<=0) a4=.005
  }
  
  
  ### Update Measurement Errors
#   parms=updateMeasErrorM(xsig,ysig,xobs,yobs,xtrue1,ytrue1,xcens,ycens)
#   xsig=parms$xsig
#   parms=updateMeasErrorD(xsig,ysig,xobs,yobs,xtrue1,ytrue1,xcens,ycens)
#   ysig=parms$ysig
    
  
  ###save
  coefMat[iter,]=coefs
  # sigmaM[iter]=xsig
  # sigmaD[iter]=ysig
  if(iter>burnin){    
    llk[iter-burnin]=sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))
    popllk[iter-burnin]=sum(getPOPLLK2(xtrue,densEst))
    prior[iter-burnin]=sum(getPriorLLK(coefs))
    post[iter-burnin]=llk[iter-burnin]+popllk[iter-burnin]+popllk[iter-burnin]+prior[iter-burnin]
    fitMat[iter-burnin,]=getYtrue(coefs,xgrid)
    coefMatSave[iter-burnin,]=coefs
    MICDens[iter-burnin,]=densEst
    kSave[iter-burnin]=k
  }
  
}
print('Finished') 
print(Sys.time()-begin)



### Diagnostics

