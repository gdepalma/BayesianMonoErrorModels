
library(tidyverse)
library(DPpackage)
library(mvtnorm)
library(reshape2)
library(ggthemes)
library(gridExtra)
library(grid)


source('logisticFunctions.R')
source('Initialize_parms_logistic.R')
a1=read_csv(file='DATA/ERTEB.csv')

xcens=rep(0,nrow(a1))
ycens=rep(0,nrow(a1))
xobs=rep(NA,nrow(a1))
yobs=rep(NA,nrow(a1))

#check for censoring

for(i in 1:nrow(a1)){
  if(grepl('>',a1[i,1])==TRUE){
    xobs[i]=substr(a1[i,1],2,nchar(a1[i,1]))
    xcens[i]=1
  }else if(grepl('<',a1[i,1])==TRUE){
    xobs[i]=substr(a1[i,1],2,nchar(a1[i,1]))
    xcens[i]=-1
  }else
    xobs[i]=a1[i,1]
  
  if(grepl('>',a1[i,2])==TRUE){
    yobs[i]=substr(a1[i,2],2,nchar(a1[i,2]))
    ycens[i]=1
  }else if(grepl('<',a1[i,2])==TRUE){
    yobs[i]=substr(a1[i,2],2,nchar(a1[i,2]))
    ycens[i]=-1
  }else
    yobs[i]=a1[i,2]
}

# xobs=as.numeric(xobs)
xobs=log(as.numeric(xobs),2)
yobs=as.numeric(yobs)

### Initialize
nobs=length(xobs)
xsig=.707; ysig=2.121
xgrid=seq(min(xobs)-3,max(xobs)+3,length=1000)
parms=initialize_parms_logistic(xobs,yobs,xcens)
xtrue=parms$xtrue
coefs=parms$coefs
ytrue=getylogtrue1(coefs,xtrue)

### Update Dens
parms=updateDens(xtrue,state=0,first=TRUE,niter=1000,xgrid)
mu=parms$mu; sigma=parms$sigma; p=parms$p; state=parms$state;
dens=parms$densEst; k=length(p); groups=parms$groups

### Set up variables
burnin=6000
numIter=10000
coefMat=matrix(nrow=numIter,ncol=length(coefs))
fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
MICDens=matrix(nrow=numIter-burnin,ncol=length(xgrid))
acceptCoef=rep(NA,numIter)
post=rep(NA,numIter)
xtrue_sav=matrix(nrow=numIter-burnin,ncol=length(xtrue))


tmpTime=Sys.time()
begin=Sys.time()

for(iter in 1:numIter){
  
  if(iter %% 500==0){
    cat('Iter: ',iter,'\n')
    print(Sys.time()-tmpTime)
    tmpTime=Sys.time()
  }
  
  ## Update Density
  if(iter %% 1000==0){
    parms=updateDens(xtrue,state,first=FALSE,niter=500,xgrid)
    mu=parms$mu; sigma=parms$sigma; p=parms$p; state=parms$state; 
    dens=parms$densEst; k=length(p); groups=parms$groups
  }
  
  ### Update xtrue (and corresponding ytrue)
  if(iter %% 25==0 | iter==1){
    parms=updateXtrueLog(xtrue,ytrue,coefs,xcens,ycens,xsig,ysig,xobs,yobs,
         nobs,mu,sigma,k,p,groups)
    xtrue=parms$xtrue; ytrue=as.numeric(parms$ytrue)
  }
  
  ### xupdate Coefs
  parms=updateLogCoef(iter,coefs,xobs,yobs,xtrue,xcens,ycens,ytrue,xsig,ysig,coefMat)
  acceptCoef[iter]=parms$accept; coefs=parms$coefs; ytrue=parms$ytrue
  
  
  ###save
  coefMat[iter,]=coefs
  post[iter]=sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))
  if(iter>burnin){    
    fit=getylogtrue1(coefs,xgrid)
    fitMat[iter-burnin,]=fit
    MICDens[iter-burnin,]=dens
    xtrue_sav[iter-burnin,]=xtrue
  }
  
}
print('Finished') 
print(Sys.time()-begin)


### Thin
MICDens=MICDens[seq(1,nrow(MICDens),by=2),]
fitMat=fitMat[seq(1,nrow(fitMat),by=2),]

coefMat=coefMat[(burnin+1):numIter,]
acceptCoef=acceptCoef[(burnin+1):numIter]






