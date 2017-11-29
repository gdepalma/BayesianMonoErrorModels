
library(mixtools)
library(reshape2)
library(DPpackage)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggthemes)


source('GenData.R')
source('splineFunctions.R')
source('initialize_parms_spline.R')
dyn.load('splineFunctionsC.so')

### Generate Data
xsig=.707; ysig=2.121
temp=genData(type=1,nobs=500)
xobs=temp$xobs; yobs=temp$yobs; xcens=temp$xcens; ycens=temp$ycens; 
trueFit=temp$true; xtrue=temp$xtrue; densTrue=temp$densTrue; xgrid=temp$xgrid
nobs=length(xobs)

parms=initialize_parms_spline(xobs,yobs,xcens,ycens,xgrid)
coefs=parms$coefs
xtrue=parms$xtrue
ytrue=parms$ytrue
bases=parms$bases
knotseq=parms$knotseq
lowept=parms$lowept
upperept=parms$upperept
designMatrixGrid=parms$designMatrixGrid

### Update Dens
parms=updateDens(xtrue,state=0,first=TRUE,niter=1000,xgrid)
k=parms$k; densEst=parms$densEst; state=parms$state


### Set up variables
burnin=5000
numIter=10000
coefMat=matrix(nrow=numIter,ncol=length(coefs))
fitMat=matrix(nrow=numIter-burnin,ncol=length(xgrid))
MICDens=matrix(nrow=numIter-burnin,ncol=length(xgrid))


llk=rep(NA,numIter-burnin)
popllk=rep(NA,numIter-burnin)
prior=rep(NA,numIter-burnin)
post=rep(NA,numIter-burnin)
coefMatSave=matrix(nrow=numIter-burnin,ncol=length(coefs))
acceptCoef=rep(NA,numIter)
smoothAccept=rep(NA,numIter)
kSave=rep(NA,numIter-burnin)
smoothSplineSave=rep(NA,numIter)
smoothSpline = .5

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
  parms=updateSplineXtrue(xtrue,ytrue,knotseq,bases,lowept,upperept,coefs,xcens,
                          ycens,xsig,ysig,xobs,yobs,nobs,densEst)
  xtrue=parms$xtrue; ytrue=as.numeric(parms$ytrue)

  
  ### Update Coefs
  parms=updateSplineCoefs(smoothSpline,coefs,iter,bases,xtrue,
                          ytrue,knotseq,xobs,yobs,xcens,ycens,xsig,ysig,coefMat)
  acceptCoef[iter]=parms$acceptCoef; coefs=parms$coefs; ytrue=as.numeric(parms$ytrue)
  
  ### update smoothness parameter
  parms=updateSmoothParm(smoothSpline,coefs)
  smoothSpline=parms$smooth; smoothAccept[iter]=parms$accept
  
  ###save
  coefMat[iter,]=log(coefs)

  if(iter>burnin){    
    llk[iter-burnin]=sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))
    popllk[iter-burnin]=sum(getPOPLLK(xtrue,densEst))
    prior[iter-burnin]=sum(getPRIORLLK(coefs,smoothSpline))
    post[iter-burnin]=llk[iter-burnin]+popllk[iter-burnin]+popllk[iter-burnin]+prior[iter-burnin]
    fitMat[iter-burnin,]=coefs%*%t(designMatrixGrid)
    coefMatSave[iter-burnin,]=coefs
    MICDens[iter-burnin,]=densEst
    kSave[iter-burnin]=k
  }
  
}
print('Finished') 
print(Sys.time()-begin)



### Diagnostics

