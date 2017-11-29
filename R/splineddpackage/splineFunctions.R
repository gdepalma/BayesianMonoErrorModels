
getLLK=function(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig)
{
  nobs=length(xobs)
  onevec <- c(rep(1,nobs))
  xcensl=rep(0,nobs)
  xcensl[xcens==-1] = 1
  xcensu=rep(0,nobs)
  xcensu[xcens==1] = 1
  
  pmic <- pnorm(xobs,xtrue,xsig)*(1-xcensu) + onevec*xcensu - 
    pnorm(xobs-1,xtrue,xsig)*(1-xcensl)
  pmic[pmic<.001]=.0001
  
  ycensl=rep(0,nobs)
  ycensl[ycens==-1] = 1
  ycensu=rep(0,nobs)
  ycensu[ycens==1] = 1
  
  pdia <- pnorm(yobs+0.5,ytrue,ysig)*(1-ycensu) + onevec*ycensu -
    pnorm(yobs-0.5,ytrue,ysig)*(1-ycensl)
  pdia[pdia<.001]=.0001
  
  llike <- log(pmic) + log(pdia)
  
  return(llike)
}


getPOPLLK=function(x,dens)
{
  ###interpolate  
  pxden=log(approx(xgrid,dens,xout=x,rule=2)$y)
  return(pxden)
}

getPRIORLLK=function(coefs,smoothSpline){
  
  priorPcoef=dunif(smoothSpline,0,3,log=T)
  
  # lpriorcoef = dlnorm(coefs[1],1,100,log=T)
  # for(i in 2:length(coefs))
  #   lpriorcoef=lpriorcoef+dnorm(log(coefs[i]),coefs[i-1],smoothSpline,log=T)
  
  lpriorcoef = dlnorm(coefs[length(coefs)],1,100,log=T)
  for(i in (length(coefs)-1):1)
    lpriorcoef=lpriorcoef+dnorm(log(coefs[i]),coefs[i+1],smoothSpline,log=T)
  
  return(priorPcoef+lpriorcoef)
}


findInitialSpline=function(xtrue,bases,knotseq,yobs,designMatrix){
  
  min=999999
  for(i in seq(.5,7,by=.1)){
    coefs=seq(.1,i,length=ncol(designMatrix))
    ytrue=coefs%*%t(designMatrix)
    if(sum(abs(yobs-ytrue))<min){ min=sum(abs(yobs-ytrue)); save=i}
  }
  coef=seq(.1,save,length=ncol(designMatrix))
  return(coef)
}


Ispline=function(intknots,lowept,upperept){
  k=3
  #determine knot sequence
  knotseq=c(rep(lowept,k+1),intknots,rep(upperept,k+1))
  numbases=length(knotseq)-k-2
  
  #create matrix of bases
  bases=matrix(NA,nrow=numbases,ncol=2)
  for(i in 1:numbases) bases[i,]=c(knotseq[i+1],knotseq[i+k+1])
  
  return(list(bases=bases,knotseq=knotseq))
}


getIsplineC=function(xtrue,knotseq,bases){
  
  numBases=nrow(bases)
  lxtrue=length(xtrue)
  
  mat=rep(0,numBases*lxtrue)
  
  storage.mode(mat) <- "double"
  storage.mode(knotseq) <- "double"
  storage.mode(xtrue) <- "double"
  storage.mode(numBases) <- "integer"
  storage.mode(lxtrue) <- "integer"
  temp=.C("getIspline",xtrue,lxtrue,knotseq,mat,numBases)
  designMatrix=matrix(temp[[4]],ncol=numBases,nrow=lxtrue)
  return(designMatrix)
  
}


updateDens=function(xtrue,state,first,niter,xgrid){
  
  nburn <- niter/2
  nsave <- niter
  nskip <- 10
  ndisplay <- 500
  mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)
  priorDP <- list(alpha=.1,m1=rep(0,1),psiinv2=solve(diag(0.5,1)),
                  nu1=4,nu2=4,tau1=1,tau2=100)
  
  if(first==TRUE){
    fit <- DPdensity(y=xtrue,prior=priorDP,mcmc=mcmc,status=TRUE,grid=xgrid)
  }else{
    fit <- DPdensity(y=xtrue,prior=priorDP,mcmc=mcmc,status=FALSE,state=state,grid=xgrid)
  }
  state=fit$state
  k=fit$state$ncluster
  densEst=fit$dens/sum(fit$dens)
  
  return(list(state=state,densEst=densEst,k=k))
}

# updateSplineXtrueDens=function(xtrue,ytrue,knotseq,bases,lowept,upperept,coefs,xcens,ycens,xsig,ysig,xobs,yobs,nobs,dens){
#   
#   pxtrue = xtrue + rnorm(nobs,0,0.5)
#   pxtrue[pxtrue<lowept] = lowept+.01
#   pxtrue[pxtrue>upperept] = upperept-.01
#   pytrue=30-2*pxtrue
#   llike = getPOPLLK(pxtrue,dens) +
#     getLLK(xobs,yobs,pxtrue,pytrue,xcens,ycens,xsig,ysig) -
#     getPOPLLK(xtrue,dens) -
#     getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig)
#   llike[is.na(llike)]=-99
#   mhrat = runif(nobs)
#   cond1=log(mhrat)<llike
#   xtrue[cond1]=pxtrue[cond1]
#   ytrue=30-2*xtrue
#   
#   return(list(xtrue=xtrue,ytrue=ytrue))
#   
# }


updateSplineXtrue=function(xtrue,ytrue,knotseq,bases,lowept,upperept,coefs,xcens,ycens,xsig,ysig,
                           xobs,yobs,nobs,dens){
  
    pxtrue = xtrue + rnorm(nobs,0,.5)
    pxtrue[pxtrue<lowept] = lowept+.01
    pxtrue[pxtrue>upperept] = upperept-.01
    designMatrix=getIsplineC(pxtrue,knotseq,bases)
    pytrue=as.numeric(coefs%*%t(designMatrix))
    llike = getPOPLLK(pxtrue,dens) +
      getLLK(xobs,yobs,pxtrue,pytrue,xcens,ycens,xsig,ysig) -
      getPOPLLK(xtrue,dens) -
      getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig)
    llike[is.na(llike)]=-99
    mhrat = runif(nobs)
    cond1=log(mhrat)<llike
    xtrue[cond1]=pxtrue[cond1]
    designMatrix=getIsplineC(xtrue,knotseq,bases)
    ytrue=coefs%*%t(designMatrix)  
  
  return(list(xtrue=xtrue,ytrue=ytrue))
  
}


updateSplineCoefs=function(smoothSpline,coefs,iter,
                 bases,xtrue,ytrue,knotseq,xobs,yobs,xcens,ycens,xsig,ysig,coefMat){
  
  if(iter<1000){
    Sigma=diag(.001,nrow=ncol(coefMat),ncol=ncol(coefMat))
  }else{
    Sigma=cov(coefMat[(iter-900):(iter-1),])
  }
  
  propCoef=tryCatch(as.numeric(exp(rmvnorm(n=1,as.numeric(log(coefs)),sigma=Sigma))),
          error=function(e){as.numeric(exp(rmvnorm(n=1,as.numeric(log(coefs)),
          sigma=diag(.001,nrow=nrow(bases),ncol=nrow(bases)))))})

  accept=0
  designMatrix=getIsplineC(xtrue,knotseq,bases)
  newY=propCoef%*%t(designMatrix)
  
  #log likelihood and priors
  oldLLK=sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))
  oldPRIORLLK=sum(getPRIORLLK(coefs,smoothSpline))
  newLLK=sum(getLLK(xobs,yobs,xtrue,newY,xcens,ycens,xsig,ysig))
  newPRIORLLK=sum(getPRIORLLK(propCoef,smoothSpline))
  
  #accept/reject
  A=newLLK+newPRIORLLK-oldLLK-oldPRIORLLK
  if(is.nan(A)==FALSE & log(runif(1)) < A){
    ytrue=newY
    coefs=propCoef
    accept=1
  }
  
  return(list(acceptCoef=accept,coefs=coefs,ytrue=ytrue))
  
}

# getCovariance=function(iter,acceptCoef,scaleSig,coefMat){
#   
#   if(iter<1000){
#     Sigma=diag(.001,nrow=ncol(coefMat),ncol=ncol(coefMat))
#   }else{
#     Sigma=cov(coefMat[(iter-900):(iter-1),])
#     if(mean(acceptCoef[(iter-500):(iter-1)])<.1) scaleSig=max(.01,scaleSig-.01)
#     if(mean(acceptCoef[(iter-500):(iter-1)])>.3) scaleSig=scaleSig+.01
#   }
#   
#   return(list(Sigma=Sigma,scaleSig=scaleSig))
#   
# }

updateSmoothParm=function(smoothSpline,coefs){
  
  accept=0
  propSmooth=rnorm(1,smoothSpline,.2)
  if(propSmooth<0 | propSmooth>2)
    return(list(smooth=smoothSpline,accept=accept))
  oldPRIORLLK=sum(getPRIORLLK(coefs,smoothSpline))
  newPRIORLLK=sum(getPRIORLLK(coefs,propSmooth))
  #accept/reject
  if(log(runif(1)) < newPRIORLLK-oldPRIORLLK){
    smoothSpline=propSmooth
    accept=1
  }
  
  return(list(smooth=smoothSpline,accept=accept))
}

# updateMeasErrorM=function(xsig,ysig,xobs,yobs,xtrue,ytrue,xcens,ycens){
#   
#   propXsig=rnorm(1,xsig,.07)
#   accept=0
#   if(propXsig<0)
#     return(list(xsig=xsig))
#   
#   oldLLK=sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))
#   newLLK=sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,propXsig,ysig))
#   
#   if(log(runif(1)) < newLLK-oldLLK)
#     xsig=propXsig
#     
#   return(list(xsig=xsig))
# }
# 
# updateMeasErrorD=function(xsig,ysig,xobs,yobs,xtrue,ytrue,xcens,ycens){
#   
#   propYsig=rnorm(1,ysig,.25)
#   if(propYsig<0)
#     return(list(ysig=ysig))
#   
#   oldLLK=sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))
#   newLLK=sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,propYsig))
#   
#   if(log(runif(1)) < newLLK-oldLLK)
#     ysig=propYsig
#   
#   return(list(ysig=ysig))
# }

findDIAC=function(yobs,gridx,weights,fit,M1,M2,xsig,ysig,minWidth,maxWidth){
  D1=0; D2=0; index=0
  lgrid=length(gridx)
  MICLowerObs=M1
  MICLowerTrue=M1-.5
  MICUpperObs=M2
  MICUpperTrue=M2-.5
  minDIA=round(quantile(yobs,probs=.1))
  maxDIA=round(quantile(yobs,probs=.9))
  storage.mode(gridx) <- "double"
  storage.mode(weights) <- "double"
  storage.mode(fit) <- "double"
  storage.mode(MICLowerObs) <- "double"
  storage.mode(MICUpperObs) <- "double"
  storage.mode(MICLowerTrue) <- "double"
  storage.mode(MICUpperTrue) <- "double"
  storage.mode(xsig) <- "double"
  storage.mode(ysig) <- "double"
  storage.mode(minDIA) <- "double"
  storage.mode(maxDIA) <- "double"
  storage.mode(lgrid) <- "integer"
  storage.mode(D1) <- "double"
  storage.mode(D2) <- "double"
  storage.mode(index) <- "double"
  storage.mode(minWidth) <- "integer"
  storage.mode(maxWidth) <- "integer"
  temp=.C("findDIATrue",gridx,weights,fit,MICLowerObs,MICUpperObs,MICLowerTrue,MICUpperTrue,
          xsig,ysig,minDIA,maxDIA,lgrid,D1,D2,index,minWidth,maxWidth)
  D1=temp[[13]]
  D2=temp[[14]]
  index=temp[[15]]
  #   print(c(D1,D2,index))
  
  return(list(D1=D1,D2=D2,index=index))
}

