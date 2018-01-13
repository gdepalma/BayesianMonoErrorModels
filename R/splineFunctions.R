

getPRIORLLK_spline=function(coefs,smoothParam){

  priorPcoef=dunif(smoothParam,0,2,log=T)

  lpriorcoef = dlnorm(coefs[length(coefs)],1,100,log=T)
  for(i in (length(coefs)-1):1)
    lpriorcoef=lpriorcoef+dnorm(log(coefs[i]),coefs[i+1],smoothParam,log=T)

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
  temp=.C("getIspline",xtrue,lxtrue,knotseq,mat,numBases,PACKAGE=BayesianMonoErrorModels)
  designMatrix=matrix(temp[[4]],ncol=numBases,nrow=lxtrue)
  return(designMatrix)

}


updateSplineXtrue=function(xtrue,ytrue,mu,sigma,p,k,groups,knotseq,bases,lowept,upperept,
                           coefs,xcens,ycens,xsig,ysig,xobs,yobs,nobs){

    pxtrue = xtrue + rnorm(nobs,0,.5)
    pxtrue[pxtrue<lowept] = lowept+.01
    pxtrue[pxtrue>upperept] = upperept-.01
    designMatrix=getIsplineC(pxtrue,knotseq,bases)
    pytrue=as.numeric(coefs%*%t(designMatrix))
    llike = getPOPLLK(pxtrue,mu,sigma,p,k,groups) +
      getLLK(xobs,yobs,pxtrue,pytrue,xcens,ycens,xsig,ysig) -
      getPOPLLK(xtrue,mu,sigma,p,k,groups) -
      getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig)
    llike[is.na(llike)]=-99
    mhrat = runif(nobs)
    cond1=log(mhrat)<llike
    xtrue[cond1]=pxtrue[cond1]
    designMatrix=getIsplineC(xtrue,knotseq,bases)
    ytrue=coefs%*%t(designMatrix)

  return(list(xtrue=xtrue,ytrue=ytrue))

}


updateSplineCoefs=function(smoothParam,coefs,iter,
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
  oldPRIORLLK=sum(getPRIORLLK_spline(coefs,smoothParam))
  newLLK=sum(getLLK(xobs,yobs,xtrue,newY,xcens,ycens,xsig,ysig))
  newPRIORLLK=sum(getPRIORLLK_spline(propCoef,smoothParam))

  #accept/reject
  A=newLLK+newPRIORLLK-oldLLK-oldPRIORLLK
  if(is.nan(A)==FALSE & log(runif(1)) < A){
    ytrue=newY
    coefs=propCoef
    accept=1
  }

  return(list(acceptCoef=accept,coefs=coefs,ytrue=ytrue))

}

updateSmoothParm=function(smoothParam,coefs){

  accept=0
  propSmooth=rnorm(1,smoothParam,.2)
  if(propSmooth<0 | propSmooth>2)
    return(list(smooth=smoothParam,accept=accept))
  oldPRIORLLK=sum(getPRIORLLK_spline(coefs,smoothParam))
  newPRIORLLK=sum(getPRIORLLK_spline(coefs,propSmooth))
  #accept/reject
  if(log(runif(1)) < newPRIORLLK-oldPRIORLLK){
    smoothParam=propSmooth
    accept=1
  }

  return(list(smooth=smoothParam,accept=accept))
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



