
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


getPOPLLK=function(x,mu,sigma,p,k,groups)
{
  
  pxden <- log(p[groups])+dnorm(x,mu[groups],sigma[groups],log=TRUE)
  if(is.na(sum(pxden))==TRUE)  stop("Prior NA")
  return(pxden)
  
}

getPOPLLK2=function(x,dens)
{
  ###interpolate  
  pxden=log(approx(xgrid,dens,xout=x,rule=2)$y)
  return(pxden)
}

getPriorLLK=function(coefs)
{
  
  priorcoefs <- dlnorm(coefs,0,10,log=TRUE)
  return(priorcoefs)
  
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
  # groups=fit$state$ss
  # 
  # mu=rep(NA,k); sigma=rep(NA,k); p=rep(NA,k);
  # for(i in 1:k){
  #   mu[i]=mean(xtrue[groups==i])
  #   sigma[i]=sd(xtrue[groups==i])
  #   p[i]=sum(groups==i)/length(groups)
  # }
  densEst=fit$dens/sum(fit$dens)
  
  return(list(state=state,densEst=densEst,k=k))
}

getYtrue=function(coefs,xtrue){
  
  mb = (2*coefs[3]*coefs[4])/(coefs[3]+coefs[4])
  fx = 1/(1+exp(-mb*(coefs[2]-xtrue)))
  ytrue=coefs[1]*(fx*exp(coefs[3]*(coefs[2]-xtrue))+(1-fx)*exp(coefs[4]*
      (coefs[2]-xtrue)))/(1+fx*exp(coefs[3]*(coefs[2]-xtrue))+(1-fx)*
      exp(coefs[4]*(coefs[2]-xtrue)))
  
  return(ytrue)  
}

updateXtrue=function(xtrue,ytrue,coefs,xcens,ycens,xsig,ysig,xobs,yobs,nobs,
                        mu,sigma,k,p,groups,densEst){
  
  pxtrue = xtrue + rnorm(nobs,0,0.5)
  pytrue=getYtrue(coefs,pxtrue)
  llike = getPOPLLK2(pxtrue,densEst) +
    getLLK(xobs,yobs,pxtrue,pytrue,xcens,ycens,xsig,ysig) -
    getPOPLLK2(xtrue,densEst) -
    getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig)
  llike[is.na(llike)]=-99
  
  a=log(runif(nobs,0,1))
  xtrue[which(a<llike)] = pxtrue[which(a<llike)]
  ytrue[which(a<llike)] = pytrue[which(a<llike)]
  
  return(list(xtrue=xtrue,ytrue=ytrue))
  
}


### three and four together
updateLogCoef_3=function(iter,coefs,xobs,yobs,xtrue,xcens,ycens,ytrue,xsig,ysig,a1,a2,a3,a4){
  
  propCoefs=coefs
  accepts=rep(0,4)
  
  for(i in 1:3){
    if(i==1) propCoefs[i]=max(.001,coefs[i]+runif(1,-a1,a1))
    else if(i==2) propCoefs[i]=max(.001,coefs[i]+runif(1,-a2,a2))
    else if (i==3) propCoefs[i]=max(.001,coefs[i]+runif(1,-a3,a3))
    else if(i==4) propCoefs[i]=max(.001,coefs[i]+runif(1,-a4,a4))
    pytrue = getYtrue(propCoefs,xtrue)
    llike = sum(getLLK(xobs,yobs,xtrue,pytrue,xcens,ycens,xsig,ysig))+
      sum(getPriorLLK(propCoefs))-
      sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))-
      sum(getPriorLLK(coefs))
    accept=0
    if(is.na(llike)) llike=-99
    
    if(log(runif(1)) < llike){
      accept=1
      ytrue = pytrue
      coefs[i]=propCoefs[i]
      accepts[i]=1
    }else{
      propCoefs[i]=coefs[i]
    }
  }
  
  return(list(accept=accept,coefs=coefs,ytrue=ytrue,accepts=accepts))
}

updateLogCoef_2=function(iter,coefs,xobs,yobs,xtrue,xcens,ycens,ytrue,xsig,ysig,a1,a2,a3,a4){
  
  propCoefs=coefs
  accepts=rep(0,4)

  for(i in 1:4){
    if(i==1) propCoefs[i]=max(.001,coefs[i]+runif(1,-a1,a1))
    else if(i==2) propCoefs[i]=max(.001,coefs[i]+runif(1,-a2,a2))
    else if (i==3) propCoefs[i]=max(.001,coefs[i]+runif(1,-a3,a3))
    else if(i==4) propCoefs[i]=max(.001,coefs[i]+runif(1,-a4,a4))
    pytrue = getYtrue(propCoefs,xtrue)
    llike = sum(getLLK(xobs,yobs,xtrue,pytrue,xcens,ycens,xsig,ysig))+
      sum(getPriorLLK(propCoefs))-
      sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))-
      sum(getPriorLLK(coefs))
    accept=0
    if(is.na(llike)) llike=-99
    
    if(log(runif(1)) < llike){
      accept=1
      ytrue = pytrue
      coefs[i]=propCoefs[i]
      accepts[i]=1
    }else{
      propCoefs[i]=coefs[i]
    }
  }

  return(list(accept=accept,coefs=coefs,ytrue=ytrue,accepts=accepts))
}


updateLogCoef=function(iter,coefs,xobs,yobs,xtrue,xcens,ycens,ytrue,xsig,ysig,coefsMat){
  
  
  if(iter<1000){
    Sigma=diag(.001,nrow=ncol(coefsMat),ncol=ncol(coefsMat))
  }else{
    Sigma=cov(coefsMat[(iter-900):(iter-1),])
  }
  
  propCoefs=mvrnorm(n=1,mu=coefs,Sigma=Sigma)
  pytrue = getYtrue(propCoefs,xtrue)
  llike = sum(getLLK(xobs,yobs,xtrue,pytrue,xcens,ycens,xsig,ysig))+
    sum(getPriorLLK(propCoefs))-
    sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))-
    sum(getPriorLLK(coefs))
  accept=0
  if(is.na(llike)) llike=-99

  if(log(runif(1)) < llike){
    accept=1
    coefs = propCoefs
    ytrue = pytrue
  }
  
  return(list(accept=accept,coefs=coefs,ytrue=ytrue))
}


# updateMeasErrorM=function(xsig,ysig,xobs,yobs,xtrue,ytrue,xcens,ycens){
#   
#   propXsig=rnorm(1,xsig,.2)
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

