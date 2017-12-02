

updateXtrueLog=function(xtrue,ytrue,coefs,xcens,ycens,xsig,ysig,xobs,yobs,nobs,
                        mu,sigma,k,p,groups){

  pxtrue = xtrue + rnorm(nobs,0,0.5)
  pytrue=getylogtrue(coefs,pxtrue)
  llike = getPOPLLK(pxtrue,mu,sigma,p,k,groups) +
    getLLK(xobs,yobs,pxtrue,pytrue,xcens,ycens,xsig,ysig) -
    getPOPLLK(xtrue,mu,sigma,p,k,groups) -
    getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig)
  llike[is.na(llike)]=-99
  mhrat = runif(nobs)
  cond1=log(mhrat)<llike
  xtrue[cond1]=pxtrue[cond1]
  ytrue=getylogtrue(coefs,xtrue)

  return(list(xtrue=xtrue,ytrue=ytrue))

}

getylogtrue=function(coef,xtrue){

  mb = (2*coef[3]*coef[4])/(coef[3]+coef[4])
  fx = 1/(1+exp(-mb*(coef[2]-xtrue)))
  ytrue=coef[1]*(fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*exp(coef[4]*
        (coef[2]-xtrue)))/(1+fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*
        exp(coef[4]*(coef[2]-xtrue)))

  return(ytrue)
}

updateLogCoef=function(iter,coefs,xobs,yobs,xtrue,xcens,ycens,ytrue,xsig,ysig,coefMat){


  if(iter<1000){
    Sigma=diag(c(.01,.001,.001,.001),nrow=ncol(coefMat),ncol=ncol(coefMat))
  }else{
    Sigma=cov(coefMat[(iter-900):(iter-1),])
  }

  coefsa=mvrnorm(n=1,mu=coefs,Sigma=Sigma)
  pytrue = getylogtrue(coefsa,xtrue)
  llike = sum(getLLK(xobs,yobs,xtrue,pytrue,xcens,ycens,xsig,ysig)) -
    sum(getLLK(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig))
  accept=0
  if(is.na(llike)) llike=-99
  #Acceptance and make sure coeffecients are positive
  if(log(runif(1)) < llike & length(which(coefsa>0))==length(coefsa))
  {
    accept=1
    coefs = coefsa
    ytrue = pytrue
  }

  return(list(accept=accept,coefs=coefs,ytrue=ytrue))
}



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

